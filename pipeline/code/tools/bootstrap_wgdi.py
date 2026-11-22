#!/usr/bin/env python3
"""Create the pinned WGDI virtual environment — but scan the machine first.

Flow: preflight (read-only) -> decide reuse/install/offline/blocked -> act.
Nothing is installed before the scan reports what is already present, what
conflicts, and whether the package index is reachable at all. Failures write a
log and print a classified diagnosis instead of raw pip noise.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import venv
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pgci_preprocess import envcheck  # noqa: E402
from pgci_preprocess.version import __version__  # noqa: E402

EXTERNAL_TOOLS = envcheck.EXTERNAL_TOOLS


def venv_python(env_dir: Path) -> Path:
    for candidate in (env_dir / "Scripts" / "python.exe", env_dir / "bin" / "python"):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"虚拟环境创建失败：{env_dir}")


def classify_failure(log_text: str) -> str:
    low = log_text.lower()
    if any(k in low for k in ("could not fetch", "connection", "timed out", "timeout",
                              "temporary failure", "name or service not known",
                              "proxy", "ssl", "certificate")):
        return "network"
    if any(k in low for k in ("no matching distribution", "requires-python",
                              "could not find a version")):
        return "wheel"
    if any(k in low for k in ("resolutionimpossible", "conflict", "incompatible")):
        return "conflict"
    return "unknown"


ADVICE = {
    "network": "索引不可达（网络/代理/防火墙）。不要反复重试：配置镜像后重试\n"
               "  bootstrap_wgdi.py --index-url https://pypi.tuna.tsinghua.edu.cn/simple\n"
               "或离线运行预处理（不需要 pip）：\n"
               "  set PYTHONPATH=<repo>\\code\\src && python -m pgci_preprocess run --project-root <项目>",
    "wheel": "钉版没有适配当前 Python 的 wheel。换用 Python 3.11–3.13，"
             "或用 --python-exe 指定一个受支持的解释器。",
    "conflict": "依赖解析冲突。查看日志中的冲突包；默认应装进独立 .wgdi-env，"
                "不要装进 conda base 或已有冲突的环境。",
    "unknown": "查看完整日志定位原因；把日志尾部和 preflight --json 输出一起反馈。",
}


def run_logged(argv: list[str], log_path: Path) -> None:
    print("+", " ".join(argv), flush=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write("+ " + " ".join(argv) + "\n")
        result = subprocess.run(argv, capture_output=True, text=True,
                                encoding="utf-8", errors="replace", check=False)
        log.write(result.stdout)
        log.write(result.stderr)
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-25:])
        kind = classify_failure(result.stdout + result.stderr)
        raise RuntimeError(
            f"安装失败（{kind}）。日志: {log_path}\n--- 日志尾部 ---\n{tail}\n"
            f"--- 建议 ---\n{ADVICE[kind]}"
        )


def bootstrap(args: argparse.Namespace) -> int:
    root = args.project_root.resolve()
    requirements = (args.requirements or root / "requirements-wgdi-validated.txt").resolve()
    python_exe = str(Path(args.python_exe).resolve()) if args.python_exe else sys.executable

    report = envcheck.scan(
        root,
        python_exe=python_exe,
        requirements=requirements,
        source_version=__version__,
        probe_index=not args.offline,
        index_url=args.index_url,
    )
    print(envcheck.render(report))
    print(flush=True)
    plan = str(report["plan"])

    if args.check_only:
        return {"reuse": 0, "install": 0, "offline": 3, "blocked": 2}.get(plan, 2)
    if plan == "blocked":
        print("错误：预扫发现阻断项，未安装任何东西。先解决上述 [blocker] 再重试。",
              file=sys.stderr)
        return 2

    env_dir = (args.env_dir or root / ".wgdi-env").resolve()
    log_path = env_dir / "bootstrap.log"

    if args.offline:
        print("[OFFLINE] 跳过全部 pip 安装。预处理可用源码模式运行：")
        setenv = "set" if os.name == "nt" else "export"
        print(f"  {setenv} PYTHONPATH={root / 'code' / 'src'}")
        print("  python -m pgci_preprocess run --project-root <项目>")
        print("或离线安装核心包（核心无第三方依赖，用本地 setuptools，不需网络）：")
        print(f"  {python_exe} -m pip install -e {root / 'code'} --no-build-isolation --no-deps")
        print("WGDI 下游运行时需要依赖就绪后再执行 bootstrap（或配置镜像）。")
        return 0

    if args.core_only:
        # Offline-safe core install: no build isolation (local setuptools), no
        # deps (core is pure stdlib). Works on locked-down machines where the
        # normal `pip install -e .[gbff]` fails at the first network call.
        code_dir = root / "code"
        run_logged(
            [python_exe, "-m", "pip", "install", "-e", str(code_dir),
             "--no-build-isolation", "--no-deps"],
            root / "_pipeline_logs" / "core-install.log",
        )
        check = subprocess.run(
            [python_exe, "-c", "import pgci_preprocess; print(pgci_preprocess.__version__)"],
            capture_output=True, text=True, encoding="utf-8", errors="replace", check=False,
        )
        if check.returncode != 0:
            raise RuntimeError(f"核心包安装后无法导入：{check.stdout}{check.stderr}")
        print(f"[CORE] 预处理核心已离线安装: {check.stdout.strip()}")
        print("GBFF 转换如需 biopython，另行离线安装 wheel；WGDI 下游另跑 bootstrap。")
        return 0

    if plan == "reuse" and not args.force:
        runtime_python = venv_python(env_dir)
        print(f"[REUSE] 已有环境满足钉版，跳过安装: {runtime_python}")
    else:
        if env_dir.exists() and args.force:
            shutil.rmtree(env_dir)
        if not env_dir.exists():
            print(f"[CREATE] 独立虚拟环境: {env_dir}")
            venv.EnvBuilder(with_pip=True).create(env_dir)
        runtime_python = venv_python(env_dir)
        if args.upgrade_pip:
            run_logged([str(runtime_python), "-m", "pip", "install", "--upgrade", "pip"], log_path)
        if not requirements.is_file():
            raise FileNotFoundError(requirements)
        pip_cmd = [str(runtime_python), "-m", "pip", "install",
                   "--retries", "2", "--timeout", "30", "-r", str(requirements)]
        if args.index_url:
            pip_cmd += ["--index-url", args.index_url]
        try:
            run_logged(pip_cmd, log_path)
        except RuntimeError:
            if not args.allow_unpinned_fallback:
                raise
            print("[WARN] 钉版安装失败，按显式请求改用未锁定依赖（pandas 仍限制 <3）。",
                  file=sys.stderr)
            fallback = [str(runtime_python), "-m", "pip", "install",
                        "--retries", "2", "--timeout", "30",
                        "wgdi", "matplotlib", "numpy", "scipy", "pandas<3",
                        "biopython", "tabulate"]
            if args.index_url:
                fallback += ["--index-url", args.index_url]
            run_logged(fallback, log_path)

    check = subprocess.run(
        [str(runtime_python), "-c",
         "import wgdi, wgdi.run, pandas; "
         "assert int(pandas.__version__.split('.')[0]) < 3, 'pandas>=3 breaks WGDI 0.75'; "
         "print('wgdi import OK, pandas', pandas.__version__)"],
        capture_output=True, text=True, encoding="utf-8", errors="replace", check=False,
    )
    if check.returncode != 0:
        raise RuntimeError(f"运行时校验失败：{check.stdout}{check.stderr}")
    print(check.stdout.strip())

    pal2nal = root / "external" / "pal2nal.v14" / "pal2nal.pl"
    if args.skip_pal2nal:
        print("[SKIP] PAL2NAL installation")
    elif pal2nal.is_file():
        print(f"PAL2NAL 已安装: {pal2nal}")
    else:
        from install_pal2nal import DEFAULT_SHA256, DEFAULT_URL, install as install_pal2nal

        installed = install_pal2nal(DEFAULT_URL, DEFAULT_SHA256, pal2nal.parent, force=False)
        print(f"PAL2NAL 已安装: {installed / 'pal2nal.pl'}")
    print(f"\nPython 运行时已就绪: {runtime_python}")
    print(f"启动方式: {root / 'bin' / 'wgdi'} <wgdi 参数>")

    if not args.skip_external_check:
        print("\n外部程序需要单独安装并可从 PATH 解析：")
        for tool in EXTERNAL_TOOLS:
            path = shutil.which(tool)
            print(f"  [{'OK' if path else 'MISS'}] {tool}" + (f" -> {path}" if path else ""))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--python-exe", default="")
    parser.add_argument("--env-dir", type=Path)
    parser.add_argument("--requirements", type=Path)
    parser.add_argument("--index-url", default="", help="PyPI mirror for restricted networks.")
    parser.add_argument("--offline", action="store_true",
                        help="Skip pip entirely; report source-mode usage for preprocessing.")
    parser.add_argument("--core-only", action="store_true",
                        help="Offline-safe editable install of the pure-stdlib core only "
                             "(--no-build-isolation --no-deps); no WGDI venv.")
    parser.add_argument("--check-only", action="store_true",
                        help="Only run the preflight scan and exit with its verdict.")
    parser.add_argument("--upgrade-pip", action="store_true",
                        help="Opt-in pip self-upgrade (needs network; off by default).")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-external-check", action="store_true")
    parser.add_argument("--skip-pal2nal", action="store_true")
    parser.add_argument("--allow-unpinned-fallback", action="store_true")
    args = parser.parse_args()
    try:
        return bootstrap(args)
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"错误：{exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
