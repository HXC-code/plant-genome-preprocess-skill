#!/usr/bin/env python3
"""Create the pinned WGDI virtual environment used by the workflow launcher."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import venv
from pathlib import Path


EXTERNAL_TOOLS = ("diamond", "blastp", "makeblastdb", "muscle", "yn00", "perl")


def run(argv: list[str]) -> None:
    print("+", " ".join(argv), flush=True)
    subprocess.run(argv, check=True)


def venv_python(env_dir: Path) -> Path:
    candidates = (env_dir / "Scripts" / "python.exe", env_dir / "bin" / "python")
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"虚拟环境创建失败：{env_dir}")


def bootstrap(args: argparse.Namespace) -> int:
    root = args.project_root.resolve()
    env_dir = (args.env_dir or root / ".wgdi-env").resolve()
    requirements = (args.requirements or root / "requirements-wgdi-validated.txt").resolve()
    python_exe = Path(args.python_exe).resolve() if args.python_exe else Path(sys.executable)
    if not python_exe.is_file():
        raise FileNotFoundError(f"找不到 Python：{python_exe}")
    version = subprocess.check_output(
        [str(python_exe), "-c", "import sys; print('.'.join(map(str, sys.version_info[:2])))"],
        text=True,
    ).strip()
    if tuple(map(int, version.split("."))) < (3, 11):
        raise RuntimeError(f"需要 Python >= 3.11，{python_exe} 为 {version}")

    if env_dir.exists() and args.force:
        shutil.rmtree(env_dir)
    if not env_dir.exists():
        venv.EnvBuilder(with_pip=True).create(env_dir)
    runtime_python = venv_python(env_dir)

    run([str(runtime_python), "-m", "pip", "install", "--upgrade", "pip"])
    try:
        if not requirements.is_file():
            raise FileNotFoundError(requirements)
        run([str(runtime_python), "-m", "pip", "install", "-r", str(requirements)])
    except (FileNotFoundError, subprocess.CalledProcessError):
        if not args.allow_unpinned_fallback:
            raise RuntimeError(
                "固定版本 WGDI 依赖安装失败。请修复 requirements，或显式使用 "
                "--allow-unpinned-fallback。"
            )
        print("[WARN] 固定版本安装失败，按显式请求改用未锁定依赖。", file=sys.stderr)
        run(
            [
                str(runtime_python),
                "-m",
                "pip",
                "install",
                "wgdi",
                "matplotlib",
                "numpy",
                "scipy",
                "pandas",
                "biopython",
                "tabulate",
            ]
        )

    run([str(runtime_python), "-c", "import wgdi, wgdi.run; print('wgdi import OK')"])
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
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-external-check", action="store_true")
    parser.add_argument("--allow-unpinned-fallback", action="store_true")
    args = parser.parse_args()
    try:
        return bootstrap(args)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f"错误：{exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
