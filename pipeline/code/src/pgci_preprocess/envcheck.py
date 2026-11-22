"""Environment preflight: scan the machine BEFORE installing anything.

The historical bootstrap created a virtualenv and ran pip blind. On a fresh or
locked-down machine that fails at the first network call and prints raw pip
noise, with no answer to "what is already here, what conflicts, what should I
do instead". This module answers those questions first:

* which interpreter is this (version, conda base vs venv, writable),
* which relevant distributions are already installed and at what version,
* which of them conflict with the pinned WGDI runtime (pandas>=3 breaks
  WGDI 0.75 string handling; a foreign wgdi version shadows the pin),
* is the package index reachable at all (and is a mirror configured),
* are the external bioinformatics tools resolvable,
* is there an existing environment that can simply be reused.

The scan is read-only. ``choose_plan`` turns the report into one of:
``reuse`` / ``install`` / ``offline`` / ``blocked``.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

MIN_PYTHON = (3, 11)
RELEVANT = (
    "pgci-preprocess",
    "wgdi",
    "pandas",
    "numpy",
    "scipy",
    "matplotlib",
    "biopython",
    "tabulate",
    "pip",
    "setuptools",
)
EXTERNAL_TOOLS = ("diamond", "blastp", "makeblastdb", "muscle", "yn00", "perl", "bash")
DEFAULT_INDEX = "https://pypi.org/simple/"


def _probe(python_exe: str, code: str, timeout: int = 30) -> str:
    try:
        result = subprocess.run(
            [python_exe, "-c", code],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    return result.stdout.strip() if result.returncode == 0 else ""


def interpreter_info(python_exe: str) -> dict[str, Any]:
    code = (
        "import sys,os,json;print(json.dumps({"
        "'version':'.'.join(map(str,sys.version_info[:3])),"
        "'major':sys.version_info[0],'minor':sys.version_info[1],"
        "'executable':sys.executable,'prefix':sys.prefix,"
        "'base_prefix':getattr(sys,'base_prefix',sys.prefix),"
        "'conda':bool(os.environ.get('CONDA_PREFIX')) or 'conda' in sys.prefix.lower()"
        "}))"
    )
    raw = _probe(python_exe, code)
    if not raw:
        return {"version": "", "major": 0, "minor": 0, "executable": python_exe,
                "prefix": "", "base_prefix": "", "conda": False, "reachable": False}
    info = json.loads(raw)
    info["reachable"] = True
    info["is_venv"] = info["prefix"] != info["base_prefix"]
    info["is_conda_base"] = bool(info["conda"]) and not info["is_venv"]
    return info


def installed_packages(python_exe: str) -> dict[str, str]:
    code = (
        "import importlib.metadata as m,json;"
        "print(json.dumps({(d.metadata['Name'] or '').lower():d.version "
        "for d in m.distributions()}))"
    )
    raw = _probe(python_exe, code)
    if not raw:
        return {}
    try:
        return {k: v for k, v in json.loads(raw).items() if k in RELEVANT}
    except ValueError:
        return {}


def parse_requirements(path: Path | str) -> dict[str, str]:
    pins: dict[str, str] = {}
    path = Path(path)
    if not path.is_file():
        return pins
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        match = re.match(r"^([A-Za-z0-9_.-]+)==([^\s;]+)", line)
        if match:
            pins[match.group(1).lower()] = match.group(2)
    return pins


def _major(version: str) -> int:
    match = re.match(r"^(\d+)", version or "")
    return int(match.group(1)) if match else 0


def evaluate_conflicts(
    installed: dict[str, str],
    pins: dict[str, str],
    source_version: str = "",
) -> list[dict[str, str]]:
    """Pure conflict evaluation; unit-testable without touching the machine."""

    problems: list[dict[str, str]] = []
    pandas = installed.get("pandas", "")
    if pandas and _major(pandas) >= 3:
        problems.append({
            "level": "blocker",
            "item": "pandas",
            "detail": f"pandas {pandas}：WGDI 0.75 的字符串赋值在 pandas>=3 下报错；"
                      f"目标环境必须用 pandas {pins.get('pandas', '2.x')}",
        })
    wgdi = installed.get("wgdi", "")
    if wgdi and pins.get("wgdi") and wgdi != pins["wgdi"]:
        problems.append({
            "level": "warn",
            "item": "wgdi",
            "detail": f"已装 wgdi {wgdi}，钉版为 {pins['wgdi']}；"
                      "在独立 venv 中安装会覆盖，复用 base 环境则会冲突",
        })
    own = installed.get("pgci-preprocess", "")
    if own and source_version and own != source_version:
        problems.append({
            "level": "warn",
            "item": "pgci-preprocess",
            "detail": f"site-packages 里已有 pgci-preprocess {own}，源码为 {source_version}；"
                      "可能遮蔽源码树（用 PYTHONPATH=src 或卸载旧包）",
        })
    for name, pin in pins.items():
        have = installed.get(name, "")
        if have and have != pin and name not in {"wgdi", "pandas"}:
            problems.append({
                "level": "info",
                "item": name,
                "detail": f"已装 {have}，钉版 {pin}（独立 venv 中会按钉版安装）",
            })
    return problems


def probe_network(index_url: str = "", timeout: int = 6) -> dict[str, Any]:
    url = index_url or os.environ.get("PIP_INDEX_URL") or DEFAULT_INDEX
    try:
        request = urllib.request.Request(url, method="GET")
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return {"reachable": True, "index": url, "status": response.status}
    except (urllib.error.URLError, urllib.error.HTTPError, OSError, ValueError) as exc:
        return {"reachable": False, "index": url, "error": str(exc)}


def external_tools() -> dict[str, str]:
    return {tool: (shutil.which(tool) or "") for tool in EXTERNAL_TOOLS}


def existing_env(root: Path) -> dict[str, Any]:
    env_dir = root / ".wgdi-env"
    for rel in ("Scripts/python.exe", "bin/python"):
        candidate = env_dir / rel
        if candidate.is_file():
            info = interpreter_info(str(candidate))
            packages = installed_packages(str(candidate))
            return {
                "present": True,
                "python": str(candidate),
                "version": info.get("version", ""),
                "wgdi": packages.get("wgdi", ""),
                "pandas": packages.get("pandas", ""),
            }
    return {"present": False}


def choose_plan(report: dict[str, Any]) -> str:
    """reuse / install / offline / blocked — decided from the scan, not guesswork."""

    info = report["interpreter"]
    if not info.get("reachable") or (info.get("major", 0), info.get("minor", 0)) < MIN_PYTHON:
        return "blocked"
    blockers = [p for p in report["conflicts"] if p["level"] == "blocker"]
    env = report["existing_env"]
    pins = report["pins"]
    if env.get("present") and env.get("wgdi") == pins.get("wgdi") and not blockers:
        pandas = str(env.get("pandas", ""))
        if not pandas or _major(pandas) < 3:
            return "reuse"
    if blockers:
        return "blocked"
    if not report["network"].get("reachable"):
        # An unprobed index is treated like an unreachable one: never start
        # pip when we have no evidence that a package index can be reached.
        return "offline"
    return "install"


def scan(
    project_root: Path | str,
    python_exe: str = "",
    requirements: Path | str = "",
    source_version: str = "",
    probe_index: bool = True,
    index_url: str = "",
) -> dict[str, Any]:
    root = Path(project_root).resolve()
    python_exe = python_exe or sys.executable
    req_path = Path(requirements) if requirements else root / "requirements-wgdi-validated.txt"
    pins = parse_requirements(req_path)
    info = interpreter_info(python_exe)
    installed = installed_packages(python_exe) if info.get("reachable") else {}
    conflicts = evaluate_conflicts(installed, pins, source_version)
    if info.get("is_conda_base"):
        conflicts.append({
            "level": "warn",
            "item": "conda-base",
            "detail": "当前解释器是 conda base 环境；默认会另建 .wgdi-env，"
                      "不要把钉版依赖装进 base",
        })
    network = (
        probe_network(index_url)
        if probe_index
        else {"reachable": False, "skipped": True, "index": "(skipped)"}
    )
    report: dict[str, Any] = {
        "project_root": str(root),
        "interpreter": info,
        "installed": installed,
        "pins": pins,
        "conflicts": conflicts,
        "network": network,
        "external_tools": external_tools(),
        "existing_env": existing_env(root),
        "writable": os.access(root, os.W_OK),
        "path_length": len(str(root)),
    }
    report["plan"] = choose_plan(report)
    return report


def render(report: dict[str, Any]) -> str:
    info = report["interpreter"]
    lines = [
        "== 环境预扫（只读，未安装任何东西）==",
        f"解释器   : {info.get('executable')}  Python {info.get('version') or '不可用'}"
        + ("  [conda base]" if info.get("is_conda_base") else "")
        + ("  [venv]" if info.get("is_venv") else ""),
        "已装相关包: " + (", ".join(f"{k}={v}" for k, v in sorted(report["installed"].items())) or "(无)"),
        "网络/索引 : "
        + (
            f"{report['network'].get('index')}  未探测（--no-network，无法判断安装可行性）"
            if report["network"].get("skipped")
            else f"{report['network'].get('index')}  "
            + ("可达" if report["network"].get("reachable")
               else f"不可达 ({report['network'].get('error', '')})")
        ),
        "已有环境 : " + (
            f".wgdi-env  Python {report['existing_env'].get('version')}  wgdi={report['existing_env'].get('wgdi') or '无'}"
            if report["existing_env"].get("present") else "无 .wgdi-env"
        ),
        "外部工具 : " + ", ".join(
            f"{tool}{'✓' if path else '✗'}" for tool, path in report["external_tools"].items()
        ),
        f"目录可写 : {'是' if report['writable'] else '否'}   路径长度: {report['path_length']}",
    ]
    if report["conflicts"]:
        lines.append("冲突/提示 :")
        for problem in report["conflicts"]:
            lines.append(f"  [{problem['level']}] {problem['item']}: {problem['detail']}")
    plan = report["plan"]
    advice = {
        "reuse": "结论: 已有 .wgdi-env 满足钉版 → 直接复用，无需安装。",
        "install": "结论: 环境干净且索引可达 → 创建/使用 .wgdi-env 并按钉版安装。",
        "offline": (
            "结论: 索引未探测（--no-network）→ 安装可行性未知；要安装请去掉 --no-network 重扫。"
            if report["network"].get("skipped")
            else "结论: 索引不可达 → 不要反复重试 pip。改用离线源码模式运行预处理"
                 "（PYTHONPATH=code/src），或配置镜像后重试："
                 "bootstrap_wgdi.py --index-url https://pypi.tuna.tsinghua.edu.cn/simple"
        ),
        "blocked": "结论: 存在阻断项（Python 版本/解释器不可用/pandas>=3 冲突）→ 先解决再安装，"
                   "不要往当前环境里装。",
    }
    lines.append(advice.get(plan, f"结论: {plan}"))
    return "\n".join(lines)
