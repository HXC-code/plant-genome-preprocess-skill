"""Launch pinned upstream WGDI with the compatibility patches managed here."""
import sys
import os
import shutil
from pathlib import Path

# 兼容 matplotlib >= 3.9：wgdi 0.75 的 block_ks.py 用了已移除的 cm.get_cmap
try:
    import matplotlib
    import matplotlib.cm as _cm
    if not hasattr(_cm, "get_cmap"):
        if hasattr(matplotlib, "colormaps"):
            _cm.get_cmap = matplotlib.colormaps.get_cmap
        else:
            from matplotlib import cm as _cm2
            _cm.get_cmap = getattr(_cm2, "get_cmap", None)
except Exception:
    pass

import wgdi.base as _wgdi_base

_project_root = Path(__file__).resolve().parent.parent
_original_config = _wgdi_base.config


def _resolve_tool(environment_name: str, legacy_environment: str, command: str):
    configured = os.environ.get(environment_name) or os.environ.get(legacy_environment)
    return configured or shutil.which(command)


def _config_with_portable_tools():
    items = dict(_original_config())
    overrides = {
        "mafft_path": _resolve_tool("WGDI_MAFFT", "MAFFT", "mafft"),
        "muscle_path": _resolve_tool("WGDI_MUSCLE", "MUSCLE", "muscle"),
        "yn00_path": _resolve_tool("WGDI_YN00", "YN00", "yn00"),
        "iqtree_path": _resolve_tool("WGDI_IQTREE", "IQTREE", "iqtree"),
        "trimal_path": _resolve_tool("WGDI_TRIMAL", "TRIMAL", "trimal"),
        "fasttree_path": _resolve_tool("WGDI_FASTTREE", "FASTTREE", "fasttree"),
        "divvier_path": _resolve_tool("WGDI_DIVVIER", "DIVVIER", "divvier"),
    }
    pal2nal = (
        os.environ.get("WGDI_PAL2NAL")
        or os.environ.get("PAL2NAL")
        or _project_root / "code" / "tools" / "pal2nal.v14" / "pal2nal.pl"
    )
    if Path(pal2nal).is_file():
        overrides["pal2nal_path"] = str(pal2nal)
    for key, value in overrides.items():
        if value:
            items[key] = value
    return list(items.items())


_wgdi_base.config = _config_with_portable_tools

_original_str_to_bool = _wgdi_base.str_to_bool


def _str_to_bool(value):
    if isinstance(value, bool):
        return value
    return _original_str_to_bool(value)


_wgdi_base.str_to_bool = _str_to_bool

# Workflow extension: ``-r`` as the final argument means "remove tandem
# duplicates during tandem-aware WGDI stages".  Upstream ``-r <config>`` remains
# reserved for its original retention analysis.  ``--remove-tandem`` is the
# unambiguous long form and is accepted anywhere.
_raw_args = list(sys.argv[1:])
_remove_tandem = "--remove-tandem" in _raw_args or (
    bool(_raw_args) and _raw_args[-1] == "-r"
)
if _remove_tandem:
    sys.argv = [
        item
        for item in sys.argv
        if item not in {"--remove-tandem"}
    ]
    if sys.argv[-1:] == ["-r"]:
        sys.argv.pop()

    _original_load_conf = _wgdi_base.load_conf

    def _load_conf_remove_tandem(file, section):
        options = _original_load_conf(file, section)
        if section in {"correspondence", "blockks", "kspeaks"}:
            options = [
                (key, "false" if key == "tandem" else value)
                for key, value in options
            ]
            if not any(key == "tandem" for key, _value in options):
                options.append(("tandem", "false"))
        return options

    _wgdi_base.load_conf = _load_conf_remove_tandem

from wgdi.run import main

if __name__ == "__main__":
    sys.argv[0] = "wgdi"
    main()
