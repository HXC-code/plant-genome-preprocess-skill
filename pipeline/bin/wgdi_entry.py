"""WGDI 启动入口：先注入 matplotlib 兼容补丁，再运行 wgdi.run。

规避 sitecustomize/PYTHONPATH 在 Git Bash + Windows Python 混合环境下的路径问题：
补丁在进程内直接执行，100% 生效。不影响 Ks 等非绘图流程。
"""
import sys

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
