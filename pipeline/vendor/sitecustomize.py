"""WGDI 兼容性补丁：适配 matplotlib >= 3.9（cm.get_cmap 已移除）。

wgdi 0.75 的 block_ks.py 使用 plt.cm.get_cmap('gist_rainbow')，
在 matplotlib 3.11 中会抛 AttributeError。这里在 Python 启动时注入一个
向后兼容的 shim，把 matplotlib.cm.get_cmap 映射到 matplotlib.colormaps.get_cmap。
其余行为不变。
"""
try:
    import matplotlib
    import matplotlib.cm as _cm
    if not hasattr(_cm, "get_cmap"):
        if hasattr(matplotlib, "colormaps"):
            _cm.get_cmap = matplotlib.colormaps.get_cmap
        else:  # 更老版本兜底
            from matplotlib import cm as _cm2
            _cm.get_cmap = getattr(_cm2, "get_cmap", None)
except Exception:
    # 补丁失败不应阻断任何非绘图流程（如 Ks 计算）
    pass
