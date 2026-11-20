from __future__ import annotations

from pathlib import Path


WORKFLOW_ROOT = Path(__file__).resolve().parents[2]


def test_generated_wgdi_configs_remove_tandem_duplicates_by_default() -> None:
    for name in ("4run_wgdi_different.sh", "5run_wgdi_same.sh"):
        text = (WORKFLOW_ROOT / name).read_text(encoding="utf-8")
        assert 'wgdi_tandem="${WGDI_TANDEM:-false}"' in text
        assert "WGDI_TANDEM_ARGS=(-r)" in text
        assert text.count("tandem = ${wgdi_tandem}") == 3
        assert '"\\${WGDI_TANDEM_ARGS[@]}" -c total.conf' in text
        assert '"\\${WGDI_TANDEM_ARGS[@]}" -bk total.conf' in text
        assert '"\\${WGDI_TANDEM_ARGS[@]}" -kp total.conf' in text


def test_wgdi_wrapper_treats_final_r_as_tandem_removal_extension() -> None:
    text = (WORKFLOW_ROOT / "bin" / "wgdi_entry.py").read_text(encoding="utf-8")
    assert '_raw_args[-1] == "-r"' in text
    assert '"--remove-tandem"' in text
    assert '("tandem", "false")' in text


def test_blockks_parses_tandem_as_boolean() -> None:
    text = (WORKFLOW_ROOT / "vendor" / "wgdi" / "block_ks.py").read_text(
        encoding="utf-8"
    )
    assert "self.tandem = base.str_to_bool(self.tandem)" in text
    assert "self.remove_tandem =" not in text
