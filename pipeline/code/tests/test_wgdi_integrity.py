from __future__ import annotations

from pathlib import Path

from pgci_preprocess.wgdi_integrity import PNG_MAGIC, check_comparison_outputs


def _valid_folder(tmp_path: Path, prefix: str = "A") -> Path:
    folder = tmp_path / "cmp"
    folder.mkdir()
    for name in (".png", ".bk.png", ".kp.png", ".pf.png"):
        (folder / f"{prefix}{name}").write_bytes(PNG_MAGIC + b"\x00\x00IEND")
    (folder / f"{prefix}.ks").write_text(
        "id1\tid2\tka_NG86\tks_NG86\tka_YN00\tks_YN00\n"
        "A1g0001\tA2g0001\t0.01\t0.1\t0.01\t0.1\n",
        encoding="utf-8",
    )
    (folder / f"{prefix}.icl").write_text(
        "# Alignment 1: score=1 pvalue=0.0 N=2 1&2 plus\nA1g0001 1 A2g0001 1 1\n",
        encoding="utf-8",
    )
    (folder / f"{prefix}.bi.csv").write_text(
        "id,chr1,chr2,start1,end1\n1,1,2,1,2\n", encoding="utf-8"
    )
    # header-only is legitimate (tandem removal can yield zero blocks)
    (folder / f"{prefix}.bi.c.csv").write_text("id,chr1,chr2\n", encoding="utf-8")
    return folder


def test_valid_outputs_pass(tmp_path: Path) -> None:
    folder = _valid_folder(tmp_path)
    assert check_comparison_outputs(folder, "A") == []


def test_truncated_png_fails(tmp_path: Path) -> None:
    folder = _valid_folder(tmp_path)
    (folder / "A.png").write_bytes(b"\x89PN")
    problems = check_comparison_outputs(folder, "A")
    assert "bad_png_magic:A.png" in problems


def test_garbage_ks_fails(tmp_path: Path) -> None:
    folder = _valid_folder(tmp_path)
    (folder / "A.ks").write_text("<html>error</html>\n", encoding="utf-8")
    assert "bad_ks_header:A.ks" in check_comparison_outputs(folder, "A", ks_enabled=True)


def test_comment_only_icl_passes(tmp_path: Path) -> None:
    folder = _valid_folder(tmp_path)
    (folder / "A.icl").write_text("# no alignments\n", encoding="utf-8")
    assert check_comparison_outputs(folder, "A") == []


def test_missing_file_reported(tmp_path: Path) -> None:
    folder = _valid_folder(tmp_path)
    (folder / "A.ks").unlink()
    assert "missing_or_empty:A.ks" in check_comparison_outputs(folder, "A", ks_enabled=True)


def test_core_mode_ignores_ks_products(tmp_path: Path) -> None:
    folder = _valid_folder(tmp_path)
    for name in ("A.ks", "A.bi.csv", "A.bi.c.csv", "A.bk.png", "A.kp.png", "A.pf.png"):
        (folder / name).unlink()
    assert check_comparison_outputs(folder, "A") == []
    problems = check_comparison_outputs(folder, "A", ks_enabled=True)
    assert "missing_or_empty:A.ks" in problems
    assert "missing_or_empty:A.bk.png" in problems