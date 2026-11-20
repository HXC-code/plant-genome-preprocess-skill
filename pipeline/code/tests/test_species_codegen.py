"""Unit tests for the 3-letter species code rule (genus initial + 2 epithet letters)."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from pgci_preprocess.species_codegen import derive_code, derive_codes, split_latin_name
from pgci_preprocess.species_names import resolve_species_code

PROJECT_ROOT = Path(__file__).resolve().parents[1]
TABLE = PROJECT_ROOT / "resources" / "species_abbreviations.tsv"


def _table_codes() -> set[str]:
    with TABLE.open(encoding="utf-8-sig", newline="") as handle:
        return {
            (row["abbreviation"] or "").strip().casefold()
            for row in csv.DictReader(handle, delimiter="\t")
            if (row["abbreviation"] or "").strip()
        }


def test_split_latin_name_drops_authority() -> None:
    assert split_latin_name("Cercis chinensis Bunge") == ("Cercis", "chinensis")
    assert split_latin_name("Pinus taeda L.") == ("Pinus", "taeda")


def test_default_code_is_genus_initial_plus_first_two_epithet_letters() -> None:
    assert derive_code("Cercis chinensis") == "Cch"
    assert derive_code("Pinus taeda") == "Pta"


def test_collision_moves_to_another_two_letter_combination() -> None:
    assert derive_code("Cercis chinensis", {"Cch"}) == "Cci"
    assert derive_code("Cercis chinensis", {"Cch", "Cci"}) == "Ccn"


def test_collision_widens_to_four_letters_when_three_is_exhausted() -> None:
    epithet = "taeda"
    pairs = {f"P{epithet[first]}{epithet[second]}" for first in range(5) for second in range(5) if first < second}
    code = derive_code("Pinus taeda", pairs)
    assert code == "Ptae"
    assert len(code) == 4


def test_short_epithet_falls_back_to_numeric_suffix() -> None:
    code = derive_code("Pinus ta", {"Pta"})
    assert code == "Pta2"


def test_sequence_assignment_never_duplicates() -> None:
    codes = derive_codes(["Cercis chinensis", "Cercis chinensis", "Cinnamomum chinensis"])
    assert len({code.casefold() for code in codes}) == 3
    assert codes[0] == "Cch"


def test_rule_reproduces_existing_table_code_for_known_species() -> None:
    """The published table was produced by this rule, so re-deriving must match."""

    taken = _table_codes() - {"cci"}
    assert derive_code("Cercis chinensis", taken) == "Cci"


def test_table_has_no_duplicate_codes() -> None:
    codes = [code for code in _table_codes()]
    assert len(codes) == len(set(codes))


def test_new_species_code_avoids_every_table_code() -> None:
    code = derive_code("Cercidopsis novissima", _table_codes())
    assert code.casefold() not in _table_codes()


def test_resolve_species_code_derives_instead_of_falling_back_to_s_id(tmp_path: Path) -> None:
    (tmp_path / "resources").mkdir()
    (tmp_path / "resources" / "species_abbreviations.tsv").write_bytes(TABLE.read_bytes())
    code = resolve_species_code(tmp_path, "999999", latin_name="Testus septimus")
    assert code == derive_code("Testus septimus", _table_codes())
    assert code != "s999999"


def test_known_species_still_uses_the_table_row(tmp_path: Path) -> None:
    (tmp_path / "resources").mkdir()
    (tmp_path / "resources" / "species_abbreviations.tsv").write_bytes(TABLE.read_bytes())
    assert resolve_species_code(tmp_path, "1126", latin_name="Cercis chinensis") == "Cci"


@pytest.mark.parametrize("name", ["", "Cercis", "123"])
def test_invalid_latin_name_is_rejected(name: str) -> None:
    with pytest.raises(ValueError):
        derive_code(name)
