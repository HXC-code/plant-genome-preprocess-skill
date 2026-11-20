from __future__ import annotations

from dataclasses import dataclass

DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN")
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBXZUOJ*")
BASES = "TCAG"
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STANDARD_CODON_TABLE = {
    first + second + third: AMINO_ACIDS[index]
    for index, (first, second, third) in enumerate(
        (a, b, c) for a in BASES for b in BASES for c in BASES
    )
}
CODON_TABLE = STANDARD_CODON_TABLE

# NCBI translation-table differences from table 1.  These cover the tables
# commonly attached to plant nuclear, plastid and mitochondrial GFF/GBFF data.
# Translation is deliberately independent of Biopython so GFF+genome extraction
# remains usable in the standard-library-only preprocessing environment.
TABLE_OVERRIDES: dict[int, dict[str, str]] = {
    1: {},
    2: {"TGA": "W", "ATA": "M", "AGA": "*", "AGG": "*"},
    3: {
        "ATA": "M",
        "CTT": "T",
        "CTC": "T",
        "CTA": "T",
        "CTG": "T",
        "TGA": "W",
    },
    4: {"TGA": "W"},
    5: {"TGA": "W", "ATA": "M"},
    6: {"TAA": "Q", "TAG": "Q"},
    9: {"TGA": "W", "AAA": "N", "AGA": "S", "AGG": "S"},
    10: {"TGA": "C"},
    11: {"TGA": "W"},
    12: {"CTG": "S", "TGA": "W"},
    13: {"TGA": "W", "AGA": "G", "AGG": "G", "ATA": "M"},
    14: {
        "TAA": "Y",
        "TAG": "Y",
        "TGA": "W",
        "AAA": "N",
        "AGA": "S",
        "AGG": "S",
        "ATA": "M",
    },
    16: {"TGA": "W", "TAG": "L"},
    21: {"TGA": "W", "ATA": "M", "AAA": "N", "AGA": "S", "AGG": "S"},
    22: {"TCA": "*", "TAG": "L"},
    23: {"TTA": "*", "TGA": "W", "AGA": "S", "AGG": "S"},
    24: {"TGA": "W", "AGA": "S", "AGG": "K"},
    25: {"TGA": "G"},
    26: {"TGA": "G"},
    27: {"TAA": "Q", "TAG": "Q", "TGA": "W"},
    28: {"TAA": "Q", "TAG": "Q", "TGA": "W"},
    29: {"TAA": "Y", "TAG": "Y"},
    30: {"TAA": "E", "TAG": "E"},
    31: {"TAA": "E", "TAG": "E"},
    33: {"TAG": "L", "TGA": "W"},
}
SUPPORTED_TRANSLATION_TABLES = frozenset(TABLE_OVERRIDES)


@dataclass(frozen=True, slots=True)
class TranslationCheck:
    valid: bool
    reason: str
    exact: bool
    mismatches: int
    translated: str


def get_codon_table(translation_table: int = 1) -> dict[str, str]:
    if translation_table not in TABLE_OVERRIDES:
        raise ValueError(f"unsupported NCBI translation table: {translation_table}")
    table = dict(STANDARD_CODON_TABLE)
    table.update(TABLE_OVERRIDES[translation_table])
    return table


def translate_cds(sequence: str, translation_table: int = 1) -> str:
    sequence = sequence.upper()
    table = get_codon_table(translation_table)
    return "".join(
        table.get(sequence[index : index + 3], "X")
        for index in range(0, len(sequence) - 2, 3)
    )


def _strip_terminal_stop(sequence: str) -> str:
    return sequence[:-1] if sequence.endswith("*") else sequence


def validate_cds_pep(
    cds: str,
    pep: str,
    *,
    translation_table: int = 1,
) -> TranslationCheck:
    cds = cds.strip().upper().replace("-", "")
    pep = pep.strip().upper().replace("-", "")
    if not cds:
        return TranslationCheck(False, "empty_cds", False, 0, "")
    if not pep:
        return TranslationCheck(False, "empty_protein", False, 0, "")
    if not set(cds) <= DNA_ALPHABET:
        return TranslationCheck(False, "invalid_cds_alphabet", False, 0, "")
    if not set(pep) <= PROTEIN_ALPHABET:
        return TranslationCheck(False, "invalid_protein_alphabet", False, 0, "")
    if len(cds) % 3:
        return TranslationCheck(False, "cds_length_not_multiple_of_3", False, 0, "")

    try:
        translated = translate_cds(cds, translation_table)
    except ValueError:
        return TranslationCheck(False, "unsupported_translation_table", False, 0, "")
    translated_clean = _strip_terminal_stop(translated)
    protein_clean = _strip_terminal_stop(pep)
    if "*" in translated_clean:
        return TranslationCheck(False, "internal_stop", False, 0, translated)
    if "*" in protein_clean:
        return TranslationCheck(False, "internal_stop", False, 0, translated)
    if len(translated_clean) != len(protein_clean):
        return TranslationCheck(False, "translation_length_mismatch", False, 0, translated)

    mismatches = sum(
        1
        for translated_aa, protein_aa in zip(translated_clean, protein_clean)
        if not (
            translated_aa == protein_aa
            or translated_aa == "X"
            or protein_aa in {"X", "B", "Z", "J"}
        )
    )
    if mismatches:
        return TranslationCheck(False, "translation_mismatch", False, mismatches, translated)
    return TranslationCheck(
        True,
        "",
        translated_clean == protein_clean,
        mismatches,
        translated,
    )

