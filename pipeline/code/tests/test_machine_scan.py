"""Tests for machine/data discovery (``pgci_preprocess.machine_scan``).

The scanner is what makes a copied skill work on a new computer, so these
tests pin the behaviours that previously produced wrong inventories: organelle
false positives, reference databases treated as PEP, container folders treated
as species, and mixed-assembly bundles.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from pgci_preprocess import machine_scan
from pgci_preprocess.project_paths import downloads_root

SPECIES = "123_Testales_Testaceae(Fakea testa)"


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _gff(path: Path, seqid: str = "Chr1", gene: str = "gene1") -> Path:
    rows = [
        "##gff-version 3",
        f"{seqid}\ttest\tgene\t10\t200\t.\t+\t.\tID={gene};Name={gene}",
        f"{seqid}\ttest\tmRNA\t10\t200\t.\t+\t.\tID={gene}.1;Parent={gene}",
        f"{seqid}\ttest\texon\t10\t100\t.\t+\t.\tID={gene}.1.exon1;Parent={gene}.1",
        f"{seqid}\ttest\tCDS\t10\t100\t.\t+\t.\tID=cds-{gene}.1;Parent={gene}.1",
    ]
    return _write(path, "\n".join(rows) + "\n")


def _fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    body = "".join(f">{header}\n{sequence}\n" for header, sequence in records)
    return _write(path, body)


def build_species(root: Path, seqid: str = "Chr1", gene: str = "gene1") -> Path:
    directory = root / SPECIES
    _gff(directory / "annotation" / "asm.genes.gff", seqid=seqid, gene=gene)
    _fasta(
        directory / "assembly" / "asm.fa",
        [(seqid, "ACGT" * 60)],
    )
    _fasta(directory / "annotation" / "asm.cds.fa", [(f"{gene}.1", "ATGAAATAG")])
    _fasta(directory / "annotation" / "asm.pep.fa", [(f"{gene}.1", "MK")])
    return directory


# ---------------------------------------------------------------------------
# naming / ids
# ---------------------------------------------------------------------------


def test_parse_species_dirname() -> None:
    parsed = machine_scan.parse_species_dirname(SPECIES)
    assert parsed == {
        "show_id": "123",
        "order": "Testales",
        "family": "Testaceae",
        "species": "Fakea testa",
    }
    assert machine_scan.parse_species_dirname("not-a-species") is None


def test_assembly_key_keeps_versions_apart() -> None:
    key = machine_scan._assembly_key
    assert key("ppa_v2.asm.genes.gff") == key("ppa_v2.asm.CDS.fasta")
    assert key("ppa_v1.asm.transcriptome.gtf") != key("ppa_v2.asm.genes.gff")
    assert key("GCF_003935025.1_Abrus_2018_cds_from_genomic.fna.gz") == "gcf_003935025.1"
    assert key("B.blakeana_mat.genome.normalized.gff.gz") == key("B.blakeana_mat.genome.cds.fa")


def test_normalise_id_folds_common_spellings() -> None:
    normalise = machine_scan._normalise_id
    assert normalise("cds-XP_054824436.1") == normalise("XP_054824436.1")
    assert normalise("lcl|NW_026555332.1_cds_XP_054824436.1_3") == "xp_054824436"
    assert normalise("NPAL_000001_R1") == normalise("NPAL_000001")


# ---------------------------------------------------------------------------
# file classification
# ---------------------------------------------------------------------------


def test_nuclear_cds_with_organelle_product_is_not_organelle(tmp_path: Path) -> None:
    """Product descriptions mention 'chloroplastic'; that is not a genome."""

    path = _fasta(
        tmp_path / "asm.cds.fa",
        [
            ("g1", "ATGAAATAG"),
            ("g2", "ATGGGGTAG"),
        ],
    )
    _write(
        path,
        path.read_text(encoding="utf-8").replace(
            ">g2", ">g2 biotin carboxylase 1, chloroplastic"
        ),
    )
    record = machine_scan.classify_file(path, tmp_path)
    assert record.role == "cds"


def test_true_organelle_fasta_is_excluded(tmp_path: Path) -> None:
    path = _fasta(
        tmp_path / "cp.genome.fa",
        [("CP1 Chloroplast Complete=T Circular=T Len=126696", "ACGT" * 40)],
    )
    record = machine_scan.classify_file(path, tmp_path)
    assert record.role == "organelle"
    assert record.score < 0


def test_reference_database_is_not_pep(tmp_path: Path) -> None:
    path = _fasta(
        tmp_path / "Annotation_SwissProtPlants_27032024.fasta",
        [("P12345", "MKV")],
    )
    record = machine_scan.classify_file(path, tmp_path)
    assert record.role == "reference"


def test_simplified_five_column_gff_is_a_valid_annotation(tmp_path: Path) -> None:
    path = _write(tmp_path / "manual.gff3", "Chr1\t10\t200\t+\tgene1\nChr1\t300\t400\t-\tgene2\n")
    record = machine_scan.classify_file(path, tmp_path)
    assert record.role == "gff"
    assert record.score >= 0, record.note
    assert "gene_models" in record.note


def test_repeat_track_ranks_below_real_annotation(tmp_path: Path) -> None:
    _gff(tmp_path / "asm.genes.gff")
    _write(tmp_path / "asm.repeats.gff3", "Chr1\ttest\tmatch\t1\t50\t.\t+\t.\tID=r1\n")
    real = machine_scan.classify_file(tmp_path / "asm.genes.gff", tmp_path)
    repeat = machine_scan.classify_file(tmp_path / "asm.repeats.gff3", tmp_path)
    assert real.score > repeat.score


# ---------------------------------------------------------------------------
# bundle selection
# ---------------------------------------------------------------------------


def test_bundle_keeps_the_assembly_that_matches_the_annotation(tmp_path: Path) -> None:
    directory = build_species(tmp_path)
    # An unrelated chromosome assembly with the same filenames pattern.
    _fasta(
        directory / "other" / "asm.fa",
        [("scaffold_1", "ACGT" * 60), ("scaffold_2", "ACGT" * 60)],
    )
    files = [machine_scan.classify_file(p, directory) for p in directory.rglob("*") if p.is_file()]
    chosen, note = machine_scan.select_bundle(files)

    assert Path(chosen["genome"].path).name == "asm.fa"
    assert Path(chosen["gff"].path).name == "asm.genes.gff"
    assert "seqid-overlap=100%" in note
    assert set(chosen) >= {"gff", "genome", "cds", "pep"}


def test_bundle_drops_unrelated_transcriptome(tmp_path: Path) -> None:
    directory = build_species(tmp_path)
    _fasta(directory / "annotation" / "other.rna.fa", [("other1", "ACGU")])
    files = [machine_scan.classify_file(p, directory) for p in directory.rglob("*") if p.is_file()]
    chosen, _note = machine_scan.select_bundle(files)
    assert "rna" not in chosen


def test_bundle_leaves_genome_unpaired_when_seqids_disagree(tmp_path: Path) -> None:
    directory = tmp_path / SPECIES
    _gff(directory / "asm.genes.gff", seqid="Chr1")
    _fasta(directory / "asm.fa", [("scaffold_99", "ACGT" * 60)])
    files = [machine_scan.classify_file(p, directory) for p in directory.rglob("*") if p.is_file()]
    chosen, note = machine_scan.select_bundle(files)
    assert "genome" not in chosen
    assert "UNPAIRED" in note


# ---------------------------------------------------------------------------
# discovery of folders
# ---------------------------------------------------------------------------


def test_container_directory_is_not_a_species(tmp_path: Path) -> None:
    data = tmp_path / "data"
    build_species(data)
    found = machine_scan.find_species_dirs([data])
    assert [item.show_id for item in found] == ["123"]

    # Scanning the project root must not swallow data/ as one fake species.
    from_root = machine_scan.find_species_dirs([tmp_path])
    assert [item.show_id for item in from_root] == ["123"]


def test_unparseable_folder_falls_back_to_content(tmp_path: Path) -> None:
    directory = tmp_path / "999_weird_name"
    _gff(directory / "a.gff")
    _fasta(directory / "a.fa", [("Chr1", "ACGT" * 60)])
    found = machine_scan.find_species_dirs([tmp_path])
    assert len(found) == 1
    assert found[0].strict_name is False
    assert found[0].warnings


# ---------------------------------------------------------------------------
# config / pointer / data root
# ---------------------------------------------------------------------------


def test_build_config_writes_inventory_pointer_and_data_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    data = tmp_path / "project" / "data"
    build_species(data)
    pointer_dir = tmp_path / "pointer"
    monkeypatch.setenv("PGCI_POINTER_DIR", str(pointer_dir))

    config = machine_scan.build_config([data])
    project_root = Path(config["paths"]["project_root"])
    assert project_root == tmp_path / "project"
    assert Path(config["paths"]["data_root"]) == data

    machine_scan.seed_resources(project_root)
    machine_scan.write_discovery_report(config, project_root)
    machine_scan.write_sequence_mode_map(config, project_root)
    inventory = machine_scan.write_inventory(config, project_root, overwrite=True)
    config_path = machine_scan.write_config(
        config, project_root / "resources" / machine_scan.CONFIG_FILENAME
    )
    pointer = machine_scan.write_pointer(project_root, config_path)

    assert pointer.is_file()
    assert machine_scan.read_pointer()["project_root"] == str(project_root)
    rows = inventory.read_text(encoding="utf-8").splitlines()
    assert rows[0].startswith("showid\t")
    assert len(rows) == 2
    assert (project_root / "resources" / "sequence_mode.tsv").is_file()
    assert machine_scan.read_sequence_mode_map(project_root)["123"] == "chromosomes"


def test_downloads_root_honours_recorded_data_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    project = tmp_path / "project"
    data = tmp_path / "elsewhere" / "genomes"
    _gff(data / SPECIES / "a.gff")
    (project / "resources").mkdir(parents=True)
    (project / "resources" / "machine_paths.json").write_text(
        json.dumps(
            {
                "schema_version": machine_scan.SCHEMA_VERSION,
                "paths": {"project_root": str(project), "data_root": str(data)},
            }
        ),
        encoding="utf-8",
    )
    assert downloads_root(project) == data.resolve()


def test_write_inventory_does_not_clobber_reviewed_rows(tmp_path: Path) -> None:
    project = tmp_path / "project"
    target = project / "resources" / "download_inventory.tsv"
    _write(target, "showid\tspecies_code\n1\thand-reviewed\n")
    config = {"suggested_inventory": [{"showid": "2"}]}
    machine_scan.write_inventory(config, project)
    assert "hand-reviewed" in target.read_text(encoding="utf-8")


def test_discover_cli_creates_project_from_bare_data(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from pgci_preprocess.cli import main

    project = tmp_path / "project"
    build_species(project / "data")
    monkeypatch.setenv("PGCI_POINTER_DIR", str(tmp_path / "pointer"))
    monkeypatch.delenv("PGCI_PROJECT_ROOT", raising=False)
    monkeypatch.chdir(project)

    assert main(["discover"]) == 0
    assert (project / "resources" / "download_inventory.tsv").is_file()
    assert (project / "resources" / "machine_paths.json").is_file()
    assert (project / "resources" / "sequence_mode.tsv").is_file()
    rows = (project / "resources" / "download_inventory.tsv").read_text(
        encoding="utf-8"
    ).splitlines()
    assert len(rows) == 2
    assert machine_scan.read_pointer()["project_root"] == str(project.resolve())

    from pgci_preprocess.doctor import run_doctor

    report = run_doctor(project)
    checks = {item["name"]: item for item in report["checks"]}
    assert checks["machine_paths"]["status"] == "PASS", checks["machine_paths"]


def test_discover_cli_reports_when_nothing_matches(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from pgci_preprocess.cli import main

    monkeypatch.setenv("PGCI_POINTER_DIR", str(tmp_path / "pointer"))
    empty = tmp_path / "empty"
    empty.mkdir()
    monkeypatch.chdir(empty)
    assert main(["discover"]) == 1
