from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from pgci_preprocess.annotation_convert import (
    ConversionError,
    convert_annotation,
    detect_conversion_format,
)
from pgci_preprocess.fasta import open_text, read_cds_fasta, read_pep_fasta
from pgci_preprocess.gff import parse_gff
from pgci_preprocess.pipeline import PipelineOptions, SpeciesInput, process_species

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if not SRC_DIR.is_dir():
    SRC_DIR = PROJECT_ROOT / "scripts"


def _sequence() -> str:
    return "A" * 100 + "ATGGCTTAA" + "A" * 91 + "TTAAGCCAT" + "A" * 91


def _write_gtf(path: Path) -> None:
    rows = [
        "SEQ1\tEMBL\tgene\t101\t109\t.\t+\t.\tgene_id \"G1\"; gene_biotype \"protein_coding\";",
        "SEQ1\tEMBL\ttranscript\t101\t109\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\"; transcript_biotype \"mRNA\";",
        "SEQ1\tEMBL\texon\t101\t109\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TX1\"; exon_number \"1\";",
        "SEQ1\tEMBL\tCDS\t101\t109\t.\t+\t0\tgene_id \"G1\"; transcript_id \"TX1\"; protein_id \"P1\";",
        "SEQ1\tGnomon\tgene\t201\t209\t.\t-\t.\tgene_id \"G2\"; gene_biotype \"protein_coding\";",
        "SEQ1\tGnomon\texon\t201\t209\t.\t-\t.\tgene_id \"G2\"; transcript_id \"TX2\"; exon_number \"1\";",
        "SEQ1\tGnomon\tCDS\t201\t209\t.\t-\t0\tgene_id \"G2\"; transcript_id \"TX2\"; protein_id \"P2\";",
    ]
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write("#!genome-build fixture\n")
        handle.write("\n".join(rows) + "\n")


def _write_fastas(directory: Path) -> None:
    with gzip.open(directory / "fixture.genome.fa.gz", "wt", encoding="utf-8") as handle:
        handle.write(">SEQ1\n" + _sequence() + "\n")
    (directory / "fixture.cds.fa").write_text(">P1\nATGGCTTAA\n>P2\nATGGCTTAA\n", encoding="utf-8")
    (directory / "fixture.pep.fa").write_text(">P1\nMA\n>P2\nMA\n", encoding="utf-8")


def _gbff_fixture() -> str:
    features = """FEATURES             Location/Qualifiers
     source          1..300
                     /organism="Testus plantus"
                     /mol_type="genomic DNA"
                     /db_xref="taxon:12345"
     gene            101..109
                     /locus_tag="G1"
     mRNA            101..109
                     /locus_tag="G1"
                     /transcript_id="TX1"
     CDS             101..109
                     /locus_tag="G1"
                     /protein_id="P1"
                     /codon_start="1"
                     /translation="MA"
     gene            complement(201..209)
                     /locus_tag="G2"
     mRNA            complement(201..209)
                     /locus_tag="G2"
                     /transcript_id="TX2"
     CDS             complement(201..209)
                     /locus_tag="G2"
                     /protein_id="P2"
                     /codon_start="1"
                     /translation="MA"
"""
    sequence = _sequence()
    origin = "ORIGIN\n"
    for offset in range(0, len(sequence), 60):
        origin += f"{offset + 1:>9} {sequence[offset:offset + 60]}\n"
    return (
        "LOCUS       TESTSEQ1                 300 bp    DNA     linear   PLN 19-NOV-2026\n"
        "DEFINITION  converter fixture.\n"
        "ACCESSION   TESTSEQ1\n"
        "VERSION     TESTSEQ1.1\n"
        "KEYWORDS    .\n"
        + features
        + origin
        + "//\n"
    )


def _run_species(output_dir: Path, prefix: str, project: Path, genome_name: str) -> dict[str, object]:
    (project / "resources").mkdir(parents=True, exist_ok=True)
    species = SpeciesInput(
        show_id="901",
        code="s901",
        species="Testus plantus",
        order="Testales",
        family="Testaceae",
        directory=output_dir,
        gff_name=f"{prefix}.gff3.gz",
        cds_name=f"{prefix}.cds.fa.gz" if (output_dir / f"{prefix}.cds.fa.gz").exists() else "fixture.cds.fa",
        pep_name=f"{prefix}.pep.fa.gz" if (output_dir / f"{prefix}.pep.fa.gz").exists() else "fixture.pep.fa",
        genome_name=genome_name,
        accession="GCA_000000901.1",
    )
    return process_species(
        species,
        PipelineOptions(
            project_root=project,
            min_genes_per_sequence=1,
            accession_policy="legacy",
            overwrite=False,
        ),
    )


def test_gtf_conversion_and_pipeline_roundtrip(tmp_path: Path) -> None:
    source = tmp_path / "fixture.gtf.gz"
    _write_gtf(source)
    _write_fastas(tmp_path)
    assert detect_conversion_format(source) == "gtf"
    result = convert_annotation(source, tmp_path, prefix="fixture")
    assert result.complete is True
    assert result.counts["genes"] == 2
    assert result.counts["transcripts"] == 2
    assert result.counts["transcripts_with_protein_id"] == 2
    assert result.manifest_path.is_file()
    parsed = parse_gff(result.outputs["gff"])
    assert set(parsed.transcripts) == {"TX1", "TX2"}
    assert {next(iter(item.protein_ids)) for item in parsed.transcripts.values()} == {"P1", "P2"}

    run = _run_species(tmp_path, "fixture", tmp_path / "project", "fixture.genome.fa.gz")
    assert run["status"] == "OK", run
    assert run["n_representative_genes"] == 2
    assert run["n_wgdi_genes"] == 2


def test_gtf_output_is_not_overwritten_by_default(tmp_path: Path) -> None:
    source = tmp_path / "fixture.gtf.gz"
    _write_gtf(source)
    output = tmp_path / "converted"
    first = convert_annotation(source, output, prefix="fixture")
    before = first.outputs["gff"].read_bytes()
    with pytest.raises(ConversionError, match="no overwrite"):
        convert_annotation(source, output, prefix="fixture")
    assert first.outputs["gff"].read_bytes() == before


def test_raw_gtf_run_failure_names_the_conversion_route(tmp_path: Path) -> None:
    source = tmp_path / "fixture.gtf.gz"
    _write_gtf(source)
    _write_fastas(tmp_path)
    species = SpeciesInput(
        show_id="902",
        code="s902",
        species="Testus plantus",
        order="Testales",
        family="Testaceae",
        directory=tmp_path,
        gff_name="fixture.gtf.gz",
        cds_name="fixture.cds.fa",
        pep_name="fixture.pep.fa",
        genome_name="fixture.genome.fa.gz",
        accession="GCA_000000902.1",
    )
    result = process_species(
        species,
        PipelineOptions(
            project_root=tmp_path / "project",
            min_genes_per_sequence=1,
            accession_policy="legacy",
        ),
    )
    assert result["status"] == "FAIL"
    assert result["reason"] == "conversion_required"
    assert "convert-annotation" in result["detail"]


def test_transcript_only_gff_with_explicit_protein_id_runs_end_to_end(
    tmp_path: Path,
) -> None:
    with gzip.open(tmp_path / "fixture.gff3.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "SEQ1\tsource\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Name=alpha\n"
            "SEQ1\tsource\tCDS\t1\t9\t.\t+\t0\t"
            "Parent=tx1;protein_id=P1;transl_table=11\n"
        )
    _write_fastas(tmp_path)
    (tmp_path / "fixture.cds.fa").write_text(">P1\nATGGCTTAA\n", encoding="utf-8")
    (tmp_path / "fixture.pep.fa").write_text(">P1\nMA\n", encoding="utf-8")
    result = _run_species(
        tmp_path,
        "fixture",
        tmp_path / "project",
        "fixture.genome.fa.gz",
    )
    assert result["status"] == "OK", result
    assert result["n_wgdi_genes"] == 1


def test_gbff_conversion_and_pipeline_roundtrip(tmp_path: Path) -> None:
    source = tmp_path / "fixture.gbff.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(_gbff_fixture())
    assert detect_conversion_format(source) == "gbff"
    result = convert_annotation(source, tmp_path / "gbff", prefix="fixture")
    assert result.counts["proteins"] == 2
    assert result.counts["transcripts"] == 2
    assert result.counts["genome_records"] == 1
    parsed = parse_gff(result.outputs["gff"])
    assert set(parsed.transcripts) == {"TX1", "TX2"}
    assert parsed.transcripts["TX2"].strand == "-"
    cds = read_cds_fasta(result.outputs["cds"])
    pep = read_pep_fasta(result.outputs["pep"])
    assert cds.sequences == {"P1": "ATGGCTTAA", "P2": "ATGGCTTAA"}
    assert pep.sequences == {"P1": "MA", "P2": "MA"}

    run = _run_species(result.output_dir, "fixture", tmp_path / "project", "fixture.genome.fa.gz")
    assert run["status"] == "OK", run
    assert run["n_representative_genes"] == 2
    assert run["n_wgdi_genes"] == 2


def test_gbff_partial_conversion_is_explicitly_incomplete(tmp_path: Path) -> None:
    source = tmp_path / "fixture.gbff.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(_gbff_fixture())
    result = convert_annotation(source, tmp_path / "partial", prefix="fixture", max_records=1)
    assert result.complete is False
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert manifest["complete"] is False
    assert manifest["counts"]["partial_input"] == 1


def test_gbff_partial_records_without_cds_are_not_full_file_failure(tmp_path: Path) -> None:
    source = tmp_path / "empty.gbff"
    source.write_text(
        "LOCUS       EMPTY1                 9 bp    DNA     linear   PLN 19-NOV-2026\n"
        "DEFINITION  empty fixture.\n"
        "ACCESSION   EMPTY1\n"
        "VERSION     EMPTY1.1\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..9\n"
        "                     /organism=\"Testus plantus\"\n"
        "                     /mol_type=\"genomic DNA\"\n"
        "ORIGIN\n"
        "        1 atggcttaa\n"
        "//\n",
        encoding="utf-8",
    )
    partial = convert_annotation(source, tmp_path / "partial", prefix="empty", max_records=1)
    assert partial.complete is False
    assert partial.counts["proteins"] == 0
    assert partial.counts["partial_without_coding_sequences"] == 1
    with pytest.raises(ConversionError, match="no usable coding"):
        convert_annotation(source, tmp_path / "full", prefix="empty")


def test_gbff_product_variant_disambiguates_identical_cds_parts() -> None:
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from pgci_preprocess.gbff_convert import _Gene, _match_mrna

    cds = SeqFeature(
        FeatureLocation(10, 20, strand=1),
        type="CDS",
        qualifiers={"product": ["example isoform X2"], "protein_id": ["P2"]},
    )
    gene = _Gene("G1|1", "G1", "SEQ1", 0, 100, "+", "G1")
    first = SeqFeature(
        FeatureLocation(0, 100, strand=1),
        type="mRNA",
        qualifiers={"product": ["example transcript variant X1"], "transcript_id": ["TX1"]},
    )
    second = SeqFeature(
        FeatureLocation(0, 100, strand=1),
        type="mRNA",
        qualifiers={"product": ["example transcript variant X2"], "transcript_id": ["TX2"]},
    )
    assert _match_mrna(cds, gene, [("G1|1", first), ("G1|1", second)]) is second
    assert _match_mrna(cds, gene, [("G1|1", first), ("G1|1", second)], {id(second)}) is None


def test_functional_go_table_is_normalized_not_turned_into_gff(tmp_path: Path) -> None:
    source = tmp_path / "species.annotations.csv.gz"
    rows = [
        "gene1\tGO:0000001\t1e-20\t100.5\tdescription, with comma\tTestus plantus\tGO",
        "gene2\tGO:0000002\t2e-10\t88.0\tsecond description\tTestus plantus\tGO",
    ]
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(rows) + "\n")
    assert detect_conversion_format(source) == "functional-go"
    result = convert_annotation(source, tmp_path / "functional", prefix="species")
    assert set(result.outputs) == {"functional_annotations"}
    with open_text(result.outputs["functional_annotations"]) as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    assert table[0] == ["feature_id", "annotation_id", "evalue", "score", "description", "species", "source"]
    assert len(table) == 3
    assert table[1][4] == "description, with comma"
    assert result.counts["features"] == 2


def test_functional_table_preserves_mixed_sources_and_escaped_newlines(tmp_path: Path) -> None:
    source = tmp_path / "mixed.annotations.csv.gz"
    with gzip.open(source, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write("gene1\tXP_1\t1e-5\t50\tfirst line\n")
        handle.write("continued description\tSpecies one\tNR\n")
        handle.write("gene2\tGO:0000001\t2e-5\t60\tgo description\tSpecies one\tGO\n")
    result = convert_annotation(source, tmp_path / "out", prefix="mixed")
    assert result.counts["rows"] == 2
    assert result.counts["continuation_rows"] == 1
    assert result.counts["source_NR"] == 1
    assert result.counts["source_GO"] == 1
    with open_text(result.outputs["functional_annotations"]) as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert rows[1][4] == "first line\\ncontinued description"
    assert rows[1][-1] == "NR"


def test_unknown_comma_csv_fails_closed(tmp_path: Path) -> None:
    source = tmp_path / "unknown.csv"
    source.write_text("gene,GO:0000001,1e-20\n", encoding="utf-8")
    with pytest.raises(ConversionError, match="unrecognized input"):
        detect_conversion_format(source)


def test_locus_tag_gtf_is_detected_across_data_rows(tmp_path: Path) -> None:
    source = tmp_path / "locus_tag.gtf"
    source.write_text(
        "SEQ1\tsource\tgene\t101\t109\t.\t+\t.\tlocus_tag \"G1\";\n"
        "SEQ1\tsource\ttranscript\t101\t109\t.\t+\t.\t"
        "locus_tag \"G1\"; transcript_id \"TX1\";\n"
        "SEQ1\tsource\tCDS\t101\t109\t.\t+\t0\t"
        "locus_tag \"G1\"; transcript_id \"TX1\"; protein_id \"P1\";\n",
        encoding="utf-8",
    )
    assert detect_conversion_format(source) == "gtf"


def test_gtf_with_all_cds_rows_unparented_is_not_complete(tmp_path: Path) -> None:
    source = tmp_path / "broken.gtf.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "SEQ1\tsource\tgene\t101\t109\t.\t+\t.\tgene_id \"G1\";\n"
            "SEQ1\tsource\ttranscript\t101\t109\t.\t+\t.\t"
            "gene_id \"G1\"; transcript_id \"TX1\";\n"
            "SEQ1\tsource\tCDS\t101\t109\t.\t+\t0\tgene_id \"G1\";\n"
        )
    with pytest.raises(ConversionError, match="refusing to emit a CDS-free GFF3"):
        convert_annotation(source, tmp_path / "out", prefix="broken")
    assert not (tmp_path / "out" / "broken.gff3.gz").exists()


def test_gbff_translation_qualifier_must_match_extracted_cds(tmp_path: Path) -> None:
    source = tmp_path / "bad.gbff.gz"
    fixture = _gbff_fixture().replace('/translation="MA"', '/translation="MW"', 1)
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(fixture)
    with pytest.raises(ConversionError, match="translation qualifier are inconsistent"):
        convert_annotation(source, tmp_path / "out", prefix="bad")
    assert not (tmp_path / "out" / "bad.gff3.gz").exists()


def test_convert_annotation_cli_json(tmp_path: Path) -> None:
    source = tmp_path / "fixture.gtf.gz"
    _write_gtf(source)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC_DIR)
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "pgci_preprocess",
            "convert-annotation",
            "--input",
            str(source),
            "--output-dir",
            str(tmp_path / "cli"),
            "--prefix",
            "fixture",
            "--json",
        ],
        cwd=PROJECT_ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=120,
    )
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["status"] == "OK"
    assert payload["format"] == "gtf"
    assert Path(payload["manifest"]).is_file()
