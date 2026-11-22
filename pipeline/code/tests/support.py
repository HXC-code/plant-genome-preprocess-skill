from __future__ import annotations

import csv
import json
import os
from pathlib import Path


def build_validation_project(
    tmp_path: Path,
    *,
    species_count: int = 1,
    compatible_links: bool = True,
) -> dict[str, Path]:
    workspace = tmp_path / "workspace"
    project = workspace / "plant-genome-resources"
    resources = project / "resources"
    analysis = project / "analysis_v2"
    resources.mkdir(parents=True)
    analysis.mkdir(parents=True)

    index_path = resources / "wgdi_v2_index.tsv"
    fields = [
        "code",
        "species",
        "accession",
        "wgdi_ready",
        "pep",
        "cds",
        "gff",
        "lens",
        "all_pep",
        "all_cds",
        "all_gff",
        "all_lens",
        "gene_map",
        "chr_map",
        "rejected",
        "qc",
    ]
    rows: list[dict[str, str]] = []
    for index in range(1, species_count + 1):
        code = f"s{index}"
        accession = f"GCF_00000000{index}.1"
        species_dir = analysis / code
        species_dir.mkdir()

        gff = f"1\t{code}1g0001\t10\t100\t+\t1\tG{index}\n"
        lens = "1\t1000\t1\n"
        cds = f">{code}1g0001\nATGGCTTAA\n"
        pep = f">{code}1g0001\nMA\n"

        files = {
            "pep": species_dir / f"{code}.pep",
            "cds": species_dir / f"{code}.cds",
            "gff": species_dir / f"{code}.new.gff",
            "lens": species_dir / f"{code}.lens",
            "gene_map": species_dir / f"{code}.gene_map.tsv",
            "chr_map": species_dir / f"{code}.chr_map.tsv",
            "rejected": species_dir / f"{code}.removed_list.tsv",
            "qc": species_dir / f"{code}.qc.json",
        }
        files["pep"].write_text(pep, encoding="utf-8")
        files["cds"].write_text(cds, encoding="utf-8")
        files["gff"].write_text(gff, encoding="utf-8")
        files["lens"].write_text(lens, encoding="utf-8")
        files["gene_map"].write_text("canonical_gene_id\n", encoding="utf-8")
        files["chr_map"].write_text("source_seqid\n", encoding="utf-8")
        files["rejected"].write_text("scope\tid\treason\tdetail\n", encoding="utf-8")
        files["qc"].write_text(
            json.dumps({"code": code, "status": "OK"}, ensure_ascii=False),
            encoding="utf-8",
        )

        if compatible_links:
            compatibility_dir = workspace / code
            compatibility_dir.mkdir()
            for kind in ("pep", "cds", "gff", "lens"):
                source = files[kind]
                destination = compatibility_dir / source.name
                try:
                    os.link(source, destination)
                except OSError:
                    destination.write_bytes(source.read_bytes())

        rows.append(
            {
                "code": code,
                "species": f"Testus plantus {index}",
                "accession": accession,
                "wgdi_ready": "yes",
                "pep": files["pep"].relative_to(project).as_posix(),
                "cds": files["cds"].relative_to(project).as_posix(),
                "gff": files["gff"].relative_to(project).as_posix(),
                "lens": files["lens"].relative_to(project).as_posix(),
                "all_pep": files["pep"].relative_to(project).as_posix(),
                "all_cds": files["cds"].relative_to(project).as_posix(),
                "all_gff": files["gff"].relative_to(project).as_posix(),
                "all_lens": files["lens"].relative_to(project).as_posix(),
                "gene_map": files["gene_map"].relative_to(project).as_posix(),
                "chr_map": files["chr_map"].relative_to(project).as_posix(),
                "rejected": files["rejected"].relative_to(project).as_posix(),
                "qc": files["qc"].relative_to(project).as_posix(),
            }
        )

    with index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    return {
        "workspace": workspace,
        "project": project,
        "analysis": analysis,
        "index": index_path,
    }


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _fasta(records: dict[str, str]) -> str:
    return "".join(f">{key}\n{sequence}\n" for key, sequence in records.items())


def build_multiformat_project(tmp_path: Path) -> dict[str, Path]:
    """Build four species covering every supported GFF structure.

    * s1 standard gene/mRNA/CDS with two isoforms and a mitochondrion;
    * s2 transcript-only EVM GFF with a genome FASTA used for region backfill;
    * s3 five-column simplified GFF without region rows;
    * s4 standard GFF without region rows, resolved from transcript coordinates.
    """

    project = tmp_path / "plant-genome-resources"
    downloads = project / "data" / "downloads"
    resources = project / "resources"
    resources.mkdir(parents=True)

    # s1: standard NCBI-style hierarchy.
    d1 = downloads / "1_Testales_Testaceae(Testus primus)"
    d1.mkdir(parents=True)
    gff1 = "\n".join(
        [
            "##gff-version 3",
            "1\tRefSeq\tregion\t1\t2000\t.\t+\t.\tID=1;genome=genomic;chromosome=1",
            "1\tRefSeq\tgene\t100\t120\t.\t+\t.\tID=gene1;Name=gene1;gene_biotype=protein_coding",
            "1\tRefSeq\tmRNA\t100\t120\t.\t+\t.\tID=rna-tx1;Parent=gene1;transcript_id=tx1",
            "1\tRefSeq\tCDS\t100\t108\t.\t+\t0\tID=cds-tx1;Parent=rna-tx1;protein_id=prot_tx1",
            "1\tRefSeq\tmRNA\t100\t120\t.\t+\t.\tID=rna-tx1long;Parent=gene1;transcript_id=tx1long",
            "1\tRefSeq\tCDS\t100\t117\t.\t+\t0\tID=cds-tx1long;Parent=rna-tx1long;protein_id=prot_tx1long",
            "1\tRefSeq\tgene\t300\t320\t.\t+\t.\tID=gene2;Name=gene2;gene_biotype=protein_coding",
            "1\tRefSeq\tmRNA\t300\t320\t.\t+\t.\tID=rna-tx2;Parent=gene2;transcript_id=tx2",
            "1\tRefSeq\tCDS\t300\t317\t.\t+\t0\tID=cds-tx2;Parent=rna-tx2;protein_id=prot_tx2",
            "1\tRefSeq\tgene\t600\t620\t.\t+\t.\tID=gene3;Name=gene3;gene_biotype=protein_coding",
            "1\tRefSeq\tmRNA\t600\t620\t.\t+\t.\tID=rna-tx3;Parent=gene3;transcript_id=tx3",
            "1\tRefSeq\tCDS\t600\t608\t.\t+\t0\tID=cds-tx3;Parent=rna-tx3;protein_id=prot_tx3",
            "MT\tRefSeq\tregion\t1\t300\t.\t+\t.\tID=MT;genome=mitochondrion",
            "MT\tRefSeq\tgene\t10\t20\t.\t+\t.\tID=gene4;Name=gene4;gene_biotype=protein_coding",
            "MT\tRefSeq\tmRNA\t10\t20\t.\t+\t.\tID=rna-tx4;Parent=gene4;transcript_id=tx4",
            "MT\tRefSeq\tCDS\t10\t18\t.\t+\t0\tID=cds-tx4;Parent=rna-tx4;protein_id=prot_tx4",
            "",
        ]
    )
    _write_text(d1 / "GCF_000000001.1_test1_genomic.gff", gff1)
    _write_text(
        d1 / "GCF_000000001.1_test1_cds_from_genomic.fna",
        _fasta(
            {
                "prot_tx1": "ATGGCTTAA",
                "prot_tx1long": "ATGGCTGCTGCTGCTTAA",
                "prot_tx2": "ATGGCTGCTGCTGCTTAA",
                "prot_tx3": "ATGGCTTAA",
                "prot_tx4": "ATGAAATAA",
            }
        ),
    )
    _write_text(
        d1 / "GCF_000000001.1_test1_protein.faa",
        _fasta(
            {
                "prot_tx1": "MA",
                "prot_tx1long": "MAAAA",
                "prot_tx2": "MAAAA",
                "prot_tx3": "MA",
                "prot_tx4": "MK",
            }
        ),
    )

    # s2: transcript-only GFF, region rows supplied by a genome FASTA.
    d2 = downloads / "2_Testales_Testaceae(Testus secundus)"
    d2.mkdir(parents=True)
    gff2 = "\n".join(
        [
            "##gff-version 3",
            "1\tEVM\tmRNA\t100\t117\t.\t+\t.\tID=txA;Name=txA",
            "1\tEVM\tCDS\t100\t117\t.\t+\t0\tID=cds.txA;Parent=txA",
            "1\tEVM\tmRNA\t300\t308\t.\t+\t.\tID=txB;Name=txB",
            "1\tEVM\tCDS\t300\t308\t.\t+\t0\tID=cds.txB;Parent=txB",
            "",
        ]
    )
    _write_text(d2 / "GCF_000000002.1_test2_genomic.gff", gff2)
    _write_text(d2 / "GCF_000000002.1_test2_genome.fa", ">1\n" + "N" * 2000 + "\n")
    _write_text(
        d2 / "GCF_000000002.1_test2_cds.fa",
        _fasta({"txA": "ATGGCTGCTGCTGCTTAA", "txB": "ATGAAATAA"}),
    )
    _write_text(
        d2 / "GCF_000000002.1_test2_protein.faa",
        _fasta({"txA": "MAAAA", "txB": "MK"}),
    )

    # s3: five-column simplified GFF.
    d3 = downloads / "3_Testales_Testaceae(Testus tertius)"
    d3.mkdir(parents=True)
    gff3 = "\n".join(
        [
            "1\t100\t114\t+\tgA",
            "1\t300\t308\t+\tgB",
            "2\t100\t108\t-\tgC",
            "",
        ]
    )
    _write_text(d3 / "GCF_000000003.1_test3.gff3", gff3)
    _write_text(
        d3 / "GCF_000000003.1_test3.cds.fa",
        _fasta({"gA": "ATGGCTGCTGCTTAA", "gB": "ATGGCTTAA", "gC": "ATGAAATAA"}),
    )
    _write_text(
        d3 / "GCF_000000003.1_test3_protein.faa",
        _fasta({"gA": "MAAA", "gB": "MA", "gC": "MK"}),
    )

    # s4: standard hierarchy without region rows or genome.
    d4 = downloads / "4_Testales_Testaceae(Testus quartus)"
    d4.mkdir(parents=True)
    gff4 = "\n".join(
        [
            "##gff-version 3",
            "1\tsource\tgene\t100\t120\t.\t+\t.\tID=g1",
            "1\tsource\tmRNA\t100\t120\t.\t+\t.\tID=rna-t1;Parent=g1;transcript_id=t1",
            "1\tsource\tCDS\t100\t108\t.\t+\t0\tID=cds-t1;Parent=rna-t1;protein_id=p1",
            "1\tsource\tmRNA\t100\t120\t.\t+\t.\tID=rna-t1long;Parent=g1;transcript_id=t1long",
            "1\tsource\tCDS\t100\t117\t.\t+\t0\tID=cds-t1long;Parent=rna-t1long;protein_id=p1long",
            "2\tsource\tgene\t50\t58\t.\t-\t.\tID=g2",
            "2\tsource\tmRNA\t50\t58\t.\t-\t.\tID=rna-t2;Parent=g2;transcript_id=t2",
            "2\tsource\tCDS\t50\t58\t.\t-\t0\tID=cds-t2;Parent=rna-t2;protein_id=p2",
            "",
        ]
    )
    _write_text(d4 / "GCF_000000004.1_test4_genomic.gff", gff4)
    _write_text(
        d4 / "GCF_000000004.1_test4_cds.fa",
        _fasta({"p1": "ATGGCTTAA", "p1long": "ATGGCTGCTGCTGCTTAA", "p2": "ATGAAATAA"}),
    )
    _write_text(
        d4 / "GCF_000000004.1_test4_protein.faa",
        _fasta({"p1": "MA", "p1long": "MAAAA", "p2": "MK"}),
    )

    inventory_rows = [
        {
            "showid": "1",
            "Order": "Testales",
            "Family": "Testaceae",
            "Species": "Testus primus",
            "gff": "GCF_000000001.1_test1_genomic.gff",
            "genome": "",
            "cds": "GCF_000000001.1_test1_cds_from_genomic.fna",
            "pep": "GCF_000000001.1_test1_protein.faa",
            "rna": "",
            "gbff": "",
            "gtf": "",
        },
        {
            "showid": "2",
            "Order": "Testales",
            "Family": "Testaceae",
            "Species": "Testus secundus",
            "gff": "GCF_000000002.1_test2_genomic.gff",
            "genome": "GCF_000000002.1_test2_genome.fa",
            "cds": "GCF_000000002.1_test2_cds.fa",
            "pep": "GCF_000000002.1_test2_protein.faa",
            "rna": "",
            "gbff": "",
            "gtf": "",
        },
        {
            "showid": "3",
            "Order": "Testales",
            "Family": "Testaceae",
            "Species": "Testus tertius",
            "gff": "GCF_000000003.1_test3.gff3",
            "genome": "",
            "cds": "GCF_000000003.1_test3.cds.fa",
            "pep": "GCF_000000003.1_test3_protein.faa",
            "rna": "",
            "gbff": "",
            "gtf": "",
        },
        {
            "showid": "4",
            "Order": "Testales",
            "Family": "Testaceae",
            "Species": "Testus quartus",
            "gff": "GCF_000000004.1_test4_genomic.gff",
            "genome": "",
            "cds": "GCF_000000004.1_test4_cds.fa",
            "pep": "GCF_000000004.1_test4_protein.faa",
            "rna": "",
            "gbff": "",
            "gtf": "",
        },
    ]
    fields = ["showid", "Order", "Family", "Species", "gff", "genome", "cds", "pep", "rna", "gbff", "gtf"]
    with (resources / "download_inventory.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(inventory_rows)
    return {
        "project": project,
        "downloads": downloads,
        "directories": {"1": d1, "2": d2, "3": d3, "4": d4},
    }


