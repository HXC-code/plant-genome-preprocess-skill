"""Derive Golden-109 rebuild inputs from the frozen source manifest.

Outputs (all new files, never overwriting existing deliverables):

- ``download_inventory_golden109.tsv``: canonical inventory for the 109-species
  scope.  File names are the basenames of the verified manifest sources, so an
  isolated project can mirror exactly those files flat into each species
  directory and run the standard ``run`` flow without suffix guessing.
- ``golden109_v3_assembly_map.tsv``: evidence-backed assembly labels for the
  strict accession policy (every species mapped to its frozen index accession,
  evidenced by the verified GFF source hash).
- ``golden109_v3_format_exceptions.tsv``: copy of the evidence-backed format
  switches required by the pipeline (merge_split_genes / gene_span_only).

Deterministic: identical manifest input yields byte-identical outputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import re
import sys
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_ROOT / "src"))

from pgci_preprocess.species_names import load_species_names  # noqa: E402

INVENTORY_FIELDS = ["showid", "species_code", "Order", "Family", "Species", "gff", "genome", "cds", "pep", "rna", "gbff", "gtf"]
ASSEMBLY_FIELDS = ["show_id", "assembly_id", "source", "evidence_path", "evidence_sha256"]
SPECIES_FOLDER_RE = re.compile(r"^\d+_")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    data = buffer.getvalue().encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return _sha256_bytes(data)


def species_folder_name(source_path: str) -> str:
    """The ``<showid>_<Order>_<Family>(<Species>)`` ancestor of a source file."""
    path = Path(source_path)
    for ancestor in path.parents:
        if SPECIES_FOLDER_RE.match(ancestor.name) and ancestor.parent.name in {"downloads", "data"}:
            return ancestor.name
    raise ValueError(f"no species folder ancestor for {source_path}")


def assembly_source_label(accession: str, origin: str) -> str:
    if accession.startswith(("GCF_", "GCA_")):
        return "NCBI/public"
    if accession.startswith("GWH"):
        return "CNCB GWH"
    if "Figshare" in accession:
        return "Figshare/public"
    if "GigaDB" in accession:
        return "GigaDB/public"
    if accession.startswith("CNP"):
        return "CNP/public"
    if "Phytozome" in accession:
        return "Phytozome/public"
    if accession == "JGI" or accession.startswith("JGI_"):
        return "JGI/public"
    return f"curated_label ({origin})"


def derive(manifest_path: Path, output_dir: Path, format_exceptions: Path | None) -> dict[str, object]:
    with manifest_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_species: dict[str, tuple[str, dict[str, dict[str, str]]]] = {}
    for row in rows:
        show_id = row["show_id"].strip()
        if show_id not in by_species:
            by_species[show_id] = (row["code"], {})
        by_species[show_id][1][row["role"]] = row

    inventory: list[dict[str, str]] = []
    assembly: list[dict[str, str]] = []
    problems: list[str] = []
    species_names = load_species_names(PACKAGE_ROOT)
    for show_id in sorted(by_species, key=int):
        manifest_code, roles = by_species[show_id]
        code = (
            species_names.code_for(show_id)
            if species_names is not None
            else manifest_code
        )
        gff = roles["gff"]
        if not gff.get("source_path"):
            problems.append(f"{code}: gff source unresolved")
            continue
        selected_names: dict[str, str] = {}
        for role in ("gff", "cds", "pep", "genome"):
            entry = roles.get(role, {})
            if entry.get("status") == "OMITTED_OPTIONAL" or not entry.get("source_path"):
                selected_names[role] = ""
                continue
            selected_names[role] = Path(entry["source_path"]).name
        if gff.get("content_format") == "GBFF":
            # The gbff adapter extracts all roles from one original file; the
            # isolated rebuild mirrors the canonical derived names instead.
            selected_names = {
                "gff": f"{code}.gff3.gz",
                "cds": f"{code}.cds.fa.gz",
                "pep": f"{code}.pep.fa.gz",
                "genome": f"{code}.genome.fa.gz" if selected_names.get("genome") else "",
            }
        used = [name for name in selected_names.values() if name]
        if len(set(used)) != len(used):
            problems.append(f"{code}: basename collision among selected sources: {used}")
            continue
        inventory.append({
            "showid": show_id,
            "species_code": code,
            "Order": gff["order"],
            "Family": gff["family"],
            "Species": gff["species"],
            "gff": selected_names["gff"],
            "genome": selected_names["genome"],
            "cds": selected_names["cds"],
            "pep": selected_names["pep"],
            "rna": "",
            "gbff": "",
            "gtf": "",
        })
        assembly.append({
            "show_id": show_id,
            "assembly_id": gff["assembly_accession"],
            "source": assembly_source_label(gff["assembly_accession"], gff["source_origin"]),
            "evidence_path": gff["source_path"],
            "evidence_sha256": gff["sha256"],
        })

    hashes = {
        "download_inventory_golden109.tsv": _write_tsv(
            output_dir / "download_inventory_golden109.tsv", INVENTORY_FIELDS, inventory
        ),
        "golden109_v3_assembly_map.tsv": _write_tsv(
            output_dir / "golden109_v3_assembly_map.tsv", ASSEMBLY_FIELDS, assembly
        ),
    }
    if format_exceptions is not None and format_exceptions.is_file():
        data = format_exceptions.read_bytes()
        (output_dir / "golden109_v3_format_exceptions.tsv").write_bytes(data)
        hashes["golden109_v3_format_exceptions.tsv"] = _sha256_bytes(data)
    return {
        "species": len(inventory),
        "assembly_entries": len(assembly),
        "problems": problems,
        "sha256": hashes,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--format-exceptions", type=Path)
    args = parser.parse_args(argv)
    result = derive(args.manifest, args.output_dir, args.format_exceptions)
    for key, value in result.items():
        print(f"{key}={value}")
    return 2 if result["problems"] else 0


if __name__ == "__main__":
    sys.exit(main())
