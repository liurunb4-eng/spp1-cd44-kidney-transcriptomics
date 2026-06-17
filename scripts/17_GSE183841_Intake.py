#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
from html.parser import HTMLParser
import csv
import gzip
import shutil
import urllib.request


GSE183841_SUPPL_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183841/suppl/"


class HrefParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.hrefs.append(value)


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8-sig")


def fetch_supplementary_file_names() -> list[str]:
    with urllib.request.urlopen(GSE183841_SUPPL_URL, timeout=60) as response:
        html_text = response.read().decode("utf-8", errors="replace")
    parser = HrefParser()
    parser.feed(html_text)
    return [href for href in parser.hrefs if href.endswith(".gz")]


def build_download_plan() -> list[dict[str, str]]:
    return [
        {
            "accession": "GSE183841",
            "file_name": "GSE183841_EXPORT_BulkRNAseq_TPM_Counts.txt.gz",
            "file_type": "bulk_tpm",
            "url": GSE183841_SUPPL_URL + "GSE183841_EXPORT_BulkRNAseq_TPM_Counts.txt.gz",
        },
        {
            "accession": "GSE183841",
            "file_name": "GSE183841_EXPORT_BulkRNAseq_Raw_Counts.txt.gz",
            "file_type": "bulk_raw_counts",
            "url": GSE183841_SUPPL_URL + "GSE183841_EXPORT_BulkRNAseq_Raw_Counts.txt.gz",
        },
        {
            "accession": "GSE183841",
            "file_name": "GSE183841_EXPORT_BulkRNAseq_metadata.txt.gz",
            "file_type": "bulk_metadata",
            "url": GSE183841_SUPPL_URL + "GSE183841_EXPORT_BulkRNAseq_metadata.txt.gz",
        },
    ]


def normalize_treatment(value: str) -> str:
    text = value.strip().lower()
    if text == "control":
        return "control"
    if text == "doca":
        return "doca"
    if "finerenone" in text:
        return "finerenone"
    if "spironolactone" in text:
        return "spironolactone"
    return text.replace(" ", "_")


def build_target_gene_map() -> list[dict[str, str]]:
    return [
        {
            "module": "vascular_microcirculation",
            "gene": "NOS3",
            "rat_ensembl_id": "ENSRNOG00000009348",
            "mapping_source": "Ensembl rat symbol lookup",
        },
        {
            "module": "vascular_microcirculation",
            "gene": "VEGFA",
            "rat_ensembl_id": "ENSRNOG00000019598",
            "mapping_source": "Ensembl rat symbol lookup",
        },
        {
            "module": "vascular_microcirculation",
            "gene": "HIF1A",
            "rat_ensembl_id": "ENSRNOG00000008292",
            "mapping_source": "Ensembl rat symbol lookup",
        },
        {
            "module": "vascular_microcirculation",
            "gene": "AKT1",
            "rat_ensembl_id": "ENSRNOG00000021497",
            "mapping_source": "Ensembl rat symbol lookup (preferred canonical ID)",
        },
        {
            "module": "fibrosis_ecm",
            "gene": "TGFB1",
            "rat_ensembl_id": "ENSRNOG00000020652",
            "mapping_source": "Ensembl rat symbol lookup",
        },
        {
            "module": "fibrosis_ecm",
            "gene": "COL1A1",
            "rat_ensembl_id": "ENSRNOG00000003897",
            "mapping_source": "Ensembl rat symbol lookup",
        },
        {
            "module": "fibrosis_ecm",
            "gene": "FN1",
            "rat_ensembl_id": "ENSRNOG00000014288",
            "mapping_source": "Ensembl rat symbol lookup",
        },
        {
            "module": "fibrosis_ecm",
            "gene": "ACTA2",
            "rat_ensembl_id": "ENSRNOG00000058039",
            "mapping_source": "Ensembl rat symbol lookup",
        },
        {
            "module": "inflammation",
            "gene": "TLR4",
            "rat_ensembl_id": "ENSRNOG00000010522",
            "mapping_source": "Ensembl rat symbol lookup",
        },
        {
            "module": "inflammation",
            "gene": "NLRP3",
            "rat_ensembl_id": "ENSRNOG00000003170",
            "mapping_source": "Ensembl rat symbol lookup",
        },
        {
            "module": "metabolic_stress",
            "gene": "HAO2",
            "rat_ensembl_id": "ENSRNOG00000019470",
            "mapping_source": "Ensembl rat symbol lookup",
        },
    ]


def download_file(url: str, out_path: Path) -> str:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        return "already_exists"
    with urllib.request.urlopen(url, timeout=180) as response, out_path.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    return "downloaded"


def parse_metadata(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            treatment = row.get("treatment", "").strip()
            rows.append(
                {
                    "sample_id": row.get("Samples", "").strip(),
                    "sample_title": row.get("title", "").strip(),
                    "treatment": treatment,
                    "normalized_group": normalize_treatment(treatment),
                }
            )
    return rows


def inspect_expression_file(path: Path) -> dict[str, str]:
    gene_rows = 0
    header: list[str] = []
    first_gene = ""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle):
            parts = line.rstrip("\n").split("\t")
            if line_no == 0:
                header = parts
                continue
            if parts and not first_gene:
                first_gene = parts[0]
            gene_rows += 1
    return {
        "file_name": path.name,
        "gene_rows": str(gene_rows),
        "column_count": str(len(header)),
        "first_column": header[0] if header else "",
        "sample_columns": ";".join(header[1:]),
        "first_gene": first_gene,
    }


def load_gene_ids(path: Path) -> set[str]:
    gene_ids: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            if not line.strip():
                continue
            gene_ids.add(line.split("\t", 1)[0].strip().upper())
    return gene_ids


def build_gene_presence(raw_dir: Path, target_map: list[dict[str, str]]) -> list[dict[str, str]]:
    tpm_ids = load_gene_ids(raw_dir / "GSE183841_EXPORT_BulkRNAseq_TPM_Counts.txt.gz")
    raw_ids = load_gene_ids(raw_dir / "GSE183841_EXPORT_BulkRNAseq_Raw_Counts.txt.gz")
    rows: list[dict[str, str]] = []
    for target in target_map:
        ensembl_id = target["rat_ensembl_id"].upper()
        rows.append(
            {
                "module": target["module"],
                "gene": target["gene"],
                "rat_ensembl_id": target["rat_ensembl_id"],
                "present_in_tpm": "yes" if ensembl_id in tpm_ids else "no",
                "present_in_raw_counts": "yes" if ensembl_id in raw_ids else "no",
            }
        )
    return rows


def build_summary_text(
    manifest_rows: list[dict[str, str]],
    inventory_rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
    presence_rows: list[dict[str, str]],
) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "data" / "public_datasets" / "GSE183841"
    raw_dir = data_dir / "processed"
    result_dir = project_root / "results" / "public_data_gse183841"
    result_dir.mkdir(parents=True, exist_ok=True)

    plan_rows = build_download_plan()
    supplementary_names = fetch_supplementary_file_names()
    manifest_rows: list[dict[str, str]] = []
    for row in plan_rows:
        status = download_file(row["url"], raw_dir / row["file_name"])
        manifest_rows.append(
            {
                **row,
                "download_status": status,
                "found_in_supplementary_listing": "yes" if row["file_name"] in supplementary_names else "no",
                "local_path": str(raw_dir / row["file_name"]),
            }
        )

    metadata_rows = parse_metadata(raw_dir / "GSE183841_EXPORT_BulkRNAseq_metadata.txt.gz")
    inventory_rows = [
        inspect_expression_file(raw_dir / "GSE183841_EXPORT_BulkRNAseq_TPM_Counts.txt.gz"),
        inspect_expression_file(raw_dir / "GSE183841_EXPORT_BulkRNAseq_Raw_Counts.txt.gz"),
    ]
    target_map = build_target_gene_map()
    presence_rows = build_gene_presence(raw_dir, target_map)

    write_csv(result_dir / "download_manifest.csv", manifest_rows, list(manifest_rows[0].keys()))
    write_csv(result_dir / "sample_metadata.csv", metadata_rows, list(metadata_rows[0].keys()))
    write_csv(result_dir / "expression_file_inventory.csv", inventory_rows, list(inventory_rows[0].keys()))
    write_csv(result_dir / "target_gene_ensembl_map.csv", target_map, list(target_map[0].keys()))
    write_csv(result_dir / "module_gene_presence.csv", presence_rows, list(presence_rows[0].keys()))
    write_text(
        result_dir / "gse183841_intake_summary.txt",
        build_summary_text(manifest_rows, inventory_rows, metadata_rows, presence_rows),
    )
    write_text(result_dir / "status.txt", "GSE183841 intake completed.\n")

    print(f"GSE183841 data written to: {data_dir}")
    print(f"GSE183841 results written to: {result_dir}")


if __name__ == "__main__":
    main()


