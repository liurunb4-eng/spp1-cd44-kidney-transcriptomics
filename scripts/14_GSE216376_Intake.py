#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
from html.parser import HTMLParser
import csv
import gzip
import shutil
import urllib.request


GSE216376_SUPPL_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216376/suppl/"


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
    with urllib.request.urlopen(GSE216376_SUPPL_URL, timeout=60) as response:
        html_text = response.read().decode("utf-8", errors="replace")
    parser = HrefParser()
    parser.feed(html_text)
    return [href for href in parser.hrefs if href.endswith(".gz")]


def build_download_plan() -> list[dict[str, str]]:
    # GEO currently exposes two processed FPKM matrices for this series.
    return [
        {
            "accession": "GSE216376",
            "file_name": "GSE216376_genes_fpkm_expression_Control_and_adenine_samples_.txt.gz",
            "file_type": "processed_fpkm",
            "model_arm": "adenine",
            "url": GSE216376_SUPPL_URL + "GSE216376_genes_fpkm_expression_Control_and_adenine_samples_.txt.gz",
        },
        {
            "accession": "GSE216376",
            "file_name": "GSE216376_genes_fpkm_expression_Sham_and_UUO_samples_.txt.gz",
            "file_type": "processed_fpkm",
            "model_arm": "uuo",
            "url": GSE216376_SUPPL_URL + "GSE216376_genes_fpkm_expression_Sham_and_UUO_samples_.txt.gz",
        },
    ]


def build_sample_metadata() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for i in range(1, 4):
        rows.append(
            {
                "sample_id": f"adenine_control_{i}",
                "model_arm": "adenine",
                "group": "control",
                "expected_matrix_file": "GSE216376_genes_fpkm_expression_Control_and_adenine_samples_.txt.gz",
            }
        )
        rows.append(
            {
                "sample_id": f"adenine_model_{i}",
                "model_arm": "adenine",
                "group": "adenine",
                "expected_matrix_file": "GSE216376_genes_fpkm_expression_Control_and_adenine_samples_.txt.gz",
            }
        )
    for i in range(1, 4):
        rows.append(
            {
                "sample_id": f"uuo_sham_{i}",
                "model_arm": "uuo",
                "group": "sham",
                "expected_matrix_file": "GSE216376_genes_fpkm_expression_Sham_and_UUO_samples_.txt.gz",
            }
        )
        rows.append(
            {
                "sample_id": f"uuo_model_{i}",
                "model_arm": "uuo",
                "group": "uuo",
                "expected_matrix_file": "GSE216376_genes_fpkm_expression_Sham_and_UUO_samples_.txt.gz",
            }
        )
    return rows


def download_file(url: str, out_path: Path) -> str:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists() and out_path.stat().st_size > 0:
        return "already_exists"
    with urllib.request.urlopen(url, timeout=180) as response, out_path.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    return "downloaded"


def inspect_fpkm_file(path: Path) -> dict[str, str]:
    gene_count = 0
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
            gene_count += 1
    return {
        "file_name": path.name,
        "gene_rows": str(gene_count),
        "column_count": str(len(header)),
        "first_column": header[0] if header else "",
        "sample_columns": ";".join(header[1:]),
        "first_gene": first_gene,
    }


def load_gene_symbols_from_fpkm(path: Path) -> set[str]:
    genes: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            if not line.strip():
                continue
            gene = line.split("\t", 1)[0].strip()
            if gene:
                genes.add(gene.upper())
    return genes


def target_genes() -> list[dict[str, str]]:
    return [
        {"module": "vascular_microcirculation", "gene": "NOS3", "rat_alias": "NOS3"},
        {"module": "vascular_microcirculation", "gene": "VEGFA", "rat_alias": "VEGFA"},
        {"module": "vascular_microcirculation", "gene": "HIF1A", "rat_alias": "HIF1A"},
        {"module": "vascular_microcirculation", "gene": "AKT1", "rat_alias": "AKT1"},
        {"module": "fibrosis_ecm", "gene": "TGFB1", "rat_alias": "TGFB1"},
        {"module": "fibrosis_ecm", "gene": "COL1A1", "rat_alias": "COL1A1"},
        {"module": "fibrosis_ecm", "gene": "FN1", "rat_alias": "FN1"},
        {"module": "fibrosis_ecm", "gene": "ACTA2", "rat_alias": "ACTA2"},
        {"module": "inflammation", "gene": "TLR4", "rat_alias": "TLR4"},
        {"module": "inflammation", "gene": "NLRP3", "rat_alias": "NLRP3"},
        {"module": "metabolic_stress", "gene": "HAO2", "rat_alias": "HAO2"},
    ]


def build_gene_presence(download_dir: Path, plan_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    file_gene_sets: dict[str, set[str]] = {}
    for row in plan_rows:
        file_path = download_dir / row["file_name"]
        if file_path.exists():
            file_gene_sets[row["model_arm"]] = load_gene_symbols_from_fpkm(file_path)

    rows: list[dict[str, str]] = []
    for target in target_genes():
        presence = {}
        for model_arm, genes in file_gene_sets.items():
            presence[model_arm] = "yes" if target["rat_alias"].upper() in genes else "no"
        rows.append(
            {
                "module": target["module"],
                "gene": target["gene"],
                "rat_alias": target["rat_alias"],
                "present_in_adenine": presence.get("adenine", "not_checked"),
                "present_in_uuo": presence.get("uuo", "not_checked"),
            }
        )
    return rows


def build_summary_text(manifest_rows: list[dict[str, str]], inventory_rows: list[dict[str, str]], presence_rows: list[dict[str, str]]) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "data" / "public_datasets" / "GSE216376"
    raw_dir = data_dir / "processed"
    result_dir = project_root / "results" / "public_data_gse216376"
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

    metadata_rows = build_sample_metadata()
    inventory_rows = [inspect_fpkm_file(raw_dir / row["file_name"]) for row in plan_rows]
    presence_rows = build_gene_presence(raw_dir, plan_rows)

    write_csv(result_dir / "download_manifest.csv", manifest_rows, list(manifest_rows[0].keys()))
    write_csv(result_dir / "sample_metadata.csv", metadata_rows, list(metadata_rows[0].keys()))
    write_csv(result_dir / "expression_file_inventory.csv", inventory_rows, list(inventory_rows[0].keys()))
    write_csv(result_dir / "module_gene_presence.csv", presence_rows, list(presence_rows[0].keys()))
    write_text(result_dir / "gse216376_intake_summary.txt", build_summary_text(manifest_rows, inventory_rows, presence_rows))
    write_text(result_dir / "status.txt", "GSE216376 intake completed.\n")

    print(f"GSE216376 data written to: {data_dir}")
    print(f"GSE216376 results written to: {result_dir}")


if __name__ == "__main__":
    main()


