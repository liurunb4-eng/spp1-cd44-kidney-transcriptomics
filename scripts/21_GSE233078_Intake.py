#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

from pathlib import Path
from html.parser import HTMLParser
import csv
import gzip
import shutil
import urllib.request


GSE233078_SUPPL_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233078/suppl/"


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
    with urllib.request.urlopen(GSE233078_SUPPL_URL, timeout=60) as response:
        html_text = response.read().decode("utf-8", errors="replace")
    parser = HrefParser()
    parser.feed(html_text)
    return [href for href in parser.hrefs if href.endswith(".gz")]


def build_download_plan() -> list[dict[str, str]]:
    return [
        {
            "accession": "GSE233078",
            "file_name": "GSE233078_EXPORT_GEO_meta.data.csv.gz",
            "file_type": "scRNA_metadata",
            "url": GSE233078_SUPPL_URL + "GSE233078_EXPORT_GEO_meta.data.csv.gz",
        },
        {
            "accession": "GSE233078",
            "file_name": "GSE233078_EXPORT_GEO_counts_postfilter.rds.gz",
            "file_type": "scRNA_counts_rds",
            "url": GSE233078_SUPPL_URL + "GSE233078_EXPORT_GEO_counts_postfilter.rds.gz",
        },
    ]


def normalize_condition(value: str) -> str:
    text = value.strip().lower()
    if text == "lean":
        return "lean"
    if text == "obese":
        return "obese"
    if text == "obese+enalapril":
        return "obese_enalapril"
    return text.replace("+", "_plus_").replace(" ", "_")


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
        reader = csv.DictReader(handle)
        for row in reader:
            condition = row.get("exp.cond", "").strip()
            rows.append(
                {
                    "cell_id": row.get("", "").strip(),
                    "sample_id": row.get("orig.ident", "").strip(),
                    "condition": condition,
                    "normalized_condition": normalize_condition(condition),
                    "cluster": row.get("clusters3", "").strip(),
                }
            )
    return rows


def inspect_metadata(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["normalized_condition"]] = counts.get(row["normalized_condition"], 0) + 1
    return [
        {"normalized_condition": condition, "cell_count": str(count)}
        for condition, count in sorted(counts.items())
    ]


def build_summary_text(manifest_rows: list[dict[str, str]], metadata_rows: list[dict[str, str]]) -> str:
    return 'Summary written by the analysis script. See generated CSV files and figures for details.\n'

def copy_if_exists(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.exists():
        shutil.copy2(src, dst)


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    data_dir = project_root / "data" / "public_datasets" / "GSE233078"
    raw_dir = data_dir / "processed"
    result_dir = project_root / "results" / "public_data_gse233078"
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

    metadata_rows = parse_metadata(raw_dir / "GSE233078_EXPORT_GEO_meta.data.csv.gz")
    condition_rows = inspect_metadata(metadata_rows)

    write_csv(result_dir / "download_manifest.csv", manifest_rows, list(manifest_rows[0].keys()))
    write_csv(result_dir / "cell_metadata.csv", metadata_rows, list(metadata_rows[0].keys()))
    write_csv(result_dir / "condition_cell_counts.csv", condition_rows, list(condition_rows[0].keys()))
    write_text(result_dir / "gse233078_intake_summary.txt", build_summary_text(manifest_rows, metadata_rows))
    write_text(result_dir / "status.txt", "GSE233078 intake completed.\n")

    print(f"GSE233078 data written to: {data_dir}")
    print(f"GSE233078 results written to: {result_dir}")


if __name__ == "__main__":
    main()


