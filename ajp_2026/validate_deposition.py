from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parent

REQUIRED_SCRIPTS = [
    "scripts/analysis/06_frozen_signature_selection_audit.py",
    "scripts/analysis/07_candidate_sensitivity_random_control.py",
    "scripts/analysis/12_clinical_egfr_stratification.py",
    "scripts/analysis/13_public_outcome_screen.py",
    "scripts/analysis/16_pan_ckd_specificity_audit.py",
    "scripts/analysis/35_analyze_gse137570_ckd_progression.py",
    "scripts/figures/20_main_figure2_clinical_pan_ckd_audit.py",
    "scripts/figures/31_main_figure1_signature_audit_reframed.py",
    "scripts/figures/37_main_figure3_gse137570_patient_level_support.py",
    "scripts/figures/71_rebuild_figure5_marker_excluded.py",
    "scripts/figures/45_main_figure4_acoba_progression_context_bhq.py",
]

REQUIRED_DOCS = [
    "README.md",
    "docs/run_order.md",
    "docs/input_data_manifest.md",
    "requirements_python.txt",
    "requirements_r.txt",
    "source_data/SHA256SUMS.txt",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    errors: list[str] = []

    for relative in REQUIRED_SCRIPTS + REQUIRED_DOCS:
        if not (ROOT / relative).is_file():
            errors.append(f"missing required file: {relative}")

    for folder in [f"source_data/figure{i}" for i in range(1, 6)]:
        path = ROOT / folder
        if not path.is_dir() or not any(item.is_file() for item in path.iterdir()):
            errors.append(f"missing or empty source-data folder: {folder}")

    figure5 = ROOT / "source_data" / "figure5"
    stale_figure5 = sorted(path.name for path in figure5.glob("43_*"))
    if stale_figure5:
        errors.append(f"stale Figure 5 source data remain: {', '.join(stale_figure5)}")

    panel_b = figure5 / "71_figure5_panelB_marker_excluded_correlations.csv"
    if panel_b.is_file():
        with panel_b.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != 8:
            errors.append(f"Figure 5 panel B should contain 8 dataset-marker rows, found {len(rows)}")
        coordinate_types = {row.get("coordinate_type", "") for row in rows}
        if coordinate_types != {"fully_marker_excluded"}:
            errors.append(f"unexpected Figure 5 coordinate types: {sorted(coordinate_types)}")
        if any("marker_removed_from_injury_score" not in row for row in rows):
            errors.append("Figure 5 panel B is missing marker-removal audit fields")

    for relative in REQUIRED_SCRIPTS:
        path = ROOT / relative
        if path.is_file():
            try:
                compile(path.read_text(encoding="utf-8-sig"), str(path), "exec")
            except SyntaxError as exc:
                errors.append(f"syntax error in {relative}: {exc}")

    manifest = ROOT / "source_data" / "SHA256SUMS.txt"
    if manifest.is_file():
        for number, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), start=1):
            if not line.strip():
                continue
            try:
                expected, relative = line.split("  ", 1)
            except ValueError:
                errors.append(f"invalid checksum line {number}")
                continue
            target = ROOT / "source_data" / relative
            if not target.is_file():
                errors.append(f"checksum target missing: {relative}")
            elif sha256(target) != expected:
                errors.append(f"checksum mismatch: {relative}")

    if errors:
        print("AJP 2026 deposition validation: FAIL")
        for error in errors:
            print(f"- {error}")
        raise SystemExit(1)

    print(f"Required scripts: {len(REQUIRED_SCRIPTS)}")
    print("Figure source-data folders: 5")
    print("AJP 2026 deposition validation: PASS")


if __name__ == "__main__":
    main()
