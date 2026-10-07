# AJP 2026 Code Deposition Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a self-contained, versioned AJP 2026 supplement to the existing public kidney-transcriptomics repository without changing the historical v1.0.1 release.

**Architecture:** Keep the legacy SPP1/CD44 deposition intact and add an `ajp_2026/` subtree containing the AJP-specific analysis scripts, figure source data, environment records, and run documentation. Work on a dedicated `ajp-2026-update` branch, validate the deposited files locally, then push only that branch for user review and release creation.

**Tech Stack:** Git, GitHub, Python 3, pandas, NumPy, SciPy, scikit-learn, matplotlib, R where retained by the original analysis scripts.

---

### Task 1: Create the isolated update branch

**Files:**
- Modify: Git branch metadata only

**Step 1: Verify the checkout is clean and synchronized**

Run: `git status --short --branch && git fetch origin && git rev-parse HEAD && git rev-parse origin/main`

Expected: clean `main` checkout and matching commit IDs.

**Step 2: Create the branch**

Run: `git switch -c ajp-2026-update`

Expected: branch `ajp-2026-update` checked out from `origin/main`.

### Task 2: Add the AJP-specific analysis and figure scripts

**Files:**
- Create: `ajp_2026/scripts/analysis/06_frozen_signature_selection_audit.py`
- Create: `ajp_2026/scripts/analysis/07_candidate_sensitivity_random_control.py`
- Create: `ajp_2026/scripts/analysis/12_clinical_egfr_stratification.py`
- Create: `ajp_2026/scripts/analysis/13_public_outcome_screen.py`
- Create: `ajp_2026/scripts/analysis/16_pan_ckd_specificity_audit.py`
- Create: `ajp_2026/scripts/analysis/35_analyze_gse137570_ckd_progression.py`
- Create: `ajp_2026/scripts/figures/31_main_figure1_signature_audit_reframed.py`
- Create: `ajp_2026/scripts/figures/20_main_figure2_clinical_pan_ckd_audit.py`
- Create: `ajp_2026/scripts/figures/37_main_figure3_gse137570_patient_level_support.py`
- Create: `ajp_2026/scripts/figures/45_main_figure4_acoba_progression_context_bhq.py`
- Create: `ajp_2026/scripts/figures/43_rebuild_figure5_non_circular.py`

**Step 1: Copy the exact scripts used for the AJP analysis**

Use the current AJP manuscript working directory as the source and preserve the original filenames.

**Step 2: Add a deposition header**

Document that data-root and output-root paths may need to be configured locally and point readers to `ajp_2026/docs/run_order.md`.

**Step 3: Compile-check all Python scripts**

Run: `python -m compileall -q ajp_2026/scripts`

Expected: exit code 0.

### Task 3: Add figure source data and an auditable manifest

**Files:**
- Create: `ajp_2026/source_data/figure1/*`
- Create: `ajp_2026/source_data/figure2/*`
- Create: `ajp_2026/source_data/figure3/*`
- Create: `ajp_2026/source_data/figure4/*`
- Create: `ajp_2026/source_data/figure5/*`
- Create: `ajp_2026/source_data/SHA256SUMS.txt`

**Step 1: Copy the final source-data tables used by Figures 1-5**

Include only derived tabular outputs required to audit the displayed values; do not add raw third-party sequencing matrices.

**Step 2: Generate SHA-256 hashes**

Generate `SHA256SUMS.txt` with paths relative to `ajp_2026/source_data/`.

**Step 3: Verify the manifest**

Recompute every hash and require zero mismatches.

### Task 4: Add AJP documentation and environment records

**Files:**
- Create: `ajp_2026/README.md`
- Create: `ajp_2026/docs/run_order.md`
- Create: `ajp_2026/docs/input_data_manifest.md`
- Create: `ajp_2026/requirements_python.txt`
- Create: `ajp_2026/requirements_r.txt`
- Modify: `README.md`

**Step 1: Document scope and boundaries**

State the frozen 10-gene signature, public datasets, source-data layout, and the distinction between the AJP 2026 supplement and legacy v1.0.1.

**Step 2: Document inputs and run order**

List GEO accessions, external supplementary tables, expected local inputs, and the analysis-to-figure dependency chain.

**Step 3: Update the repository root README**

Add an AJP 2026 section linking to `ajp_2026/README.md` while preserving the legacy manuscript description.

### Task 5: Validate the deposition package

**Files:**
- Create: `ajp_2026/validate_deposition.py`

**Step 1: Implement structural checks**

Check that all 11 required scripts, all five figure source-data folders, documentation files, and the hash manifest are present.

**Step 2: Run validation**

Run: `python ajp_2026/validate_deposition.py`

Expected: `AJP 2026 deposition validation: PASS`.

**Step 3: Inspect Git changes**

Run: `git status --short && git diff --check`

Expected: only intended AJP additions and root README modification; no whitespace errors.

### Task 6: Commit and push the review branch

**Files:**
- Commit all validated deposition files

**Step 1: Commit**

Run: `git add README.md ajp_2026 docs/plans/2026-10-07-ajp-2026-code-deposition.md && git commit -m "Add AJP 2026 analysis deposition"`

Expected: one commit on `ajp-2026-update`.

**Step 2: Push the branch**

Run: `git push -u origin ajp-2026-update`

Expected: remote branch created successfully without modifying `main` or historical tags.

**Step 3: Hand off release actions**

Ask the user to review and merge the branch, create a new GitHub release such as `v2.0.0-ajp`, enable or confirm Zenodo archival, and return the final release URL/DOI for insertion into the AJP Code availability statement.
