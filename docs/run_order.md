# Suggested Run Order

The scripts retain the original project numbering. Run scripts from the repository root after downloading the required public datasets and checking path variables near the top of each script.

## 1. Rat Bulk Disease-Context Analyses

- `14_GSE216376_Intake.py`
- `15_GSE216376_ModuleDirection.py`
- `17_GSE183841_Intake.py`
- `18_GSE183841_ModuleDirection.py`

Purpose: define broad fibro-inflammatory and injury-associated transcriptional programs across independent rat kidney injury cohorts.

## 2. GSE233078 Single-Cell Mapping

- `21_GSE233078_Intake.py`
- `22_GSE233078_MetadataAtlas.py`
- `24_GSE233078_ModuleAtlas.py`
- `26_GSE233078_PseudobulkContrast.py`
- `69_GSE233078_DTLAnnotationAudit.py`

Purpose: map disease-associated modules to cell populations and audit the DTL-associated injured tubular source population.

## 3. Focused Ligand-Receptor Analyses

- `32_GSE233078_CellChatMVPPrep.py`
- `34_GSE233078_CellChatFormalMVP.py`
- `36_GSE233078_CellChatStrictRelaxed.py`
- `58_TargetedLRRobustness.py`
- `61_TargetedLRPermutationRobustness.py`
- `62_LIANAFormalFocusedValidation.R`

Purpose: evaluate the Spp1-Cd44 DTL-associated tubular-to-Mono/Macro candidate axis using a focused resource, multi-method support and permutation-based robustness testing.

## 4. State-Context Analyses

- `37_GSE233078_DTL_MacroStateShift.py`
- `38_GSE233078_GraphPseudotime.py`
- `40_GSE233078_Spp1Cd44ModuleSupport.py`
- `42_GSE233078_Spp1Cd44DegEnrichment.py`
- `60_GSE233078_SampleAwareTFPathway.py`

Purpose: place the candidate axis within broader condition-associated tubular and myeloid transcriptional remodeling.

## 5. Human Cohort Analyses

- `57_BulkSignatureProjection.py`
- `63_GSE175759_ClinicalBulkCorrelation.R`
- `67_GSE175759_eGFRStratificationPatch.R`
- `68_GSE30122_HumanDKDBulkValidation.R`

Purpose: evaluate the SPP1/CD44-related program against human eGFR metadata and DKD tubule endpoints.

## 6. Supplementary NicheNet Analyses

- `68_GSE233078_NicheNetLigandTarget.R`
- `70_Human_NicheNet_TargetSetValidation.R`
- `71_NicheNet_SupplementarySpecificity.R`
- `72_NicheNet_TargetSetRandomControl.R`

Purpose: generate the supplementary NicheNet receiver-target layer, human target-set projection, GSE30122 tubular-versus-glomerular specificity check and expression-matched random gene-set control.

NicheNet outputs should be interpreted as receiver-context support and target-set projection, not as functional validation of SPP1/CD44 signaling.

## 7. Human sc/snRNA Compartment-Level Corroboration

- `73_HumanScSnRNACompartmentCorroboration.R`

Purpose: summarize bounded human sc/snRNA compartment-level support from processed GSE131882, GSE195460 and GSE211785 results. This layer is intended as corroboration and boundary assessment, not as a second full single-cell discovery analysis or a ligand-receptor validation.


