# AJP 2026 input-data manifest

All analyses use public, de-identified human kidney transcriptomic resources or published supplementary differential-expression tables.

| Resource | Role in the AJP analysis | Redistribution |
| --- | --- | --- |
| GSE175759 | Candidate selection, matched random control, eGFR association | Download from GEO; raw/large matrices not included |
| GSE30122 | DKD primary-tubule audit and matched random control | Download series matrix and GPL571 annotation from GEO |
| GSE180394 | C-PROBE tubulointerstitial pan-CKD context | Download from GEO |
| GSE180393 | C-PROBE glomerular context | Download from GEO |
| GSE137570 | Patient-level progressive CKD grouping, GFR, and TIF | Download normalized workbook from GEO |
| GSE131882 | PT and myeloid cellular context; marker-excluded coordinates use an injury score originally defined from SPP1, HAVCR1, LCN2, KRT8, and KRT18 | Public source; large processed cell table not redistributed |
| GSE195460 | PT-like cellular context; marker-excluded coordinates use an injury score originally defined from SPP1, HAVCR1, LCN2, VCAM1, VIM, KRT8, and KRT18 | Public source; large processed cell table not redistributed |
| GSE211785 | Immune/myeloid CD44 context | Public source; large processed cell table not redistributed |
| Acoba/KaroKidney supplementary tables | DKD progression-associated gene-level alignment | Obtain from the article supplementary files |

The source-data tables under `../source_data/` are derived analytical outputs supporting the manuscript figures. Their SHA-256 hashes are recorded in `../source_data/SHA256SUMS.txt`.
