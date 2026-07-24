# SPP1+ TAM–myCAF immunosuppressive niche in intrahepatic cholangiocarcinoma

Analysis code for the integrative single-cell, spatial and multi-cohort survival
re-analysis of the SPP1+/TREM2+ tumour-associated macrophage (TAM) – myofibroblastic
cancer-associated fibroblast (myCAF) niche in intrahepatic cholangiocarcinoma (iCCA).

## Overview
The pipeline (Python 3.10) reproduces, in execution order:
1. Single-cell processing of GSE138709 (QC, normalization, Harmony integration, Leiden, annotation)
2. Myeloid and stromal sub-clustering and subtype annotation
3. SPP1 ligand–receptor interaction (label-permutation test)
4. Spatial co-localization in GSE316402 Visium (per-spot scoring; neighbourhood bivariate Moran's I)
5. Bulk prognostic scoring and survival (TCGA-CHOL, GSE107943)
6. E-MTAB-6389 (HTA-2.0 probe→gene mapping) and three-cohort random-effects meta-analysis

## Repository structure
```

  00_figstyle.py               plotting style (Arial) + UMAP legend helper
  01a_load_per_sample.py       load each GSE138709 sample to sparse AnnData
  01b_merge_qc.py              merge samples + QC metrics
  01c_cluster_annotate.py      filter, normalize, HVG, PCA, Harmony, Leiden, UMAP, annotation
  02_subcluster_myeloid_caf.py myeloid + CAF sub-clustering and annotation
  03_interaction_permutation.py SPP1 ligand–receptor permutation test
  04a_spatial_score.py         Visium per-spot signature scoring
  04b_spatial_moran.py         neighbourhood bivariate Moran's I
  05_bulk_survival.py          bulk scoring, Cox/KM, meta-analysis helper
  06_emtab6389_meta.py         E-MTAB-6389 + 3-cohort meta-analysis + forest plot
environment.yml / requirements.txt   software environment
```

## Data availability
All datasets are public and are NOT redistributed here:
| Dataset | Source |
|---|---|
| GSE138709 (scRNA-seq) | NCBI GEO |
| GSE316402 (spatial Visium) | NCBI GEO |
| GSE107943 (bulk RNA-seq + survival) | NCBI GEO |
| E-MTAB-6389 (microarray + survival) | EBI ArrayExpress/BioStudies |
| TCGA-CHOL (bulk + survival) | UCSC Xena GDC hub |

## Reproduce
```bash
conda env create -f environment.yml
conda activate icca-spp1-niche
# download the public datasets into a ./data folder, then run the modules in order
python 01a_load_per_sample.py 0   # ... per sample
python 01b_merge_qc.py
python 01c_cluster_annotate.py
# ...
```

## Citation
If you use this code, please cite the associated manuscript (see CITATION.cff) and the
original data-generating studies.

## License
MIT License (see LICENSE).
