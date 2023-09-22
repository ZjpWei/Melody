# CRC Real data analysis

## Data analysis

scripts:

1. `Melody_abess.R` performs meta-analysis on CRC original data by `miMeta` package

2. `train_model.R$` analyzes CRC original data by other comparison methods.

3. `train_model_batch.R$` analyzes CRC batch-corrected data by other comparison methods.

arguments:

1. scenario `s`: Data scenario.

2. file `tag`: Folder name for analysis.

examples: 

1. Run `Melody_abess.R` for data in folder `CRC_all_order` for `all` data
```console
Rscript Melody_abess.R all CRC_all_order
```

1. Run `train_model.R` for data in folder `CRC_loso_order` by leaving study `1` out.
```console
Rscript train_model.R 1 CRC_loso_order
```

## Sensitivity: meta-analysis with difference reference taxon

1. `sensitivity.R` performs meta-analysis on CRC original data by `miMeta` package with 5 different reference taxon 
