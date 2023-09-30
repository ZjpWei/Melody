# Prediction

This file is the workflow for doing Radom-split prediction and Leave-one-study-out prediction.

## Random-split prediction (RD)
Workflow:
1. Perform meta-analysis on original data using `miMeta` packages.
```console
Rscript Melody_abess_RD.R arg1
```

2. Perform other comparison methods using original data.
```console
Rscript train_models_RD.R arg1
```

3. Perform other comparison methods using batch-corrected data.
```console
Rscript train_batch_models_RD.R arg1
```

4. Summerize the top-rank taxa for all methods
```console
Rscript Get_index_RD.R
```

5. Random-split prediction for top-rank taxa
```console
Rscript Prediction_RD.R arg1 arg2 arg3
```
Arguments:
1. `arg1` (`s`): replicate number; from `1` to `100`.
2. `arg2` (`methodology`): method; `Aldex2`, `ANCOM`, `BW`, `clrlasso`, `Melody`.
3. `arg3` (`data.type`): data type; `Original`, `Batch-corrected`.

## Leave-one-study-out prediction (LOSO)
Workflow:
1. Perform meta-analysis on original data using `miMeta` packages.
```console
Rscript Melody_abess_LOSO.R arg1 arg2
```

2. Perform other comparison methods using original data.
```console
Rscript train_models_LOSO.R arg1 arg2
```

3. Perform other comparison methods using batch-corrected data.
```console
Rscript train_batch_models_LOSO.R arg1 arg2
```

4. Summerize the top-rank taxa for all methods
```console
Rscript Get_index_LOSO.R
```

5. Random-split prediction for top-rank taxa
```console
Rscript Prediction_LOSO.R arg1 arg2 arg3 arg4
```
Arguments:
1. `arg1` (`s`): replicate number; from `1` to `100`.
2. `arg2` (`ss`): study number for leaving out.
3. `arg3` (`methodology`): method; `Aldex2`, `ANCOM`, `BW`, `clrlasso`, `Melody`.
4. `arg4` (`data.type`): data type; `Original`, `Batch-corrected`.
