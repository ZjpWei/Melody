# Simulation 

## Workflow:
1. Generate simulation data.
```console
Rscript prepare_data.R arg1 arg2 arg3 arg4 arg5 arg6
```
   Arguments:

   `arg1` (`s`): replicate number; from `1` to `100`.
   
   `arg2` (`Ka`): signature sparsity; {`20`, `30`, `40`(default), `50`, `60`, `70`, `80`}.
   
   `arg3` (`pos.pt`): signature effect direction; {`0.5`, `0.6`, `0.7`(default), `0.8`, `0.9`, `1`}.
   
   `arg4` (`abd.pt`): signarure prevalence; {`0.2`, `0.3`, `0.4`, `0.5`(default), `0.6`, `0.7`, `0.8`}.
   
   `arg5` (`effect.sz`): signature effect size; {`0.5`, `1`, `1.5`, `2`(default), `2.5`, `3`}.
   
   `arg6` (`mu`): case/control sequence depth unevenness; {`1`(default), `1.25`, `1.5`, `1.75`, `2`}
   
2. Perform meta-analysis on original data using `miMeta` packages.
```console
Rscript Melody_abess.R arg1 arg2 arg3
```

3. Perform other comparison methods using original data.
```console
Rscript train_models.R arg1 arg2 arg3
```

4. Perform other comparison methods using batch-corrected data.
```console
Rscript train_batch_models.R arg1 arg2 arg3
```

5. Caluclate AUPRC
```console
Rscript AUPRC.R arg1 arg2 arg3
```
   Arguments:

   `arg1` (`s`): replicate number; from `1` to `100`.
   
   `arg2` (`scenario`): simulation scenario; {`By_taxnum`, `By_pos_prop`, `By_abd_prop`, `By_effsz`, `By_seqdepth`}.
   
   signature sparsity: "By_taxnum"
    
       signature effect direction: "By_pos_prop"
    
       signature prevalence: "By_abd_prop"
    
       signature effect size: "By_effsz"
    
       case/control sequence depth unevenness: "By_seqdepth"
    
   `arg3` (`loc`): factor of scenario;
   
        signature sparsity; {`20`, `30`, `40`, `50`, `60`, `70`, `80`}.
     
        signature effect direction; {`0.5`, `0.6`, `0.7`, `0.8`, `0.9`, `1`}.
    
        signarure prevalence; {`0.2`, `0.3`, `0.4`, `0.5`, `0.6`, `0.7`, `0.8`}.
    
        signature effect size; {`0.5`, `1`, `1.5`, `2`, `2.5`, `3`}.
    
        case/control sequence depth unevenness; {`1`, `1.25`, `1.5`, `1.75`, `2`}
