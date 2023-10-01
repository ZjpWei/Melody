# Simulation 

**1. Generate simulation data.**
```console
Rscript Prepare_data.R arg1 arg2 arg3 arg4 arg5
```
   Arguments:

   `arg1` (`s`): replicate number; from `1` to `100`.
   
   `arg2` (`Ka`): signature sparsity; {`20`, `30`, `40`(default), `50`, `60`, `70`, `80`}.
   
   `arg3` (`pos.pt`): signature effect direction; {`0.5`, `0.6`, `0.7`(default), `0.8`, `0.9`, `1`}.
   
   `arg4` (`effect.sz`): signature effect size; {`0.5`, `1`, `1.5`, `2`(default), `2.5`, `3`}.
   
   `arg5` (`mu`): case/control sequence depth unevenness; {`0`(default), `0.25`, `0.5`, `0.75`, `1`}
   
**2. Perform meta-analysis on original data using `miMeta` packages.**
```console
Rscript Melody.R arg1 arg2 arg3
```

**3. Perform other comparison methods using original data.**
```console
Rscript Compared_methods_original.R arg1 arg2 arg3
```

**4. Perform other comparison methods using batch-corrected data.**
```console
Rscript Compared_methods_batch_corrected.R arg1 arg2 arg3
```

**5. Caluclate AUPRC**
```console
Rscript AUPRC.R arg1 arg2 arg3
```

   Arguments:
  
   `arg1` (`s`): replicate number; from `1` to `100`.
   
   `arg2` (`scenario`): simulation scenario; {`Sig_number`, `Sig_effdir`, `Sig_effsz`, `Sig_depth`}.
   
   * signature sparsity: "Sig_number"
    
   * signature effect direction: "Sig_effdir"
    
   * signature effect size: "Sig_effsz"
    
   * case/control sequence depth unevenness: "Sig_depth"
    
   `arg3` (`loc`): factor of scenario;
   
   * for signature sparsity; {`20`, `30`, `40`, `50`, `60`, `70`, `80`}.
     
   * for signature effect direction; {`0.5`, `0.6`, `0.7`, `0.8`, `0.9`, `1`}.
    
   * for signature effect size; {`0.5`, `1`, `1.5`, `2`, `2.5`, `3`}.
    
   * for case/control sequence depth unevenness; {`0`, `0.25`, `0.5`, `0.75`, `1`}
