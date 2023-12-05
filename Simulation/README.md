# Simulation 

**1. Generate simulation data.**
```console
Rscript Prepare_data.R arg1 arg2 arg3 arg4 arg5 arg6
```
   Arguments:

   `arg1` (`s`): replicate number; from `1` to `100`.
   
   `arg2` (`Ka`): signature sparsity; {`20`, `30`, `40`(default), `50`, `60`, `70`, `80`}.
   
   `arg3` (`pos.pt`): signature effect direction; {`0.5`, `0.6`, `0.7`(default), `0.8`, `0.9`, `1`}.
   
   `arg4` (`effect.sz`): signature effect size; 
   
            - large sample size set:{`1`, `1.5`, `2`(default), `2.5`, `3`}.
            
            - small sample size set:{`4`, `4.5`, `5`(default), `5.5`, `6`}.
   
   `arg5` (`mu`): case/control sequence depth unevenness; {`0`(default), `0.25`, `0.5`, `0.75`, `1`}

   `arg6` (`set`): large/samll sample size set; {large, small}
   
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
Rscript Compared_methods_batch_corrected.R arg1 arg2 arg3 arg4
```
   Arguments:
  
   `arg1` (`s`): replicate number; from `1` to `100`.

   `arg2` (`set`): large/samll sample size set; {large, small}.
   
   `arg3` (`scenario`): simulation scenario; {`Sig_number`, `Sig_effdir`, `Sig_effsz`, `Sig_depth`}.
   
   * signature sparsity: "Sig_number"
    
   * signature effect direction: "Sig_effdir"
    
   * signature effect size: "Sig_effsz"
    
   * case/control sequence depth unevenness: "Sig_depth"
    
   `arg4` (`loc`): factor of scenario;
   
   * for signature sparsity; {`d20`, `d30`, `d40`, `d50`, `d60`, `d70`, `d80`}.
     
   * for signature effect direction; {`p50`, `p60`, `p70`, `p80`, `p90`, `p100`}.
    
   * for signature effect size;

            - large sample size scenario:{`sig1`, `sig1.5`, `sig2`, `sig2.5`, `sig3`}.
            
            - small sample size scenario:{`sig4`, `sig4.5`, `sig5`, `sig5.5`, `sig6`}.
    
   * for case/control sequence depth unevenness; {`seq0`, `seq0.25`, `seq0.5`, `seq0.75`, `seq1`}

**5. Caluclate AUPRC**
```console
Rscript AUPRC.R arg1 arg2
```

   Arguments:

   `arg1` (`set`): large/samll sample size set; {large, small}.
   
   `arg2` (`scenario`): simulation scenario; {`Sig_number`, `Sig_effdir`, `Sig_effsz`, `Sig_depth`}.
   
   * signature sparsity: "Sig_number"
    
   * signature effect direction: "Sig_effdir"
    
   * signature effect size: "Sig_effsz"
    
   * case/control sequence depth unevenness: "Sig_depth"
    
   `arg3` (`loc`): factor of scenario;
   
   * for signature sparsity; {`d20`, `d30`, `d40`, `d50`, `d60`, `d70`, `d80`}.
     
   * for signature effect direction; {`p50`, `p60`, `p70`, `p80`, `p90`, `p100`}.
    
   * for signature effect size; {`sig1`, `sig1.5`, `sig2`, `sig2.5`, `sig3`}.
     
            - large sample size scenario:{`sig1`, `sig1.5`, `sig2`, `sig2.5`, `sig3`}.
            
            - small sample size scenario:{`sig4`, `sig4.5`, `sig5`, `sig5.5`, `sig6`}.
    
   * for case/control sequence depth unevenness; {`seq0`, `seq0.25`, `seq0.5`, `seq0.75`, `seq1`}
