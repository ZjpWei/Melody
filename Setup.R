# =============================================== #
#                   Files setup                   #
# =============================================== #

  # Check directory and create directory
  # CRC_Real directory
  loc.vec <- c("./CRC_meta_analysis/CRC_all_K401",
               "./CRC_meta_analysis/CRC_all_K849",
               "./CRC_meta_analysis/CRC_all_order",
               "./CRC_meta_analysis/CRC_loso_K401",
               "./CRC_meta_analysis/CRC_loso_K849",
               "./CRC_meta_analysis/CRC_loso_order",
               "./CRC_meta_analysis/Prediction")

  for(loc in loc.vec){
    if(!dir.exists(loc)){
      dir.create(loc)
    }
  }

  # Prediction directory
  loc.vec <- c("./CRC_meta_analysis/Prediction/LOSO")

  for(loc in loc.vec){
    if(!dir.exists(loc)){
      dir.create(loc)
    }
  }

  # Simulation directory
  loc.vec <- c("./Simulation/large",
               "./Simulation/small",
               "./Simulation/large/Sig_number",
               "./Simulation/large/Sig_effdir",
               "./Simulation/large/Sig_effsz",
               "./Simulation/large/Sig_depth",
               "./Simulation/small/Sig_number",
               "./Simulation/small/Sig_effdir",
               "./Simulation/small/Sig_effsz",
               "./Simulation/small/Sig_depth")

  for(loc in loc.vec){
    if(!dir.exists(loc)){
      dir.create(loc)
    }
  }

  # Figures directory
  if(!dir.exists("./figures")){
    dir.create("./figures")
  }
