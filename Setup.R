  ##############################################
#                                                #
#                   Files setup                  #
#                                                #
  ##############################################

  # Check directory and create directory
  # CRC_Real directory
  loc.vec <- c("./CRC_meta_analysis/CRC_all_K401",
               "./CRC_meta_analysis/CRC_all_K849",
               "./CRC_meta_analysis/CRC_all_order",
               "./CRC_meta_analysis/CRC_loso_K401",
               "./CRC_meta_analysis/CRC_loso_K849",
               "./CRC_meta_analysis/CRC_loso_order")

  for(loc in loc.vec){
    if(!dir.exists(loc)){
      dir.create(loc)
    }
  }

  # Prediction directory
  loc.vec <- c("./CRC_meta_analysis/Prediction/LOSO",
               "./CRC_meta_analysis/Prediction/Random_split")

  for(loc in loc.vec){
    if(!dir.exists(loc)){
      dir.create(loc)
    }
  }

  # Simulation directory
  loc.vec <- c("./Simulation/Sig_number",
               "./Simulation/Sig_effdir",
               "./Simulation/Sig_effsz",
               "./Simulation/Sig_depth")

  for(loc in loc.vec){
    if(!dir.exists(loc)){
      dir.create(loc)
    }
  }

  # Figures directory
  if(!dir.exists("./figures")){
    dir.create("./figures")
  }
