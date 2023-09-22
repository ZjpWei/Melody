  ##############################################
#                                                #
#             simulation file setup              #
#                                                #
  ##############################################

  # Check directory and create directory
  # CRC_Real directory
  loc.vec <- c("./CRC_Real/CRC_all_K401",
               "./CRC_Real/CRC_all_K849",
               "./CRC_Real/CRC_all_order",
               "./CRC_Real/CRC_loso_K401",
               "./CRC_Real/CRC_loso_K849",
               "./CRC_Real/CRC_loso_order")
  
  for(loc in loc.vec){
    if(!dir.exists(loc)){
      dir.create(loc)
    }
  }
  
  # Prediction directory
  loc.vec <- c("./Prediction/LOSO",
               "./Prediction/Random_split")
  
  for(loc in loc.vec){
    if(!dir.exists(loc)){
      dir.create(loc)
    }
  }
  
  # Simulation directory
  loc.vec <- c("./Simulation/By_taxnum",
               "./Simulation/By_pos_prop",
               "./Simulation/By_abd_prop",
               "./Simulation/By_effsz",
               "./Simulation/By_seqdepth")
  
  for(loc in loc.vec){
    if(!dir.exists(loc)){
      dir.create(loc)
    }
  }
  
