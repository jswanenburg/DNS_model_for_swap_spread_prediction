# import all parameter data frames
if(settings %in% c(2,5)){
  # p.swap.nonstat.df <- read.csv(str_c(wd.result, "p_swap_nonstat.csv")) %>% dplyr::select(-X)
  # p.bond.nonstat.df <- read.csv(str_c(wd.result, "p_bond_nonstat.csv")) %>% dplyr::select(-X)
  p.double.nonstat.df <- read.csv(str_c(wd.result, "p_double_nonstat.csv")) %>% dplyr::select(-X)
} else if(settings == 4){
  p.swap.stat.df <- read.csv(str_c(wd.result, "p_swap_stat.csv")) %>% dplyr::select(-X)
  p.bond.stat.df <- read.csv(str_c(wd.result, "p_bond_stat.csv")) %>% dplyr::select(-X)
  p.double.stat.df <- read.csv(str_c(wd.result, "p_double_stat.csv")) %>% dplyr::select(-X)
}

max.iterations <- 10 # number of iterations per loop
optim.method <- "BFGS"
skip.setting <- FALSE


##### Double stationary 
if(!exists("p.double.stat.df")){
  cat("\nStatDNS model for double data is not estimated for settings", settings)
} else if(p.double.stat.df$settings[1] != 1){
  cat("\nMaximum likelihood for stationary bond rate is already converged.")
} else {
  last.colname <- tail(colnames(p.double.stat.df), 1)
  cat("\n\n=====================================================================",
      "\nStart MLE estimation for double rate with stationary factor dynamics",
      "\n<Running", max.iterations, "iterations per loop for the", 
      optim.method, "optimization algorithm>",
      "\nStarted at ", as.character(Sys.time()),
      "with", last.colname, "parameters.",
      "\n=====================================================================")
  
  # set looping parameters for MLE estimation
  p.MLE.loop <- p.double.stat.df[[last.colname]]
  convergence.indicator <- p.double.stat.df$settings[1]
  loop <- str_split(last.colname, "_")[[1]][3] %>% as.numeric() + 1
  if(is.na(loop)){loop <- 1}
  
  while(convergence.indicator == 1){
    cat("\n\nLoop:", loop, "- started at", as.character(Sys.time()),
        "\nLambda in MLE swap: ", signif(exp(p.MLE.loop[31]), 5),
        "\nLambda in MLE bond: ", signif(exp(p.MLE.loop[32]), 5))
    
    # Run MLE optimization for 'maxit' iterations
    start.time <- proc.time()
    MLE.result <- optim(par = p.MLE.loop,
                        method = optim.method,
                        fn = par_lnl_stat_double,
                        control=list("fnscale" = -1,   # for maximization (default is minimization)
                                     "maxit" = max.iterations))
    stop.time <- proc.time()
    
    # add new MLE vector to data frame
    p.double.stat.df$MLE <- MLE.result$par
    colnames(p.double.stat.df)[which(colnames(p.double.stat.df) == "MLE")] <- str_c("MLE_loop_", loop)
    
    p.double.stat.df$settings[2] <- p.double.stat.df$settings[2] + (stop.time-start.time)[3] # total runtime
    p.double.stat.df$settings[1] <- MLE.result$convergence # convergence indicator
    
    convergence.indicator <- MLE.result$convergence
    
    # write new results to csv
    write.csv(p.double.stat.df, str_c(wd.result, "p_double_stat.csv"))
    
    # print loop info
    cat("\nMLE loop", loop, "estimated in ", round((stop.time-start.time)[3]/60, 3), "minutes",
        "\n<MLE algorithm runs for settings:", settings, ">")
    if(convergence.indicator == 1){cat("\nMLE algorithm for stationary bond rate has not converged yet.")}
    
    # make all parameters ready for next loop
    loop <- loop + 1
    p.MLE.loop <- MLE.result$par
  }
  
}



##### Double nonstationary 
if(!exists("p.double.nonstat.df")){
  cat("\nNontatDNS model for double data is not estimated for settings", settings)
} else if(p.double.nonstat.df$settings[1] != 1){
  cat("\nMaximum likelihood for stationary bond rate is already converged.")
} else {
  last.colname <- tail(colnames(p.double.nonstat.df), 1)
  cat("\n\n=====================================================================",
      "\nStart MLE estimation for double rate with nonstationary factor dynamics",
      "\n<Running", max.iterations, "iterations per loop for the", 
      optim.method, "optimization algorithm>",
      "\nStarted at ", as.character(Sys.time()),
      "with", last.colname, "parameters.",
      "\n=====================================================================")
  
  # set looping parameters for MLE estimation
  p.MLE.loop <- p.double.nonstat.df[[last.colname]]
  convergence.indicator <- p.double.nonstat.df$settings[1]
  loop <- str_split(last.colname, "_")[[1]][3] %>% as.numeric() + 1
  if(is.na(loop)){loop <- 1}
  
  while(convergence.indicator == 1){
    cat("\n\nLoop:", loop, "- started at", as.character(Sys.time()),
        "\nLambda in MLE swap: ", signif(exp(p.MLE.loop[31]), 5),
        "\nLambda in MLE bond: ", signif(exp(p.MLE.loop[32]), 5))
    
    # Run MLE optimization for 'maxit' iterations
    start.time <- proc.time()
    MLE.result <- optim(par = p.MLE.loop,
                        method = optim.method,
                        fn = par_lnl_nonstat_double,
                        control=list("fnscale" = -1,   # for maximization (default is minimization)
                                     "maxit" = max.iterations))
    stop.time <- proc.time()
    
    # add new MLE vector to data frame
    p.double.nonstat.df$MLE <- MLE.result$par
    colnames(p.double.nonstat.df)[which(colnames(p.double.nonstat.df) == "MLE")] <- str_c("MLE_loop_", loop)
    
    p.double.nonstat.df$settings[2] <- p.double.nonstat.df$settings[2] + (stop.time-start.time)[3] # total runtime
    p.double.nonstat.df$settings[1] <- MLE.result$convergence # convergence indicator
    
    convergence.indicator <- MLE.result$convergence
    
    # write new results to csv
    write.csv(p.double.nonstat.df, str_c(wd.result, "p_double_nonstat.csv"))
    
    # print loop info
    cat("\nMLE loop", loop, "estimated in ", round((stop.time-start.time)[3]/60, 3), "minutes",
        "\n<MLE algorithm runs for settings:", settings, ">")
    if(convergence.indicator == 1){cat("\nMLE algorithm for nonstationary double rate has not converged yet.")}
    
    # make all parameters ready for next loop
    loop <- loop + 1
    p.MLE.loop <- MLE.result$par
    
    # if(exp(p.MLE.loop[31]) > 10 | exp(p.MLE.loop[31]) < 0.01 |
    #    exp(p.MLE.loop[32]) > 10 | exp(p.MLE.loop[32]) < 0.01){
    #   convergence.indicator <- 2
    #   cat("\nLambda's outside range. Skip this settings.")
    #   skip.setting <- TRUE
    # }
    # if(skip.setting){
    #   convergence.indicator <- 2
    # }
  }
  
}

##### Swap stationary 
if(!exists("p.swap.stat.df")){
  cat("\nStatDNS model for swap data is not estimated for settings", settings)
} else if(p.swap.stat.df$settings[1] != 1){
  cat("\nMaximum likelihood for stationary swap rate is already converged.")
} else {
  last.colname <- tail(colnames(p.swap.stat.df), 1)
  cat("\n\n=====================================================================",
      "\nStart MLE estimation for single swap rate with stationary factor dynamics",
      "\n<Running", max.iterations, "iterations per loop for the", 
      optim.method, "optimization algorithm>",
      "\nStarted at ", as.character(Sys.time()),
      "with", last.colname, "parameters.",
      "\n=====================================================================")
  
  # set looping parameters for MLE estimation
  p.MLE.loop <- p.swap.stat.df[[last.colname]]
  convergence.indicator <- p.swap.stat.df$settings[1]
  loop <- str_split(last.colname, "_")[[1]][3] %>% as.numeric() + 1
  if(is.na(loop)){loop <- 1}
  
  while(convergence.indicator == 1){
    cat("\n\nLoop:", loop, "- started at", as.character(Sys.time()),
        "\nLambda in MLE: ", signif(exp(p.MLE.loop[13]), 5))
    
    # Run MLE optimization for 'maxit' iterations
    start.time <- proc.time()
    MLE.result <- optim(par = p.MLE.loop,
                        method = optim.method,
                        fn = par_lnl_stat_swap,
                        control=list("fnscale" = -1,   # for maximization (default is minimization)
                                     "maxit" = max.iterations))
    stop.time <- proc.time()
    
    # add new MLE vector to data frame
    p.swap.stat.df$MLE <- MLE.result$par
    colnames(p.swap.stat.df)[which(colnames(p.swap.stat.df) == "MLE")] <- str_c("MLE_loop_", loop)
    
    p.swap.stat.df$settings[2] <- p.swap.stat.df$settings[2] + (stop.time-start.time)[3] # total runtime
    p.swap.stat.df$settings[1] <- MLE.result$convergence # convergence indicator
    
    convergence.indicator <- MLE.result$convergence
    
    # write new results to csv
    write.csv(p.swap.stat.df, str_c(wd.result, "p_swap_stat.csv"))
    
    # print loop info
    cat("\nMLE loop", loop, "estimated in ", round((stop.time-start.time)[3]/60, 3), "minutes",
        "\n<MLE algorithm runs for settings:", settings, ">")
    if(convergence.indicator == 1){
      cat("\nMLE algorithm for stationary swap rate has not converged yet.")
    } else{
      cat("\nMLE algorithm for stationary swap rate has converged.")
    }
    
    # make all parameters ready for next loop
    loop <- loop + 1
    p.MLE.loop <- MLE.result$par
  }
  
}

##### Swap nonstationary 
if(!exists("p.swap.nonstat.df")){
  cat("\nNontatDNS model for swap data is not estimated for settings", settings)
} else if(p.swap.nonstat.df$settings[1] != 1){
  cat("\nMaximum likelihood for nonstationary swap rate is already converged.")
} else {
  last.colname <- tail(colnames(p.swap.nonstat.df), 1)
  cat("\n\n=====================================================================",
      "\nStart MLE estimation for single swap rate with nonstationary factor dynamics",
      "\n<Running", max.iterations, "iterations per loop for the", 
      optim.method, "optimization algorithm>",
      "\nStarted at ", as.character(Sys.time()),
      "with", last.colname, "parameters.",
      "\n=====================================================================")
  
  # set looping parameters for MLE estimation
  p.MLE.loop <- p.swap.nonstat.df[[last.colname]]
  convergence.indicator <- p.swap.nonstat.df$settings[1]
  loop <- str_split(last.colname, "_")[[1]][3] %>% as.numeric() + 1
  if(is.na(loop)){loop <- 1}
  
  while(convergence.indicator == 1){
    cat("\n\nLoop:", loop, "- started at", as.character(Sys.time()),
        "\nLambda in MLE: ", signif(exp(p.MLE.loop[13]), 5))
    
    
    # Run MLE optimization for 'maxit' iterations
    start.time <- proc.time()
    MLE.result <- optim(par = p.MLE.loop,
                        method = optim.method,
                        fn = par_lnl_nonstat_swap,
                        control=list("fnscale" = -1,   # for maximization (default is minimization)
                                     "maxit" = max.iterations))
    stop.time <- proc.time()
    
    # add new MLE vector to data frame
    p.swap.nonstat.df$MLE <- MLE.result$par
    colnames(p.swap.nonstat.df)[which(colnames(p.swap.nonstat.df) == "MLE")] <- str_c("MLE_loop_", loop)
    
    p.swap.nonstat.df$settings[2] <- p.swap.nonstat.df$settings[2] + (stop.time-start.time)[3] # total runtime
    p.swap.nonstat.df$settings[1] <- MLE.result$convergence # convergence indicator
    
    convergence.indicator <- MLE.result$convergence
    
    # write new results to csv
    write.csv(p.swap.nonstat.df, str_c(wd.result, "p_swap_nonstat.csv"))
    
    # print loop info
    cat("\nMLE loop", loop, "estimated in ", round((stop.time-start.time)[3]/60, 3), "minutes",
        "\n<MLE algorithm runs for settings:", settings, ">")
    if(convergence.indicator == 1){cat("\nMLE algorithm for nonstationary swap rate has not converged yet.")}
    
    # make all parameters ready for next loop
    loop <- loop + 1
    p.MLE.loop <- MLE.result$par
    
    if(exp(p.MLE.loop[13]) > 10 | exp(p.MLE.loop[13]) < 0.01){
      convergence.indicator <- 2
      cat("\nLambda's outside range. Skip this settings.")
      skip.setting <- TRUE
    }
    if(skip.setting){
      convergence.indicator <- 2
    }
  }
  
}

##### Bond stationary 
if(!exists("p.bond.stat.df")){
  cat("\nStatDNS model for bond data is not estimated for settings", settings)
} else if(p.bond.stat.df$settings[1] != 1){
  cat("\nMaximum likelihood for stationary bond rate is already converged.")
} else {
  last.colname <- tail(colnames(p.bond.stat.df), 1)
  cat("\n\n=====================================================================",
      "\nStart MLE estimation for single bond rate with stationary factor dynamics",
      "\n<Running", max.iterations, "iterations per loop for the", 
      optim.method, "optimization algorithm>",
      "\nStarted at ", as.character(Sys.time()),
      "with", last.colname, "parameters.",
      "\n=====================================================================")
  
  # set looping parameters for MLE estimation
  p.MLE.loop <- p.bond.stat.df[[last.colname]]
  convergence.indicator <- p.bond.stat.df$settings[1]
  loop <- str_split(last.colname, "_")[[1]][3] %>% as.numeric() + 1
  if(is.na(loop)){loop <- 1}
  
  while(convergence.indicator == 1){
    cat("\n\nLoop:", loop, "- started at", as.character(Sys.time()),
        "\nLambda in MLE: ", signif(exp(p.MLE.loop[13]), 5))
    
    # Run MLE optimization for 'maxit' iterations
    start.time <- proc.time()
    MLE.result <- optim(par = p.MLE.loop,
                        method = optim.method,
                        fn = par_lnl_stat_bond,
                        control=list("fnscale" = -1,   # for maximization (default is minimization)
                                     "maxit" = max.iterations))
    stop.time <- proc.time()
    
    # add new MLE vector to data frame
    p.bond.stat.df$MLE <- MLE.result$par
    colnames(p.bond.stat.df)[which(colnames(p.bond.stat.df) == "MLE")] <- str_c("MLE_loop_", loop)
    
    p.bond.stat.df$settings[2] <- p.bond.stat.df$settings[2] + (stop.time-start.time)[3] # total runtime
    p.bond.stat.df$settings[1] <- MLE.result$convergence # convergence indicator
    
    convergence.indicator <- MLE.result$convergence
    
    # write new results to csv
    write.csv(p.bond.stat.df, str_c(wd.result, "p_bond_stat.csv"))
    
    # print loop info
    cat("\nMLE loop", loop, "estimated in ", round((stop.time-start.time)[3]/60, 3), "minutes",
        "\n<MLE algorithm runs for settings:", settings, ">")
    if(convergence.indicator == 1){cat("\nMLE algorithm for stationary bond rate has not converged yet.")}
    
    # make all parameters ready for next loop
    loop <- loop + 1
    p.MLE.loop <- MLE.result$par
  }
  
}


##### Bond nonstationary 
if(!exists("p.bond.nonstat.df")){
  cat("\nNonstatDNS model for bond data is not estimated for settings", settings)
} else if(p.bond.nonstat.df$settings[1] != 1){
  cat("\nMaximum likelihood for nonstationary bond rate is already converged.")
} else {
  last.colname <- tail(colnames(p.bond.nonstat.df), 1)
  cat("\n\n=====================================================================",
      "\nStart MLE estimation for single swap rate with nonstationary factor dynamics",
      "\n<Running", max.iterations, "iterations per loop for the", 
      optim.method, "optimization algorithm>",
      "\nStarted at ", as.character(Sys.time()),
      "with", last.colname, "parameters.",
      "\n=====================================================================")
  
  # set looping parameters for MLE estimation
  p.MLE.loop <- p.bond.nonstat.df[[last.colname]]
  convergence.indicator <- p.bond.nonstat.df$settings[1]
  loop <- str_split(last.colname, "_")[[1]][3] %>% as.numeric() + 1
  if(is.na(loop)){loop <- 1}
  
  while(convergence.indicator == 1){
    cat("\n\nLoop:", loop, "- started at", as.character(Sys.time()),
        "\nLambda in MLE: ", signif(exp(p.MLE.loop[13]), 5))
    
    
    # Run MLE optimization for 'maxit' iterations
    start.time <- proc.time()
    MLE.result <- optim(par = p.MLE.loop,
                        method = optim.method,
                        fn = par_lnl_nonstat_bond,
                        control=list("fnscale" = -1,   # for maximization (default is minimization)
                                     "maxit" = max.iterations))
    stop.time <- proc.time()
    
    # add new MLE vector to data frame
    p.bond.nonstat.df$MLE <- MLE.result$par
    colnames(p.bond.nonstat.df)[which(colnames(p.bond.nonstat.df) == "MLE")] <- str_c("MLE_loop_", loop)
    
    p.bond.nonstat.df$settings[2] <- p.bond.nonstat.df$settings[2] + (stop.time-start.time)[3] # total runtime
    p.bond.nonstat.df$settings[1] <- MLE.result$convergence # convergence indicator
    
    convergence.indicator <- MLE.result$convergence
    
    # write new results to csv
    write.csv(p.bond.nonstat.df, str_c(wd.result, "p_bond_nonstat.csv"))
    
    # print loop info
    cat("\nMLE loop", loop, "estimated in ", round((stop.time-start.time)[3]/60, 3), "minutes",
        "\n<MLE algorithm runs for settings:", settings, ">")
    if(convergence.indicator == 1){cat("\nMLE algorithm for nonstationary bond rate has not converged yet.")}
    
    # make all parameters ready for next loop
    loop <- loop + 1
    p.MLE.loop <- MLE.result$par
    
    if(exp(p.MLE.loop[13]) > 10 | exp(p.MLE.loop[13]) < 0.01){
      convergence.indicator <- 2
      cat("\nLambda's outside range. Skip this settings.")
      skip.setting <- TRUE
    }
    if(skip.setting){
      convergence.indicator <- 2
    }
  }
  
}
