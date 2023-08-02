par_ssm_stat_swap <- function(parm, mats = m.swap){
  N_b <- 3 # dimension of factor vector
  N_y <- length(mats) # dimension of observation vector
  if(length(parm) != (17+N_y-1)){return("Not correct length of maturity vector")}
  
  lamb <- exp(parm[13]) # lambda
  
  res <- list(
    "B0" = matrix(parm[1:3], N_b, 1),     # b_0
    "P0" = diag(parm[4:6]^2, N_b, N_b),   # P_0
    "Dm" = matrix(parm[7:9], N_b, 1),     # c
    "Am" = matrix(0, N_y, 1),
    "Fm" = diag(parm[10:12], N_b, N_b),   # Phi
    "Hm" = cbind(rep(1,N_y),              # Lambda
                 (1-exp(-lamb*mats))/(lamb*mats),
                 (1-exp(-lamb*mats))/(lamb*mats) - (exp(-lamb*mats))) %>% as.matrix(),
    "Qm" = diag(parm[14:16]^2, N_b, N_b), # Q
    "Rm" = diag(parm[17:(17+N_y-1)]^2, N_y,N_y)   # H
  )
  return(res)
}

par_ssm_stat_bond <- function(parm, mats = m.bond){
  N_b <- 3 # dimension of factor vector
  N_y <- length(mats) # dimension of observation vector
  if(length(parm) != (17+N_y-1)){return("Not correct length of maturity vector")}
  
  lamb <- exp(parm[13]) # lambda
  
  res <- list(
    "B0" = matrix(parm[1:3], N_b, 1),     # b_0
    "P0" = diag(parm[4:6]^2, N_b, N_b),   # P_0
    "Dm" = matrix(parm[7:9], N_b, 1),     # c
    "Am" = matrix(0, N_y, 1),
    "Fm" = diag(parm[10:12], N_b, N_b),   # Phi
    "Hm" = cbind(rep(1,N_y),              # Lambda
                 (1-exp(-lamb*mats))/(lamb*mats),
                 (1-exp(-lamb*mats))/(lamb*mats) - (exp(-lamb*mats))) %>% as.matrix(),
    "Qm" = diag(parm[14:16]^2, N_b, N_b), # Q
    "Rm" = diag(parm[17:(17+N_y-1)]^2, N_y,N_y)   # H
  )
  return(res)
}

par_ssm_nonstat_swap <- function(parm, mats = m.swap){
  N_b <- 6; N_y <- length(mats)
  if(length(parm) != (17+N_y-1)){return("Not correct length of maturity vector")}
  
  lamb <- exp(parm[13])
  Phi <- diag(parm[10:12], N_b/2, N_b/2)
  I   <- diag(1,3,3)
  
  res <- list(
    "B0" = matrix(rep(parm[1:3], 2), N_b, 1),   # b_0
    "P0" = diag(rep(parm[4:6]^2, 2), N_b, N_b),   # P_0
    "Dm" = matrix(c(parm[7:9], 0,0,0), N_b, 1), # c
    "Am" = matrix(0, N_y, 1),
    "Fm" = rbind(cbind(Phi + I, -Phi),          # Phi
                 cbind(I, matrix(0,3,3))),
    "Hm" = cbind(rep(1,N_y),                    # Lambda
                 (1-exp(-lamb*mats))/(lamb*mats),
                 (1-exp(-lamb*mats))/(lamb*mats) - (exp(-lamb*mats)),
                 matrix(0,N_y,3)),
    "Qm" = diag(c(parm[14:16]^2,0,0,0), N_b, N_b),# Q
    "Rm" = diag(parm[17:(17+N_y-1)]^2, N_y,N_y)           # H
  )
  return(res)
}

par_ssm_nonstat_bond <- function(parm, mats = m.bond){
  N_b <- 6; N_y <- length(mats)
  if(length(parm) != (17+N_y-1)){return("Not correct length of maturity vector")}
  
  lamb <- exp(parm[13])
  Phi <- diag(parm[10:12], N_b/2, N_b/2)
  I   <- diag(1,3,3)
  
  res <- list(
    "B0" = matrix(rep(parm[1:3], 2), N_b, 1),   # b_0
    "P0" = diag(rep(parm[4:6]^2, 2), N_b, N_b),   # P_0
    "Dm" = matrix(c(parm[7:9], 0,0,0), N_b, 1), # c
    "Am" = matrix(0, N_y, 1),
    "Fm" = rbind(cbind(Phi + I, -Phi),          # Phi
                 cbind(I, matrix(0,3,3))),
    "Hm" = cbind(rep(1,N_y),                    # Lambda
                 (1-exp(-lamb*mats))/(lamb*mats),
                 (1-exp(-lamb*mats))/(lamb*mats) - (exp(-lamb*mats)),
                 matrix(0,N_y,3)),
    "Qm" = diag(c(parm[14:16]^2,0,0,0), N_b, N_b),# Q
    "Rm" = diag(parm[17:(17+N_y-1)]^2, N_y,N_y)           # H
  )
  return(res)
}

par_ssm_stat_double <- function(parm, mats.s = m.swap, mats.b = m.bond){
  par <- parm
  mats <- c(mats.s, mats.b)
  # function to form a ssm object from a parameter vector for two rates
  N_b <- 6; N_y <- length(mats)
  # if(length(par) != (39+length(mats)-1)){print("Incorrect dimension of parameter of maturity vector")}
  lambda.swap <- exp(par[31]) # make sure that lambda > 0; transform back to lambda
  lambda.bond <- exp(par[32])
  
  output <- list(
    "B0" = matrix(par[1:6],N_b,1),
    "P0" = diag(par[7:12]^2,N_b,N_b),
    "Dm" = matrix(par[13:18],N_b,1),
    "Am" = matrix(0, N_y, 1),
    "Fm" = matrix(0, N_b,N_b),
    "Hm" = rbind(
      cbind(rep(1,length(mats.s)), # Lambda
            (1-exp(-lambda.swap*mats.s))/(lambda.swap*mats.s),
            (1-exp(-lambda.swap*mats.s))/(lambda.swap*mats.s) - (exp(-lambda.swap*mats.s)),
            matrix(0, length(mats.s),3)) %>% as.matrix(),
      cbind(matrix(0, length(mats.b),3),
            rep(1,length(mats.b)), # Lambda
            (1-exp(-lambda.bond*mats.b))/(lambda.bond*mats.b),
            (1-exp(-lambda.bond*mats.b))/(lambda.bond*mats.b) - (exp(-lambda.bond*mats.b)))) %>% as.matrix(),
    "Qm" = diag(par[33:38]^2,N_b,N_b),
    "Rm" = diag(c(rep(par[39], length(mats.s)), rep(par[40], length(mats.b)))^2,N_y,N_y)
  )
  
  output$Fm[1,1] <- par[19]; output$Fm[1,4] <- par[20]; output$Fm[2,2] <- par[21]
  output$Fm[2,5] <- par[22]; output$Fm[3,3] <- par[23]; output$Fm[3,6] <- par[24]
  output$Fm[4,1] <- par[25]; output$Fm[4,4] <- par[26]; output$Fm[5,2] <- par[27]
  output$Fm[5,5] <- par[28]; output$Fm[6,3] <- par[29]; output$Fm[6,6] <- par[30]
  
  return(output)
}

par_ssm_nonstat_double <- function(par, mats.s = m.swap, mats.b = m.bond){
  mats <- c(mats.s, mats.b)
  # function to form a ssm object from a parameter vector for two rates
  N_b <- 12; N_y <- length(mats)
  # if(length(par) != (39+length(mats)-1)){print("Incorrect dimension of parameter of maturity vector")}
  lambda.swap <- exp(par[31]) # make sure that lambda > 0; transform back to lambda
  lambda.bond <- exp(par[32])
  
  # Construct Phi
  Phi <- matrix(0, 6, 6)
  Phi[1,1] <- par[19]; Phi[1,4] <- par[20]; Phi[2,2] <- par[21]
  Phi[2,5] <- par[22]; Phi[3,3] <- par[23]; Phi[3,6] <- par[24]
  Phi[4,1] <- par[25]; Phi[4,4] <- par[26]; Phi[5,2] <- par[27]
  Phi[5,5] <- par[28]; Phi[6,3] <- par[29]; Phi[6,6] <- par[30]
  # Construct diagonal matrix
  I <- diag(1, 6, 6)
  LargePhi <- rbind(cbind(Phi + I, -Phi), 
                    cbind(I, matrix(0,6,6)))
  
  # Construct Lambda
  Lambda <- rbind(
    cbind(rep(1,length(mats.s)), # Lambda
          (1-exp(-lambda.swap*mats.s))/(lambda.swap*mats.s),
          (1-exp(-lambda.swap*mats.s))/(lambda.swap*mats.s) - (exp(-lambda.swap*mats.s)),
          matrix(0, length(mats.s),9)) %>% as.matrix(),
    cbind(matrix(0, length(mats.b),9),
          rep(1,length(mats.b)), # Lambda
          (1-exp(-lambda.bond*mats.b))/(lambda.bond*mats.b),
          (1-exp(-lambda.bond*mats.b))/(lambda.bond*mats.b) - (exp(-lambda.bond*mats.b))))
  
  output <- list(
    "B0" = matrix(rep(par[1:6], 2),N_b,1),
    "P0" = diag(rep(par[7:12]^2, 2),N_b,N_b),
    "Dm" = matrix(c(par[13:18], rep(0,6)),N_b,1),
    "Am" = matrix(0, N_y, 1),
    "Fm" = LargePhi %>% as.matrix(),
    "Hm" = Lambda %>% as.matrix(),
    "Qm" = diag(par[33:38]^2,N_b,N_b),
    "Rm" = diag(c(rep(par[39], length(mats.s)), rep(par[40], length(mats.b)))^2,N_y,N_y)
    # "Rm" = diag(par[39:(39+N_y-1)]^2,N_y,N_y) # 39:98
  )
  
  return(output)
}


par_lnl_stat_swap <- function(par, mat = m.swap, y = y.swap.train){
  par.ssm <- par_ssm_stat_swap(parm = par, mats = mat)
  KF <- kalman_filter(ssm = par.ssm, 
                      yt = y)
  return(KF$lnl)
}

par_lnl_nonstat_swap <- function(par, mat = m.swap, y=y.swap.train){
  par.ssm <- par_ssm_nonstat_swap(parm = par, mats = mat)
  KF <- kalman_filter(ssm = par.ssm, 
                      yt = y)
  return(KF$lnl)
}

par_lnl_stat_bond <- function(par, mat = m.bond, y = y.bond.train){
  par.ssm <- par_ssm_stat_bond(parm = par, mats = mat)
  KF <- kalman_filter(ssm = par.ssm, 
                      yt = y)
  return(KF$lnl)
}

par_lnl_nonstat_bond <- function(par, mat = m.swap, y=y.bond.train){
  par.ssm <- par_ssm_nonstat_bond(parm = par, mats = mat)
  KF <- kalman_filter(ssm = par.ssm, 
                      yt = y)
  return(KF$lnl)
}

par_lnl_stat_double <- function(par, mat.s = m.swap, mat.b = m.bond, y = rbind(y.swap.train,y.bond.train)){
  par.ssm <- par_ssm_stat_double(parm = par, mats.s = mat.s, mats.b = mat.b)
  KF <- kalman_filter(ssm = par.ssm, 
                      yt = y)
  return(KF$lnl)
}

par_lnl_nonstat_double <- function(par, mat.s = m.swap, mat.b = m.bond, y = rbind(y.swap.train,y.bond.train)){
  par.ssm <- par_ssm_nonstat_double(par = par, mats.s = mat.s, mats.b = mat.b)
  KF <- kalman_filter(ssm = par.ssm, 
                      yt = y)
  return(KF$lnl)
}
