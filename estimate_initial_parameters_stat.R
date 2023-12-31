# Initial parameters if the factor dynamics are stationary
##### Create all empty parameter vectors 
p.swap.stat <- rep(NA_real_, 17 + length(m.swap) - 1)
p.swap.stat[13] <- log(lambda.swap)
p.bond.stat <- rep(NA_real_, 17 + length(m.bond) - 1)
p.bond.stat[13] <- log(lambda.bond)
p.double.stat <- rep(NA_real_, 40)
p.double.stat[31:32] <- c(log(lambda.swap), log(lambda.bond))

##### SWAP ====================================================================
# Estimate the H matrix entries (same for stationary and nonstationary models)
p.swap.stat[17:(17 + length(m.swap) - 1)] <- apply(residuals.OLS.swap, 1, sd)
p.double.stat[39] <- mean(apply(residuals.OLS.swap, 1, sd))

# Estimate b_0 and P_0
p.swap.stat[1:3] <- b.OLS.swap[,1] # b_0 is equal to first estimated factor
p.swap.stat[4:6] <- c(sd(b.OLS.swap[1,]), sd(b.OLS.swap[2,]), sd(b.OLS.swap[3,])) # P_0
p.double.stat[1:3] <- b.OLS.swap[,1]
p.double.stat[7:9] <- c(sd(b.OLS.swap[1,]), sd(b.OLS.swap[2,]), sd(b.OLS.swap[3,]))

##### stationary DNS model 
var.level <- Arima(b.OLS.swap[1,], c(1,0,0), method="ML")
p.swap.stat[7] <- p.double.stat[13] <- coefficients(var.level)[2] # c
p.swap.stat[10] <- coefficients(var.level)[1] # Phi
p.swap.stat[14] <- p.double.stat[33] <- sd(fitted(var.level), na.rm=TRUE) # Q

var.slope <- Arima(b.OLS.swap[2,], c(1,0,0), method="ML")
p.swap.stat[8] <- p.double.stat[14] <- coefficients(var.slope)[2] # c
p.swap.stat[11] <- coefficients(var.slope)[1] # Phi
p.swap.stat[15] <- p.double.stat[34] <- sd(fitted(var.slope), na.rm=TRUE) # Q

var.curvature <- Arima(b.OLS.swap[3,], c(1,0,0), method="ML")
p.swap.stat[9] <- p.double.stat[15] <- coefficients(var.curvature)[2] # c
p.swap.stat[12] <- coefficients(var.curvature)[1] # Phi
p.swap.stat[16] <- p.double.stat[35] <- sd(fitted(var.curvature), na.rm=TRUE) # Q

##### BOND ====================================================================
# Estimate the H matrix entries (same for stationary and nonstationary models)
p.bond.stat[17:(17 + length(m.bond) - 1)] <- apply(residuals.OLS.bond, 1, sd, na.rm=TRUE)
p.double.stat[40] <- mean(apply(residuals.OLS.bond, 1, sd))

# Estimate b_0 and P_0
p.bond.stat[1:3] <- b.OLS.bond[,1] # b_0 is equal to first estimated factor
p.bond.stat[4:6] <- c(sd(b.OLS.bond[1,], na.rm = TRUE), sd(b.OLS.bond[2,], na.rm = TRUE), sd(b.OLS.bond[3,], na.rm = TRUE))
p.double.stat[4:6] <- b.OLS.bond[,1]
p.double.stat[10:12] <- c(sd(b.OLS.swap[1,]), sd(b.OLS.swap[2,]), sd(b.OLS.swap[3,]))

##### stationary DNS model 
var.level <- Arima(b.OLS.bond[1,], c(1,0,0), method="ML")
p.bond.stat[7] <- p.double.stat[16] <- coefficients(var.level)[2] # c
p.bond.stat[10] <- coefficients(var.level)[1] # Phi
p.bond.stat[14] <- p.double.stat[36] <- sd(fitted(var.level), na.rm=TRUE) # Q

var.slope <- Arima(b.OLS.bond[2,], c(1,0,0), method="ML")
p.bond.stat[8] <- p.double.stat[17] <- coefficients(var.slope)[2] # c
p.bond.stat[11] <- coefficients(var.slope)[1] # Phi
p.bond.stat[15] <- p.double.stat[37] <- sd(fitted(var.slope), na.rm=TRUE) # Q

var.curvature <- Arima(b.OLS.bond[3,], c(1,0,0), method="ML")
p.bond.stat[9] <- p.double.stat[18] <- coefficients(var.curvature)[2] # c
p.bond.stat[12] <- coefficients(var.curvature)[1] # Phi
p.bond.stat[16] <- p.double.stat[38] <- sd(fitted(var.curvature), na.rm=TRUE) # Q

##### SIMULTANEOUS ============================================================
beta.double <- rbind(b.OLS.swap, b.OLS.bond)
# estimate a restricted VAR model
suppressWarnings({
  varmodel <- VAR(t(beta.double), type="const", p=1)
  restrict <- cbind(matrix(c(1,0,0,1,0,0,
                             0,1,0,0,1,0,
                             0,0,1,0,0,1), 6, 6, byrow=TRUE),
                    rep(1,6)) # matrix containing the trend coefficients
  var.restr <- restrict(varmodel, method = "man", resmat = restrict)
})
beta.OLS.var <- fitted(var.restr) %>% t()

p.double.stat[19] <- coef(var.restr)$y1[1,1] # Phi_11
p.double.stat[20] <- coef(var.restr)$y1[2,1] # Phi_14
p.double.stat[21] <- coef(var.restr)$y2[1,1] # Phi_22
p.double.stat[22] <- coef(var.restr)$y2[2,1] # Phi_25
p.double.stat[23] <- coef(var.restr)$y3[2,1] # Phi_33
p.double.stat[24] <- coef(var.restr)$y3[1,1] # Phi_36

p.double.stat[25] <- coef(var.restr)$y4[1,1] # Phi_41
p.double.stat[26] <- coef(var.restr)$y4[2,1] # Phi_44
p.double.stat[27] <- coef(var.restr)$y5[1,1] # Phi_52
p.double.stat[28] <- coef(var.restr)$y5[2,1] # Phi_55
p.double.stat[29] <- coef(var.restr)$y6[1,1] # Phi_63
p.double.stat[30] <- coef(var.restr)$y6[2,1] # Phi_66

##### Export all results ======================================================
files <- list.files(wd.result)

# SWAP STATIONARY
if("p_swap_stat.csv" %in% files){ # if file already exists -> do not overwrite!
  csv.swap.stat <- read.csv(str_c(wd.result, "p_swap_stat.csv")) %>% dplyr::select(-X)
  # check if initial parameters are changed:
  if(!all(csv.swap.stat$OLS - p.swap.stat < 1e-14)){ # afronding dus niet exact gelijk
    cat("WARNING: initial parameters have changed for swap stationary.")
  }
} else {
  csv.swap.stat <- data.frame(settings = rep(0, length(p.swap.stat)),
                              OLS = p.swap.stat)
  csv.swap.stat$settings[1] <- 1 # convergence indicator
  write.csv(csv.swap.stat, str_c(wd.result, "p_swap_stat.csv"))
}

# BOND STATIONARY
if("p_bond_stat.csv" %in% files){ # if file already exists -> do not overwrite!
  csv.bond.stat <- read.csv(str_c(wd.result, "p_bond_stat.csv")) %>% dplyr::select(-X)
  # check if initial parameters are changed:
  if(!all(csv.bond.stat$OLS - p.bond.stat < 1e-14, na.rm = TRUE)){ # afronding dus niet exact gelijk
    cat("WARNING: initial parameters have changed for bond stationary.")
  }
} else {
  csv.bond.stat <- data.frame(settings = rep(0, length(p.bond.stat)),
                              OLS = p.bond.stat)
  csv.bond.stat$settings[1] <- 1 # convergence indicator
  write.csv(csv.bond.stat, str_c(wd.result, "p_bond_stat.csv"))
}

# DOUBLE STATIONARY
if("p_double_stat.csv" %in% files){ # if file already exists -> do not overwrite!
  csv.double.stat <- read.csv(str_c(wd.result, "p_double_stat.csv")) %>% dplyr::select(-X)
  # check if initial parameters are changed:
  if(!all(csv.double.stat$OLS - p.bond.nonstat < 1e-14, na.rm = TRUE)){ # afronding dus niet exact gelijk
    cat("WARNING: initial parameters have changed for bond nonstationary.")
  }
} else {
  csv.double.stat <- data.frame(settings = rep(0, length(p.double.stat)),
                                OLS = p.double.stat)
  csv.double.stat$settings[1] <- 1 # convergence indicator
  write.csv(csv.double.stat, str_c(wd.result, "p_double_stat.csv"))
}


