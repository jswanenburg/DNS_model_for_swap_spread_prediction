# Initial parameters if the factor dynamics are nonstationary
##### Create all empty parameter vectors 
p.swap.nonstat <- rep(NA_real_, 17 + length(m.swap) - 1)
p.swap.nonstat[13] <- log(lambda.swap)
p.bond.nonstat <- rep(NA_real_, 17 + length(m.bond) - 1)
p.bond.nonstat[13] <- log(lambda.bond)
p.double.nonstat <- rep(NA_real_, 40)
p.double.nonstat[31:32] <- c(log(lambda.swap), log(lambda.bond))


##### SWAP ====================================================================
# Estimate the H matrix entries (same for stationary and nonstationary models)
p.swap.nonstat[17:(17 + length(m.swap) - 1)] <- apply(residuals.OLS.swap, 1, sd)
p.double.nonstat[39] <- mean(apply(residuals.OLS.swap, 1, sd))

# Estimate b_0 
p.swap.nonstat[1:3]  <- b.OLS.swap[,1] # b_0 is equal to first estimated factor
p.double.nonstat[1:3] <- b.OLS.swap[,1]

# Estimate P_0
p.swap.nonstat[4:6]  <- c(sd(b.OLS.swap[1,]), sd(b.OLS.swap[2,]), sd(b.OLS.swap[3,]))
p.double.nonstat[7:9] <- c(sd(b.OLS.swap[1,]), sd(b.OLS.swap[2,]), sd(b.OLS.swap[3,]))

##### nonstationary DNS model 
var.level <- Arima(b.OLS.swap.diff[1,], c(1,0,0), method="ML")
p.swap.nonstat[7] <- coefficients(var.level)[2] # c
p.swap.nonstat[10] <- coefficients(var.level)[1] # Phi
p.swap.nonstat[14] <- sd(fitted(var.level), na.rm=TRUE) # Q

var.slope <- Arima(b.OLS.swap.diff[2,], c(1,0,0), method="ML")
p.swap.nonstat[8] <- coefficients(var.slope)[2] # c
p.swap.nonstat[11] <- coefficients(var.slope)[1] # Phi
p.swap.nonstat[15] <- sd(fitted(var.slope), na.rm=TRUE) # Q

var.curvature <- Arima(b.OLS.swap.diff[3,], c(1,0,0), method="ML")
p.swap.nonstat[9] <- coefficients(var.curvature)[2] # c
p.swap.nonstat[12] <- coefficients(var.curvature)[1] # Phi
p.swap.nonstat[16] <- sd(fitted(var.curvature), na.rm=TRUE) # Q

##### BOND ====================================================================
# Estimate the H matrix entries (same for stationary and nonstationary models)
p.bond.nonstat[17:(17 + length(m.bond) - 1)] <- apply(residuals.OLS.bond, 1, sd, na.rm=TRUE)
p.double.nonstat[40] <- mean(apply(residuals.OLS.swap, 1, sd))

# Estimate b_0 
p.bond.nonstat[1:3]  <- b.OLS.bond[,1] # b_0 is equal to first estimated factor
p.double.nonstat[4:6] <- b.OLS.bond[,1]

# Estimate P_0
p.bond.nonstat[4:6] <- c(sd(b.OLS.bond[1,], na.rm = TRUE), 
                         sd(b.OLS.bond[2,], na.rm = TRUE),
                         sd(b.OLS.bond[3,], na.rm = TRUE))
p.double.nonstat[10:12] <- c(sd(b.OLS.swap[1,]), sd(b.OLS.swap[2,]), sd(b.OLS.swap[3,]))

##### nonstationary DNS model 
var.level <- Arima(b.OLS.bond.diff[1,], c(1,0,0), method="ML")
p.bond.nonstat[7] <- coefficients(var.level)[2] # c
p.bond.nonstat[10] <- coefficients(var.level)[1] # Phi
p.bond.nonstat[14] <- sd(fitted(var.level), na.rm=TRUE) # Q

var.slope <- Arima(b.OLS.bond.diff[2,], c(1,0,0), method="ML")
p.bond.nonstat[8] <- coefficients(var.slope)[2] # c
p.bond.nonstat[11] <- coefficients(var.slope)[1] # Phi
p.bond.nonstat[15] <- sd(fitted(var.slope), na.rm=TRUE) # Q

var.curvature <- Arima(b.OLS.bond.diff[3,], c(1,0,0), method="ML")
p.bond.nonstat[9] <- coefficients(var.curvature)[2] # c
p.bond.nonstat[12] <- coefficients(var.curvature)[1] # Phi
p.bond.nonstat[16] <- sd(fitted(var.curvature), na.rm=TRUE) # Q

##### Stationary VAR model double rates
beta.double <- rbind(b.OLS.swap, b.OLS.bond)
beta.double.diff <- rbind(b.OLS.swap.diff, b.OLS.bond.diff)
# estimate a restricted VAR model
suppressWarnings({
  b.for.var <- t(beta.double.diff)[-1,]
  colnames(b.for.var) <- c("y1", "y2", "y3", "y4", "y5", "y6")
  varmodel.diff <- VAR(b.for.var, type="const", p=1)
  restrict.diff <- cbind(matrix(c(1,0,0,1,0,0,
                                  0,1,0,0,1,0,
                                  0,0,1,0,0,1), 6, 6, byrow=TRUE),
                         rep(1,6)) # matrix containing the trend coefficients
  var.restr.diff <- restrict(varmodel.diff, method = "man", resmat = restrict.diff)
})
################### TOT AAN HIER
beta.OLS.var <- fitted(var.restr.diff) %>% t()

# c
p.double.nonstat[13] <- coef(var.restr.diff)$y1[3,1]
p.double.nonstat[14] <- coef(var.restr.diff)$y2[3,1]
p.double.nonstat[15] <- coef(var.restr.diff)$y3[3,1]
p.double.nonstat[16] <- coef(var.restr.diff)$y4[3,1]
p.double.nonstat[17] <- coef(var.restr.diff)$y5[3,1]
p.double.nonstat[18] <- coef(var.restr.diff)$y6[3,1]

# Phi
p.double.nonstat[19] <- coef(var.restr.diff)$y1[1,1] # Phi_11
p.double.nonstat[20] <- coef(var.restr.diff)$y1[2,1] # Phi_14
p.double.nonstat[21] <- coef(var.restr.diff)$y2[1,1] # Phi_22
p.double.nonstat[22] <- coef(var.restr.diff)$y2[2,1] # Phi_25
p.double.nonstat[23] <- coef(var.restr.diff)$y3[2,1] # Phi_33
p.double.nonstat[24] <- coef(var.restr.diff)$y3[1,1] # Phi_36

p.double.nonstat[25] <- coef(var.restr.diff)$y4[1,1] # Phi_41
p.double.nonstat[26] <- coef(var.restr.diff)$y4[2,1] # Phi_44
p.double.nonstat[27] <- coef(var.restr.diff)$y5[1,1] # Phi_52
p.double.nonstat[28] <- coef(var.restr.diff)$y5[2,1] # Phi_55
p.double.nonstat[29] <- coef(var.restr.diff)$y6[1,1] # Phi_63
p.double.nonstat[30] <- coef(var.restr.diff)$y6[2,1] # Phi_66
 
# Q
p.double.nonstat[33] <- sd(fitted(var.restr.diff)[1,])
p.double.nonstat[34] <- sd(fitted(var.restr.diff)[2,])
p.double.nonstat[35] <- sd(fitted(var.restr.diff)[3,])
p.double.nonstat[36] <- sd(fitted(var.restr.diff)[4,])
p.double.nonstat[37] <- sd(fitted(var.restr.diff)[5,])
p.double.nonstat[38] <- sd(fitted(var.restr.diff)[6,])

##### Export all initial parameters to CSV 
files <- list.files(wd.result)

# SWAP NONSTATIONARY
if("p_swap_nonstat.csv" %in% files){ # if file already exists -> do not overwrite!
  csv.swap.nonstat <- read.csv(str_c(wd.result, "p_swap_nonstat.csv")) %>% dplyr::select(-X)
  # check if initial parameters are changed:
  if(!all(csv.swap.nonstat$OLS - p.swap.nonstat < 1e-14)){ # afronding dus niet exact gelijk
    cat("WARNING: initial parameters have changed for swap nonstationary.")
  }
} else {
  csv.swap.nonstat <- data.frame(settings = rep(0, length(p.swap.nonstat)),
                                 OLS = p.swap.nonstat)
  csv.swap.nonstat$settings[1] <- 1 # convergence indicator
  write.csv(csv.swap.nonstat, str_c(wd.result, "p_swap_nonstat.csv"))
}

# BOND NONSTATIONARY
if("p_bond_nonstat.csv" %in% files){ # if file already exists -> do not overwrite!
  csv.bond.nonstat <- read.csv(str_c(wd.result, "p_bond_nonstat.csv")) %>% dplyr::select(-X)
  # check if initial parameters are changed:
  if(!all(csv.bond.nonstat$OLS - p.bond.nonstat < 1e-14, na.rm = TRUE)){ # afronding dus niet exact gelijk
    cat("WARNING: initial parameters have changed for bond nonstationary.")
  }
} else {
  csv.bond.nonstat <- data.frame(settings = rep(0, length(p.bond.nonstat)),
                                 OLS = p.bond.nonstat)
  csv.bond.nonstat$settings[1] <- 1 # convergence indicator
  write.csv(csv.bond.nonstat, str_c(wd.result, "p_bond_nonstat.csv"))
}

# DOUBLE NONSTATIONARY
if("p_double_nonstat.csv" %in% files){ # if file already exists -> do not overwrite!
  csv.double.nonstat <- read.csv(str_c(wd.result, "p_double_nonstat.csv")) %>% dplyr::select(-X)
  # check if initial parameters are changed:
  if(!all(csv.double.nonstat$OLS - p.double.nonstat < 1e-14, na.rm = TRUE)){ # afronding dus niet exact gelijk
    cat("WARNING: initial parameters have changed for double nonstationary.")
  }
} else {
  csv.double.nonstat <- data.frame(settings = rep(0, length(p.double.nonstat)),
                                OLS = p.double.nonstat)
  csv.double.nonstat$settings[1] <- 1 # convergence indicator
  write.csv(csv.double.nonstat, str_c(wd.result, "p_double_nonstat.csv"))
}
