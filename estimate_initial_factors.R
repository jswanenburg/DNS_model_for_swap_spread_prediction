##### Estimate initial factors ================================================
##### SWAP 
# cat("\nStart estimation of swap factors with lambda =", lambda.swap)

# Create a matrix of factor loadings for lambda
pars.temp.swap <- rep(NA_real_, 17 + length(m.swap) - 1)
pars.temp.swap[13] <- log(lambda.swap)
Lambda.df <- par_ssm_stat_swap(pars.temp.swap, m.swap)$Hm %>% as.data.frame()
colnames(Lambda.df) <- c("l1", "l2", "l3")

b.OLS.swap <- c(); residuals.OLS.swap <- c() # empty vectors to store results
for(t in 1:ncol(y.swap.train)){
  Lambda.df$rate <- y.swap.train[,t]
  lm.fit <- lm(rate ~ l2 + l3, data=Lambda.df)
  b.OLS.swap <- cbind(b.OLS.swap, matrix(lm.fit$coefficients, 3,1))
  residuals.OLS.swap <- cbind(residuals.OLS.swap, residuals(lm.fit))
}

##### SWAP DIFF 
b.OLS.swap.diff <- cbind(rep(NA_real_, 3), # first element of the differences is NA, because beta_0 is not existent
                         t(diff(t(b.OLS.swap))))
rownames(b.OLS.swap.diff) <- c("y1", "y2", "y3")

##### BOND 
# cat("\nStart estimation of bond factors with lambda =", lambda.bond)

# Create a matrix of factor loadings for lambda
pars.temp.bond <- rep(NA_real_, 17 + length(m.bond) - 1)
pars.temp.bond[13] <- log(lambda.bond)
Lambda.df <- par_ssm_stat_bond(pars.temp.bond, m.bond)$Hm %>% as.data.frame()
colnames(Lambda.df) <- c("l1", "l2", "l3")

b.OLS.bond <- c(); residuals.OLS.bond <- c() # empty vectors to store results
for(t in 1:ncol(y.bond.train)){
  Lambda.df$rate <- y.bond.train[,t]
  lm.fit <- lm(rate ~ l2 + l3, data=Lambda.df)
  b.OLS.bond <- cbind(b.OLS.bond, matrix(lm.fit$coefficients, 3,1))
  res.vector <- rep(NA_real_, nrow(Lambda.df))
  res.vector[which(!is.na(Lambda.df$rate))] <- residuals(lm.fit)
  residuals.OLS.bond <- cbind(residuals.OLS.bond, res.vector) # HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

##### BOND DIFF 
b.OLS.bond.diff <- cbind(rep(NA_real_, 3), # first element of the differences is NA, because beta_0 is not existent
                         t(diff(t(b.OLS.bond))))
rownames(b.OLS.bond.diff) <- c("y1", "y2", "y3")
