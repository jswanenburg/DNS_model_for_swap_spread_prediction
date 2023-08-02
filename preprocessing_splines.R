# script to convert the bond data into fixed set of maturities
# bond data had too many different tau values --> MLE estimation did not converge in over 3 days
# this data preprocessing step reduces the number of different maturities from 308 to 30

# Import all data
wd.data <- <ENTER YOUR OWN DATA WORKING DIRECTORY HERE>
file <- <ENTER YOUR OWN FILENAME HERE>
data <- read.csv(str_c(wd.data, file)) %>% dplyr::select(-X)

# select bond data only
data.bond <- rbind(split(data, data$type)$ZCB, split(data, data$type)$CB)
data.swap <- split(data, data$type)$swap # save for later to merge with bond data after splining
data.bond.per.date <- split(data.bond, data.bond$date)

taus <- c(0.4, 0.525, 0.65, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2,
          2.4, 2.8, 3.2, 3.6, 4, 4.4, 4.8,
          5.4, 6.4, 7.4, 8.4, 9.4, 10.4, 11.5, 13.5,
          16, 19, 22, 25, 28)

data.bond.per.date[[200]] %>%
  ggplot(aes(x=tau, y=rate)) + geom_point() +
  geom_vline(xintercept = taus, alpha=.5); cat(length(taus))


# select date
for(i in names(data.bond.per.date)){
  # cat(i)
  
  df <- data.bond.per.date[[i]]
  
  add.df <- c()
  for(j in 1:(length(taus))){
    # als tau in de bond data zit, gebruik dan de true value
    if(taus[j] %in% df$tau){
      ir <- df$rate[which(df$tau == taus[j])]
    }
    
    # select data for spline: 4 closest points
    data.for.spline <- (df %>% 
                          mutate(dist = abs(tau - taus[j])) %>%
                          arrange(dist))[1:4,]
    
    # check for double values
    extra <- 1
    while(length(unique(data.for.spline$tau))<4){
      data.for.spline <- (df %>% 
                            mutate(dist = abs(tau - taus[j])) %>%
                            arrange(dist))[1:(4+extra),]
    }; rm(extra)
    
    # fit the spline
    ss <- smooth.spline(x = data.for.spline$tau, 
                        y = data.for.spline$rate)
    ir <- predict(ss, taus[j])$y
    
    add.row <- data.frame(
      date = i, ticker = NA_character_,
      rate = ir, tau = taus[j],
      type = "spline"
    ); add.df <- rbind(add.df, add.row)
    
  }
  data.bond.per.date[[i]] <- rbind(
    data.bond.per.date[[i]],
    add.df
  )
}

new_master_dataframe <- rbind(
  as.data.frame(do.call(rbind, data.bond.per.date)),
  data.swap
)

write.csv(new_master_dataframe,
          str_c(wd.data, "2023-06-15_master_with_splines.csv"),
          row.names = FALSE)
