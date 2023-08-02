##### CODE FOR EMPIRICAL DATA
rm(list = ls())
##### Load required packages ==================================================
packages <- c("MASS", "patchwork", "dplyr", "ggplot2", "tidyr", "stringr", 
              "RColorBrewer", "ggrepel", "vars", "kalmanfilter", "forecast",
              "tseries")
lapply(packages, require, character.only = TRUE)

wd.main <- <ENTER YOUR OWN WORKING DIRECTORY HERE>
wd.code <- str_c(wd.main, "Code/")
wd.data <- str_c(wd.main, "Data/")

##### Load functions
source(str_c(wd.code, "functions.R"))

##### Choose settings for estimation -> created with 'create_settings.R'
for(settings in c(2)){
  cat("\n\n=========================================================",
      "\nSettings: ", settings, "============================================",
      "\n=========================================================")
  wd.result <- str_c(wd.main, "Results/", settings,"/")
  source(str_c(wd.code, "load_settings.R"))
  
  ##### Load data
  file <- <ENTER YOUR OWN FILENAME HERE>
  source(str_c(wd.code, "load_data.R")) # file that loads data and converts it to the correct shape
  source(str_c(wd.code, "empirical_factors.R"))
  
  ##### Estimate initial parameters
  source(str_c(wd.code, "estimate_initial_factors.R"))
  if(settings %in% c(7)){ # manually choose if you want to estimate the statDNS or the nonstatDNS models
    source(str_c(wd.code, "estimate_initial_parameters_nonstat.R"))
  } else if(settings == 4){
    source(str_c(wd.code, "estimate_initial_parameters_stat.R"))
  }
  
  ##### MLE
  source(str_c(wd.code, "MLE.R"))
} # end loop for settings

