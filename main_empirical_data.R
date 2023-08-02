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
  file <- "2023-06-15_master_with_splines.csv"
  source(str_c(wd.code, "load_data.R"))
  # source(str_c(wd.code, "empirical_factors.R"))
  
  ##### Estimate initial parameters
  source(str_c(wd.code, "initial_factors.R"))
  if(settings %in% c(7)){
    source(str_c(wd.code, "initial_parameters_nonstat.R"))
  } else if(settings == 4){
    source(str_c(wd.code, "initial_parameters_stat.R"))
  }
  
  ##### MLE
  source(str_c(wd.code, "MLE.R"))
} # end loop for settings

