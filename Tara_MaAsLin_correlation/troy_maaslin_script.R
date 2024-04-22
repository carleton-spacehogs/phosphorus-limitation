# Troy's code for using MaAsLin, building on Jimmy's work in
# nutrient-COG-correlation.Rmd and .init_tara.R
# Can be run from the command line with >R troy_maaslin_script.R


# Code from Jimmy (see nutrient-COG-correlation.Rmd)

#knitr::opts_chunk$set(echo = TRUE)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
source(".init_tara.R") # important for importing Jimmy's functions (thanks Jimmy)

tara_md = init_tara_metadata()
COGabun = get_COGabundance() # only the prokaryote fraction

tara_md = tara_md[colnames(COGabun), ]

result_dir = 'MaAsLin_out/logPO4_sqrtIron_logNO2NO3_Salinity_Oxygen_temp_logDepth_size'
param = c('size_fraction','Mean_Temperature', "log_depth", 'Mean_Salinity',
           'Mean_Oxygen', 'log_PO4','log_NO2NO3', "sqrt_iron")
# end of Jimmy's code

# Jimmy's MaAsLin call:
# 
# fit_data <- Maaslin2(
#     COGabun, tara_md, result_dir, transform = "LOG",
#     fixed_effects = param,
#     random_effects = c("Mean_Lat", "Mean_Long"),
#     reference = 'size_fraction,0.22-1.6',
#     max_significance = 0.05,
#     normalization = 'NONE',
#     min_prevalence = 0.33, # at least 1/3 of the sample has it
#     standardize = FALSE)


# taking the log of the abundance data before running my own MaAsLin

my_log <- function(x) {
  return(log10(x + 1e-9)) # add a small number to avoid log(0)
}                         # choose 1e-9 because a few values are below 1e-8
COGabun_log <- sapply(COGabun[1:ncol(COGabun)], my_log)

fit_data <- Maaslin2(
    COGabun_log, tara_md, result_dir, transform = "LOG",
    fixed_effects = param,
    random_effects = c("Mean_Lat", "Mean_Long"),
    reference = 'size_fraction,0.22-1.6',
    max_significance = 0.05,
    normalization = 'NONE',
    min_prevalence = 0.33,
    standardize = FALSE,
    cores = 40, # change as necessary
    save_models = TRUE
    )