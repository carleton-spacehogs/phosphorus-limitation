# Troy's code for using MaAsLin, building on Jimmy's work in
# nutrient-COG-correlation.Rmd and .init_tara.R
# Can be run from the command line with >R troy_maaslin_script.R


# Code from Jimmy (see nutrient-COG-correlation.Rmd)

#knitr::opts_chunk$set(echo = TRUE)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path ))
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

# changing result_dir to my own
result_dir <- "MaAsLin_out/logCOGabun_against_logPO4_sqrtIron_logNO2NO3_Salinity_Oxygen_temp_logDepth_size"

# taking the log of the abundance data before running my own MaAsLin
my_log <- function(x) {
  return(log10(x + 1e-9)) # add a small number to avoid log(0)
}                         # choose 1e-9 because a few values are below 1e-8
COGabun_log <- as.data.frame(sapply(COGabun[1:ncol(COGabun)], my_log))
rownames(COGabun_log) <- rownames(COGabun)

# MaAsLin call
# fit_data <- Maaslin2(
#     COGabun_log, tara_md, result_dir,
#     fixed_effects = param,
#     random_effects = c("Mean_Lat", "Mean_Long"),
#     reference = 'size_fraction,0.22-1.6',
#     max_significance = 0.05,
#     normalization = 'NONE',
#     min_prevalence = 0.33,
#     standardize = FALSE,
#     min_abundance = log10(1e-9), # log(1e-9) = -9, doing this to filter out cases where abun = 0
#     cores = 40, # change as necessary
#     save_models = TRUE,
#     transform = 'NONE'
#     )
# gonna keep this around for now, but its output is largely the same as Jimmy's...
# see MaAsLin_out/logCOGabun_against_logPO4_sqrtIron_logNO2NO3_Salinity_Oxygen_temp_logDepth_size
# he log-transformed within MaAsLin, I did it before calling MaAsLin

tara_md$Absolute_Latitude <- abs(tara_md$Mean_Lat)
param_plusAbsLat <- c(param, "Absolute_Latitude")

# WHAT I'VE TRIED | WHAT HAPPENED WHEN I TRIED IT
# logging the abundance (with my own fxn) | same as Jimmy's (roughly)
# poisson model (keeping the log transform, in MaAsLin call) | no associations found
# poisson model (no log transform) | lots of strong correlations, likely too strong actually
# negative binomial model (no log transform) | no associations found 
# (moving forward with cplm (poisson) analysis method)
# cplm without random effects | associations look fine, notably smaller q-values (NOT a problem necessarily, earlier model may have been overfit)
# cplm, no ranefs, with absolute latitude | looks mostly the same as before, but with a strange heatmap (still need to investigate)

# now try removing random effects
result_dir <- "MaAsLin_out/cplm-logCOGabun_against_no.ranefs_AbsLat_logPO4_sqrtIron_logNO2NO3_Salinity_Oxygen_temp_logDepth_size"
fit_data <- Maaslin2(
    COGabun, tara_md, result_dir, transform = "NONE",
    fixed_effects = param_plusAbsLat,
#    random_effects = c("Mean_Lat", "Mean_Long"),
    reference = 'size_fraction,0.22-1.6',
    max_significance = 0.05,
    normalization = 'NONE',
    min_prevalence = 0.33,
    standardize = FALSE,
    min_abundance = 0.0,
    cores = 40, # change as necessary
    save_models = TRUE,
    analysis_method = "CPLM"
    )