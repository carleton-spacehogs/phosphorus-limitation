# Author: Troy Osborn
# Troy's code for using MaAsLin, building on Jimmy's work in
# nutrient-COG-correlation.Rmd and .init_tara.R
# Previously named: troy_maaslin_script.R

# Code from Jimmy (see nutrient-COG-correlation.Rmd)

source(".init_tara.R") 

tara_md = init_tara_metadata()
COGabun = get_COGabundance() # only the prokaryote fraction

tara_md = tara_md[colnames(COGabun), ]

param = c('size_fraction','Mean_Temperature', "log_depth", 'Mean_Salinity',
           'Mean_Oxygen', 'log_PO4','log_NO2NO3', "sqrt_iron")
# end of Jimmy's code


# create Absolute Latitude column
tara_md$Absolute_Latitude <- abs(tara_md$Mean_Lat)
param <- c(param, "Absolute_Latitude")


result_dir <- "MaAsLin_out/cplm-logCOGabun_against_no.ranefs_AbsLat_logPO4_sqrtIron_logNO2NO3_Salinity_Oxygen_temp_logDepth_size"
fit_data <- Maaslin2(
    COGabun, tara_md, result_dir, transform = "NONE",
    fixed_effects = param,
    reference = 'size_fraction,0.22-1.6',
    max_significance = 0.05, # alpha
    normalization = 'NONE',
    min_prevalence = 0.33,
    standardize = FALSE,
    min_abundance = 0.0,
    cores = 40, # change as necessary
    save_models = TRUE, # change this if output folder is too big
    analysis_method = "CPLM"
    )


# notes on previous iterations of the model
# WHAT I'VE TRIED | WHAT HAPPENED WHEN I TRIED IT
# logging the abundance (with my own fxn) | same as Jimmy's (roughly)
# poisson model (keeping the log transform, in MaAsLin call) | no associations found
# poisson model (no log transform) | lots of strong correlations, likely too strong actually
# negative binomial model (no log transform) | no associations found 
# (moving forward with cplm (poisson) analysis method)
# cplm without random effects | associations look fine, notably smaller q-values (NOT a problem necessarily, earlier model may have been overfit)
# cplm, no ranefs, with absolute latitude | looks mostly the same as before, but with a strange heatmap (still need to investigate)
