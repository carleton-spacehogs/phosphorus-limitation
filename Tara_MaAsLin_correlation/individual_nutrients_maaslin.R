# Author: Troy Osborn
# Date: Oct/Nov 2024
# Purpose: to run MaAsLin2 multiple times on individual nutrients
# (see troy_maaslin_script.R)


source(".init_tara.R")

tara_md <- init_tara_metadata()
COGabun <- get_COGabundance() 

tara_md <- tara_md[colnames(COGabun), ]
tara_md$Absolute_Latitude <- abs(tara_md$Mean_Lat)

nonnutrients  <-  c('size_fraction','Mean_Temperature', "log_depth", 
                    'Mean_Salinity', 'Absolute_Latitude')
nutrients <- c('Mean_Oxygen', 'log_PO4','log_NO2NO3', "sqrt_iron")

# remove the apparently problematic COG0787
COGabun <- COGabun |>
  mutate(id = rownames(COGabun)) |> 
  filter(id != "COG0787") |> 
  select(-id)


for (i in 1:length(nutrients)) {

  current_nutrient <- nutrients[i]

  predictors <- c(nonnutrients, current_nutrient)
  result_dir <- paste0("MaAsLin_out/individual_nutrients/", current_nutrient)

  fit_data <- Maaslin2(
    COGabun, tara_md, result_dir,
    fixed_effects = predictors,
    transform = "NONE",
    reference = 'size_fraction,0.22-1.6',
    max_significance = 0.05,
    normalization = 'NONE',
    min_prevalence = 0.33,
    standardize = FALSE,
    min_abundance = 0.0,
    cores = 1, # change as necessary
    save_models = F,
    analysis_method = "CPLM"
  )
}

# library("devtools")
# install_github("biobakery/maaslin3")

### TEST JUST A FEW OBSERVATIONS ###





