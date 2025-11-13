# Author: Troy Osborn
# Date: March 2025
# A script for an analysis I'm doing for Joanne looking at just relationships
# with nitrogen. First I'll run MaAsLin for an exponential relationship between
# COG abundance and log_NO2NO3 (so, a linear relationship between logAbundance 
# and log_NO2NO3), then I'll
# look at NO2 and NO3 individually and see if any COGs have opposite-parity 
# relationships with the two. Some code taken from maaslin_modelling.R

### load everything in
source(".init_tara.R") 

tara_md = init_tara_metadata()
COGabun = get_COGabundance() # only the prokaryote fraction

tara_md = tara_md[colnames(COGabun), ]

### sanity check

# run the desired lm on a random COG to see if maaslin results match it, for sanity
df.COG5598 <- 
  as.data.frame(t(COGabun["COG5598",])) |> 
  rownames_to_column(var = "measurement") |> 
  left_join(tara_md["log_NO2NO3"] |> 
              as.data.frame() |> 
              rownames_to_column(var = "measurement"),
            by = "measurement")

lm.COG5598 <- lm(data = df.COG5598, log(COG5598) ~ log_NO2NO3)
summary(lm.COG5598)


### Maaslin run
# try just changing transform = "NONE" to transform = "LOG"

output_dir <- "jo_maaslin_out/exponential_with_transformation"
fit_data <- Maaslin2(
  COGabun, tara_md, output_dir, 
  transform = "LOG",
  fixed_effects = "log_NO2NO3",
  max_significance = 0.05, # alpha
  normalization = 'NONE',
  min_prevalence = 0.33,
  standardize = FALSE,
  min_abundance = 0.0,
  cores = 40, # change as necessary
  save_models = TRUE, # change this if output folder is too big
  analysis_method = "LM" # must be LM because CPLM will log already-logged response
)

### import the models RDS file and sanity check
#exp_models[["COG5598"]] |> summary()
#summary(lm.COG5598)
# slightly different coefficients, but *same exact t-values*, so the models are 
# essentially the same. I predict the difference is due to Maaslin adding a small
# number to the response before logging it, which I did not do.

# Maaslin's graphs are on the original scale for the response, so I'll make my own
# see request 3 of collate_graphs.Rmd