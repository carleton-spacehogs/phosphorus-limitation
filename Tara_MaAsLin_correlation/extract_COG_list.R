# Author: Troy Osborn
# Date: Sept 2024

library(tidyverse)
library(tidymodels)

# Load the data
results <- read_tsv("MaAsLin_out/cplm-logCOGabun_against_no.ranefs_AbsLat_logPO4_sqrtIron_logNO2NO3_Salinity_Oxygen_temp_logDepth_size/significant_results.tsv")

COGs <- results$feature

# concatenate "metadata" and "value" if they are not the same
# this is just so size_fraction indicates the size fraction
# metadata <- rep(NA, nrow(results))
# for (i in 1:nrow(results)) {
#     if (results$metadata[i] == results$value[i]) {
#         metadata[i] = results$metadata[i]
#     } else {
#         metadata[i] = paste(results$metadata[i], results$value[i], sep = "")
#     }
# }

metadata <- results$metadata

# notate positive/negative correlation
posnegs <- ifelse(results$coef > 0, "positive", "negative")


COG_df <- data.frame(COG = COGs, metadata = metadata, polarity = posnegs)


## spearman correlations ##


# first, grab the data
source(".init_tara.R")
tara_md = init_tara_metadata()
COGabun = get_COGabundance()

tara_md = tara_md[colnames(COGabun), ] # exclude a bunch of rows from tara_md so that cor.test() works (needs vectors of same length) 

# make the corr function
get_corr <- function(COG, metadata_column, threshold = 0.6) {
  cog_abundances = as.numeric(COGabun[COG,]) # row of COGabun corresponding to given COG
  
  if (metadata_column != "Absolute_Latitude") {
    station_metadata_values = tara_md[,metadata_column] # column of tara_md corresponding to given metadata
  } else {
    station_metadata_values = abs(tara_md[,"Mean_Lat"]) # transform Mean_Lat appropriately
  }
  
  cor_obj <- cor.test(cog_abundances, station_metadata_values, 
                      method = "spearman",
                      exact = F)
  cor_obj$estimate
}

# loop through COG_df and store correlations
corrs <- rep(NA, nrow(COG_df))
for (i in 1:nrow(COG_df)) {
  if (COG_df$metadata[i] != "size_fraction") {
    corrs[i] <- get_corr(COG_df$COG[i], COG_df$metadata[i], 0.6)
  }
}

# add correlations and T/F column to COG_df
COG_df <- COG_df |> 
  mutate(correlation = corrs) |> 
  mutate(corr_above_0.6 = ifelse(abs(corrs) > 0.6, T, F))

# sanity checking
sum(is.na(COG_df$corr_above_0.6)) # should be only size_fraction rows
table(COG_df$metadata) # yep

sum(!is.na(COG_df$corr_above_0.6)) # valid values
table(COG_df$polarity, COG_df$corr_above_0.6) # number of "strong" correlations is just over 1000, much closer to "a three-digit number" than before


# export final COG_df
#write_tsv(COG_df, "significant_COG_list.tsv")




######### DIAGNOSTICS #########

# plot COG3850 (Signal transduction histidine kinase NarQ, nitrate/nitrite-specific) vs. nitrate/nitrite
COG3850 <- as.numeric(COGabun[rownames(COGabun) == "COG3850",])
ggplot(mapping = aes(y = COG3850, x = tara_md$log_NO2NO3)) +
  geom_point() +
  geom_smooth()
# jump

# plot COG5013 (Nitrate reductase alpha subunit) vs. nitrate/nitrite
COG5013 <- as.numeric(COGabun[rownames(COGabun) == "COG5013",])
ggplot(mapping = aes(y = COG5013, x = tara_md$log_NO2NO3)) +
  geom_point() +
  geom_smooth(se = F)
# big jump

# run one model by hand
model_5013 <- lm(COG5013 ~ tara_md$log_NO2NO3 + tara_md$size_fraction + tara_md$Mean_Temperature + 
             tara_md$log_depth + tara_md$Mean_Salinity + tara_md$Mean_Oxygen + 
             tara_md$log_PO4 + tara_md$log_NO2NO3 + tara_md$sqrt_iron)
tidy(model_5013)
# check variance inflation factors
car::vif(model_5013) # yep, looks like multicolinearity is at play here
# check correlations between predictors
X <- model.matrix(model_5013)
cor(X[,2:9])[1,] # yep, NO2NO3 is strongly correlated with PO4
# check the plot
ggplot(mapping = aes(y = COG5013, x = tara_md$log_PO4)) +
  geom_point() +
  geom_smooth(se = F) # nearly identical to the log_NO2NO3 graph...

# try quadratic terms
temp <- tibble(COG5013 = COG5013, NO2NO3 = tara_md$log_NO2NO3) |> drop_na()
model_sparse_5013 <- lm(data = temp, COG5013 ~ poly(NO2NO3, 2))
summary(model_sparse_5013) # significant quadratic term, this is expected

# try a quadratic term where we see a linear relationship
COG3958 <- as.numeric(COGabun[rownames(COGabun) == "COG3958",])
temp2 <- tibble(COG3958 = COG3958, NO2NO3 = tara_md$log_NO2NO3) |> drop_na()
model_sparse_3958 <- lm(data = temp2, COG3958 ~ poly(NO2NO3, 2))
summary(model_sparse_3958)
# quadratic term still significant, this is not expected...

# plot COG0166 (Glucose-6-phosphate isomerase) vs. phosphate
COG0166 <- as.numeric(COGabun[rownames(COGabun) == "COG0166",])
ggplot(mapping = aes(y = COG0166, x = tara_md$log_PO4)) +
  geom_point() +
  geom_smooth(se = F)
# bad

# plot COG0128 (5-enolpyruvylshikimate-3-phosphate synthase) vs. phosphate
COG0128 <- as.numeric(COGabun[rownames(COGabun) == "COG0128",])
ggplot(mapping = aes(y = COG0128, x = tara_md$log_PO4)) +
  geom_point() +
  geom_smooth(se = F)
# bad


# plot random COG vs random nutrient
COG1327 <- as.numeric(COGabun[rownames(COGabun) == "COG1327",])
ggplot(mapping = aes(y = COG1327, x = tara_md$log_PO4)) +
  geom_point() +
  geom_smooth(se = F)


COG2222 <- as.numeric(COGabun[rownames(COGabun) == "COG2222",])
ggplot(mapping = aes(y = COG2222, x = tara_md$log_PO4)) +
  geom_point() +
  geom_smooth(se = F)

COG0243 <- as.numeric(COGabun[rownames(COGabun) == "COG0243",])
ggplot(mapping = aes(y = COG0243, x = tara_md$log_NO2NO3)) +
  geom_point()

COG0243 <- as.numeric(COGabun[rownames(COGabun) == "COG0243",])
ggplot(mapping = aes(y = COG0243, x = tara_md$log_PO4)) +
  geom_point()

temp3 <- tibble(COG0243 = COG0243, NO2NO3 = tara_md$log_NO2NO3) |> drop_na()
model_sparse_COG0243 <- lm(data = temp3, COG0243 ~ poly(NO2NO3, 2))
summary(model_sparse_COG0243)
model_sparse_COG0243 <- lm(data = temp3, COG0243 ~ poly(NO2NO3, 2))
summary(model_sparse_COG0243)


######### END DIAGNOSTICS #########


# with individual nutrient data, remake COG list

nonnutrients  <-  c('size_fraction','Mean_Temperature', "log_depth", 
                    'Mean_Salinity', 'Absolute_Latitude')
nutrients <- c('Mean_Oxygen', 'log_PO4','log_NO2NO3', "sqrt_iron")

ind_nutrient_COGs <- tibble()

for (i in 1:(length(nutrients))) {
  current_nutrient <- nutrients[i]
  dir <- paste0("MaAsLin_out/individual_nutrients/", current_nutrient)
  nut_data <- read_tsv(paste0(dir, "/significant_results.tsv"))
  
  ind_nutrient_COGs <- ind_nutrient_COGs |> rbind(nut_data)
}

# make get_corr() play nice
get_corr <- Vectorize(get_corr)

ind_nutrient_COGs <- ind_nutrient_COGs |> 
  filter(metadata %in% nutrients) |> 
  select(-c(N, value, stderr, N.not.0, pval)) |> 
  rename("COG" = "feature") |> 
  mutate(polarity = ifelse(coef > 0, "positive", "negative")) |> 
  mutate(correlation = get_corr(COG, metadata)) |> 
  mutate(corr_above_0.6 = ifelse(abs(correlation) > 0.6, T, F))

# export final list
#write_tsv(ind_nutrient_COGs, "ind_nutrient_COGs.tsv")




######### DIAGNOSTICS #########



table(ind_nutrient_COGs$corr_above_0.6)
table(COG_df$corr_above_0.6)

table(ind_nutrient_COGs$corr_above_0.6, ind_nutrient_COGs$metadata)
#prop.test(t(table(ind_nutrient_COGs$corr_above_0.6, ind_nutrient_COGs$metadata)))

table(ind_nutrient_COGs$polarity)

ind_nutrient_COGs |> filter(COG == "COG0787")

ind_nutrient_COGs |> filter(COG == "COG5013") # no nitrogen???
temp.lm <- lm(data = tara_md, paste0("COG5013 ~ ", 
                                     str_flatten(nonnutrients, collapse = " + "), 
                                     " + ",
                                     "log_NO2NO3"))
tidy(temp.lm)
car::vif(temp.lm)
tidy(lm(data = tara_md, paste0("COG5013 ~ ", 
                               str_flatten(nonnutrients, collapse = " + "), 
                               " + ",
                               "log_NO2NO3",
                               " - Mean_Temperature"))) 
# so, the nitrate reductase COG is significant only when Temperature isn't controlled for..

cor(tara_md$Mean_Temperature, tara_md$log_NO2NO3, use = "complete.obs")

ggplot(data = tara_md, aes(x = Mean_Temperature, y = log_NO2NO3,
                           color = layer)) +
  geom_point()

nrow(tara_md)
nrow(tara_md |> drop_na())

tara_md[which(is.na(tara_md$log_NO2NO3)),]


# compare COG_df and ind_nutrient_COGs
COG_df_onlynutrients <- COG_df |> 
  filter(metadata %in% nutrients)

by <- c("COG", "metadata")

both <- inner_join(COG_df_onlynutrients, ind_nutrient_COGs, by = by)

nrow(COG_df_onlynutrients) - nrow(both) # so there's 856 correlations found in the first run but NOT the individual nutrient runs...?

first.not.individual <- anti_join(COG_df_onlynutrients, ind_nutrient_COGs, by = by)
table(first.not.individual$corr_above_0.6, first.not.individual$polarity) # mostly weaker correlations...

individual.not.first <- anti_join(ind_nutrient_COGs, COG_df_onlynutrients, by = by)
table(individual.not.first$corr_above_0.6, individual.not.first$polarity)

nrow(first.not.individual) # 856
nrow(both) # 1360
nrow(individual.not.first) # 2029
