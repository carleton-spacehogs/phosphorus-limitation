# Author: Troy Osborn
# Date: Jan 2025

# Create a table showing all the COGs and whether they had significant 
# correlations with each nutrient. For each COG and each nutrient, include
# a column for positive/negative/no correlation and a column for q-value

library(tidyverse)

results <- read_tsv("MaAsLin_out/cplm-logCOGabun_against_no.ranefs_AbsLat_logPO4_sqrtIron_logNO2NO3_Salinity_Oxygen_temp_logDepth_size/all_results.tsv")

# assign "positive", "negative", or "not_significant" based on the q-value
# and the sign of the coefficient by splitting and rejoining
sig_cors <- results |> 
  filter(qval < 0.05) |> 
  mutate(direction = ifelse(coef >= 0, "positive", "negative"))
nonsig_cors <- results |>
  filter(qval >= 0.05) |> 
  mutate(direction = "not_significant")

joined <- rbind(sig_cors, nonsig_cors) |> 
  select(feature, metadata, qval, direction) |> 
  arrange(feature)

# pivot to a "wide" dataframe (one row per COG, instead of multiple), once
# again by splitting and rejoining
pivoted.metadata <- joined |> 
  select(-qval) |> 
  pivot_wider(names_from = metadata,
              values_from = direction)
pivoted.qvals <- joined |> 
  select(-direction) |> 
  pivot_wider(names_from = metadata,
              values_from = qval)

colnames(pivoted.qvals) <- paste0(colnames(pivoted.qvals), "_qval") 
colnames(pivoted.qvals)[[1]] <- "feature"

pivoted.joined <- left_join(pivoted.metadata, pivoted.qvals, by = "feature")


# reorder columns to pair metadata columns with their q-value columns
intended.order <- c("Mean_Oxygen", "log_PO4", "log_NO2NO3", "sqrt_iron", 
                    "log_depth", "Mean_Temperature", "Mean_Salinity", "size_fraction", "Absolute_Latitude")
pivoted.joined <- pivoted.joined |> 
  select(1, 2, 11, 4, 13, 10, 19, 3, 12, 9, 18, 8, 17, 6, 15, 5, 14, 7, 16)

# save final table as csv
write_csv(pivoted.joined, "tables/COG_correlation_table.csv")



### Make figure (not table) of rho values among significant correlations

sigList <- read_tsv("significant_COG_list.tsv")
ggplot(sigList, aes(x = correlation)) +
  geom_density() +
  xlim(c(-1, 1)) +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed") +
  labs(y = "Density", x = "Rho", title = "Rho values of significant correlations")

ggsave("graphs/RhoDensity.png", width = 6, height = 4)
