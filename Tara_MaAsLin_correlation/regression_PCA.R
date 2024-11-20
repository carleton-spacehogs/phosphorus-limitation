library(tidyverse)

models <- readRDS("/workspace/data/Space_Hogs_shared_workspace/phosphorus-limitation/Tara_MaAsLin_correlation/MaAsLin_out/cplm-logCOGabun_against_no.ranefs_AbsLat_logPO4_sqrtIron_logNO2NO3_Salinity_Oxygen_temp_logDepth_size/fits/models.rds")

coefs <- numeric(length(models[[1]]@coefficients))
names(coefs) <- names(models[[1]]@coefficients)

coglist <- numeric()

for (i in 1:length(models)) {
  
  if (class(models[[i]]) != "logical") {
    coefs <- rbind(coefs, models[[i]]@coefficients)
    coglist <- append(coglist, names(models[i]))
  } else {
    # 39 NA models
  }
}

coefs <- coefs[-1,] # remove row of zeros

pca <- prcomp(coefs, scale = T, center = T)

ggplot(pca$x, aes(PC1, PC2)) +
  geom_point(alpha = 0.2)
# womp womp

screeplot(pca, type = "l", npcs = 21, main = "Screeplot of  PCs")
