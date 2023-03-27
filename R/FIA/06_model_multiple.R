#'------------------------------------------------------------------------------
#' @title Integration of PLSR and NMDS
#'------------------------------------------------------------------------------

#' @description It model the of species PLSR with NMDS.
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(caret)
library(pls)
library(plsVarSel)
library(sf)
library(parallel)
library(doSNOW)
library(metrica)
library(CAST)
library(viridis)
library(ggplot2)
library(rTLS)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param spectra
#' @param ordination
#' @param out_path

#root_path <- "E:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme Pro/Projects/DESIS_mapping"
meta <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv"))
ordination <- readRDS(paste0(root_path, "/05_ordination/dimensions_plots_taxonomic_DRYBIO.rds"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra_VN.csv"))
PD <- fread(paste0(root_path, "/03_diversity/master_phylogenetic-diversity_DRYBIO.txt"))
out_path <- paste0(root_path, "/06_model")

#'------------------------------------------------------------------------------
#' Load source code

source("R/FIA/grid_selection.R") 
source("R/FIA/optimal_number_multiple.R") 
source("R/FIA/run_models_multiple.R")
source("R/FIA/performance.R")

#'------------------------------------------------------------------------------
#' @Analysis  

#Meta file
meta <- meta[, c("PLT_CN", "INVYR", "STATECD", "COUNTYCD", "LAT", "LON", "ELEV")]

#Merge phylogenetic diversity
PD <- PD[, c("PLT_CN", "ntaxa", "pd.obs", "PSR")]
PD <- merge(meta, 
            PD,
            by = "PLT_CN", all.x = FALSE, all.y = TRUE)

#Merge meta with ordination
ordination <- merge(PD, 
                    ordination,
                    by = "PLT_CN", all.x = FALSE, all.y = TRUE)

#Merge ordination with spectra 
data <- merge(ordination, 
              spectra[, .SD, .SDcols = c(1, 7:ncol(spectra))],
              by = "PLT_CN", all.x = FALSE, all.y = TRUE)

#New frames
meta <- data[, .SD, .SDcols = 1:10]
ordination <- data[, .SD, .SDcols = 11:13]
spectra <- data[, .SD, .SDcols = 15:ncol(data)]
weigths <- data$goodness

#------------------------------------------
# Data-split
set.seed(1107)
split <- createDataPartition(meta$COUNTYCD, 
                             p = 0.70,
                             list = FALSE, 
                             times = 1)
split <- split[,1]

#------------------------------------------
# Estimation of the optimal number of components

ncomponents <- optimal_number(meta = meta,
                              ordination = ordination,
                              spectra = spectra,
                              weigths = weigths,
                              split = split,
                              int_ncomp = 50, 
                              repeated = 100, 
                              threads = 26)

saveRDS(ncomponents, paste0(out_path, "/PRESS_taxonomic.csv"))

mean_PRESS <- ncomponents[, lapply(.SD, mean), by= MDS]

plot(as.numeric(mean_PRESS[2, 3:ncol(mean_PRESS)]))
plot(as.numeric(mean_PRESS[3, 3:19]))

ncomp <- 19

#------------------------------------------
# Run models

models <- run_models_multiple(meta, 
                              ordination, 
                              spectra, 
                              split, 
                              int_ncomp = ncomp, 
                              repeated = 100,
                              threads = 26)

saveRDS(models, paste0(out_path, "/models_taxonomic.rds"))

#------------------------------------------
# Model performance

model <- model_performance_multiple(meta = meta,
                                    ordination = ordination, 
                                    spectra = spectra,
                                    models = models,
                                    ncomp = 19)

saveRDS(final_performance, paste0(out_path, "/models_performance.rds"))
saveRDS(final_predicted, paste0(out_path, "/models_predicted.rds"))

mean_summary <- final_predicted[, lapply(.SD, mean, na.rm=TRUE), by= MDS, .SDcols= 3:ncol(final_predicted)] 
sd_summary <- final_predicted[, lapply(.SD, sd, na.rm=TRUE), by= MDS, .SDcols= 3:ncol(final_predicted)] 

mean_summary <- melt(mean_summary, id.vars = c("MDS"), variable.name = "ID")
sd_summary <- melt(sd_summary, id.vars = c("MDS"), variable.name = "ID")
mean_summary$ID <- as.numeric(gsub('V', '', mean_summary$ID))
sd_summary$ID <- as.numeric(gsub('V', '', sd_summary$ID))
mean_summary <- dcast(mean_summary, ID ~ MDS)
sd_summary <- dcast(sd_summary, ID ~ MDS)

plot(ordination$MDS1, sd_summary$MDS1)
abline(lm(ordination$MDS3 ~ mean_summary$MDS3))
summary(lm(ordination$MDS3 ~ mean_summary$MDS3))

delta <- (mean_summary$MDS3 - ordination$MDS3) / ordination$MDS3

#------------------------------------------
# Evaluate performance

model$Performance

ggplot(model$Performance, aes(x= factor(MDS), y= as.numeric(Rsq), fill = factor(MDS))) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) 

saveRDS(model, paste0(root_path, "/06_model/model_cwt.rds"))
model <- readRDS(paste0(root_path, "/06_model/model_cwt.rds"))



evaluation <- model_performance(meta, ordination, spectra, split, model)
saveRDS(evaluation, paste0(root_path, "/06_model/evaluation_cwt.rds"))

ggplot(evaluation$Performance, aes(x= factor(MDS), y= as.numeric(R2), fill = factor(MDS))) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) 


