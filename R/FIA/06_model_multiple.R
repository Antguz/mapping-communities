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
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param spectra
#' @param ordination
#' @param out_path

#root_path <- "E:/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
root_path <- "/media/antonio/Extreme Pro/Projects/DESIS_mapping"
meta <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv"))
ordination <- readRDS(paste0(root_path, "/05_ordination/dimensions_plots_taxonomic_DRYBIO.rds"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra_VN.csv"))
PD <- fread(paste0(root_path, "/03_diversity/master_phylogenetic-diversity_DRYBIO.txt"))
out_path <- paste0(root_path, "/06_model")

#'------------------------------------------------------------------------------
#' Load source code

source("R/FIA/optimal_number.R") 
source("R/FIA/pls-ordination_model.R")
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
                              repeated = 4, 
                              threads = 4)

#------------------------------------------
# Final Model of prediction

model <- pls_ordination_model(meta = meta,
                              ordination = ordination, 
                              spectra = spectra, 
                              split = split, 
                              repeated = 5,
                              opt_components = c(13, 5, 2),
                              threads = 28)

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

#------------------------------------------
# Final Model of prediction

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


