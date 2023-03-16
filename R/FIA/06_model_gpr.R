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
library(blockCV)
library(sf)
library(parallel)
library(doSNOW)
library(metrica)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param spectra
#' @param ordination
#' @param out_path

#root_path <- "E:/DESIS"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
meta <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv"))
ordination <- readRDS(paste0(root_path, "/05_ordination/dimensions_plots_taxonomic.rds"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra_VN.csv"))
out_path <- paste0(root_path, "/06_model")

#'------------------------------------------------------------------------------
#' @Analysis  

#Meta file
meta <- meta[, c("PLT_CN", "INVYR", "STATECD", "COUNTYCD", "LAT", "LON", "ELEV")]

#Merge meta with ordination
ordination <- merge(meta, 
                    ordination,
                    by = "PLT_CN", all.x = FALSE, all.y = TRUE)

#Merge ordination with spectra 
data <- merge(ordination, 
              spectra[, .SD, .SDcols = c(1, 7:ncol(spectra))],
              by = "PLT_CN", all.x = FALSE, all.y = TRUE)

#New frames
meta <- data[, .SD, .SDcols = 1:7]
ordination <- data[, .SD, .SDcols = 8:10]
spectra <- data[, .SD, .SDcols = 11:ncol(data)]

#------------------------------------------
# Data-split
set.seed(1107)
split <- createDataPartition(meta$PLT_CN, 
                             p = 0.70,
                             list = FALSE, 
                             times = 1)

#Create folds and split -----------------------
#Model collectors
data_split <- list()
data_folds <- list()

#N of iterations
repeated <- 100

#Proportion of data split
prop <- 0.8

for(i in 1:100) {
  
  #Split
  split <- datasplit_tile(meta[split,], prop)
  
  #Folds
  folds <- CreateSpacetimeFolds(x = condition_training[split],
                                spacevar = "tile", 
                                k = 8,
                                seed = sample(1:10000, 1))
  
  #Collect
  data_split[[i]] <- split
  data_folds[[i]] <- folds
  
  #Residuals
  rm(list = c("split", "folds"))
  
}

ctrl <- trainControl(method = "cv", repeats = 3)

model_train <- train(Class ~ .,
                     data = training,
                     method = "gaussprLinear",
                     preProc = c("center", "scale"),
                     trControl = ctrl)

ncomponents <- optimal_number(meta = meta,
                              ordination = ordination,
                              spectra = spectra,
                              split = split,
                              int_ncomp = 60, 
                              repeated = 50, 
                              threads = 28)

saveRDS(ncomponents, paste0(root_path, "/06_model/ncomponents_cwt.rds"))

hist(ncomponents[ordination == "MDS1"]$CV)
hist(ncomponents[ordination == "MDS3"]$Onesigma)
hist(ncomponents[ordination == "MDS1"]$PRESS)


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

ncomponents[, getmode(PRESS), by = "ordination"]
ncomponents[, mean(PRESS), by = "ordination"]

#------------------------------------------
# Final Model of prediction

model <- pls_ordination_model(meta = meta,
                              ordination = ordination, 
                              spectra = spectra, 
                              split = split, 
                              repeated = 100,
                              opt_components = c(33, 9, 7),
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

ggplot(evaluation$Performance, aes(x= factor(MDS), y= as.numeric(RMSE), fill = factor(MDS))) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) 


