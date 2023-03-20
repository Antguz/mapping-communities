################################################################################
##### 07 - Model testing
################################################################################

#' @description Batch for prediction of the models using the LSP scenes.

#-------------------------------------------------------------------------------
# Libraries

library(data.table)
library(caret)
library(parallel)
library(metrica)

#-------------------------------------------------------------------------------
# Root path
root_path <- "E:/Projects/DESIS_mapping"
meta <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv"))
ordination <- readRDS(paste0(root_path, "/05_ordination/dimensions_plots_taxonomic.rds"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra_VN.csv"))
models <- readRDS(paste0(root_path, "/06_model/model_final.rds"))
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

#-------------------------------------------------------------------------------
# Model Predict Function
model_performance <- function(meta = meta,
                              ordination = ordination, 
                              spectra = spectra,
                              type = "Teesting",
                              models = models) {
  
  #Performance descriptors
  my.metrics <- c("R2", "MAE", "RMAE", "RMSE", "RRMSE")
  
  #Number of models
  nMDS <- length(models)
  
  #Collectors
  final_performance <- data.table()
  final_predicted <- data.table()
  
  for(q in 1:nMDS) {
    
    #Repeated
    repeated <- length(models[[q]])
    
    #Observed
    observed <- ordination[, .SD, .SDcols = q]
    colnames(observed)[1] <- "MDS"
    observed <- observed$MDS
    
    
    #Collectors
    result_performance <- data.table()
    result_predicted <- data.table()
    
    for(i in 1:repeated) {
      
      grp <- models[[q]][[i]]
      predicted <- predict(grp, spectra)
      
      metrics <- metrics_summary(obs = observed, 
                                 pred = predicted, 
                                 type = "regression",
                                 metrics_list = my.metrics)
      metrics <- data.table(MDS = paste0("MDS", q),
                            repeated = i,
                            R2 = metrics[1,2],
                            MAE = metrics[2,2],
                            RMAE = metrics[3,2],
                            RMSE = metrics[4,2],
                            RRMSE = metrics[5,2])
      result_performance <- rbind(result_performance, metrics)
      
      predicted <- c(paste0("MDS", q), i, predicted)
      predicted <- as.data.table(matrix(predicted, nrow = 1))
      colnames(predicted)[1:2] <- c("MDS", "repeated")
      result_predicted <- rbind(result_predicted, predicted)
      
    }
    
    final_performance <- rbind(final_performance, result_performance)
    final_predicted <- rbind(final_predicted, result_predicted)
    
  }
  
  return(list(performance = final_performance,
              predicted = final_predicted))
  
}

testing <- model_performance(meta = meta[!split[,1]],
                             ordination = ordination[!split[,1]], 
                             spectra = spectra[!split[,1]],
                             type = "Testing",
                             models = models)
