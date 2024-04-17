
min_rmse <- models[type == "PRESS" & estimate == "CV" & Axes == "Axis3",]
plot(colMeans(min_rmse[type == "PRESS" & estimate == "CV" & Axes == "Axis3", 6:ncol(models)], na.rm = TRUE))
opt_min_rmse <- min_rmse[, .SD, .SDcols = c(6:ncol(min_rmse))]

pvalue <- NA

for(i in 1:(ncol(opt_min_rmse)-1)) {
  
  test <- t.test(as.matrix(opt_min_rmse)[,i], 
                 as.matrix(opt_min_rmse)[,i+1], 
                 paired = TRUE,
                 alternative = c("less"))
  
  pvalue[i+1] <- test$p.value
  
}

plot(1:length(pvalue), pvalue)
plot(colMeans(models[type == "PRESS" & estimate == "CV" & Axes == "Axis3", 6:ncol(models)], na.rm = TRUE))
plot(colMeans(models[type == "RMSEP" & estimate == "CV" & Axes == "Axis3", 6:ncol(models)], na.rm = TRUE))
plot(colMeans(models[type == "R2" & estimate == "train" & Axes == "Axis3", 6:ncol(models)], na.rm = TRUE))
plot(colMeans(models[type == "RMSEP" & estimate == "train" & Axes == "Axis3", 6:ncol(models)], na.rm = TRUE))

25, 25, 20
