#'------------------------------------------------------------------------------
#' @title Codes to visualize the predicted values
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(ggplot2)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Load source code

source("R/data_split.R")

#'------------------------------------------------------------------------------
#' Read information

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme Pro/Projects/DESIS_mapping"
meta <- fread(paste0(root_path, "/00_FIA/plots_WI_MN_20221102.csv"))
ordination <- readRDS(paste0(root_path, "/05_ordination/dimensions_plots_beta_DRYBIO.rds"))
spectra <- fread(paste0(root_path, "/04_match_spectra-stand/spectra_SW.csv"))
PD <- fread(paste0(root_path, "/03_diversity/master_phylogenetic-diversity_DRYBIO.txt"))

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

# New frames
data <- na.exclude(data)
meta <- data[, .SD, .SDcols = 1:10]
ordination <- data[, .SD, .SDcols = 11:13]
spectra <- data[, .SD, .SDcols = 15:ncol(data)]
weigths <- data$goodness

#------------------------------------------
# Data-split
set.seed(1107)
split <- data_split(ordination)

#------------------------------------------
# Read predicted values

performance <- fread(paste0(root_path, "/06_model/SW_beta_PLS/performance.csv"))
predicted_training <- fread(paste0(root_path, "/06_model/SW_beta_PLS/predicted_training.csv"))
predicted_testing <- fread(paste0(root_path, "/06_model/SW_beta_PLS/predicted_testing.csv"))

#------------------------------------------
# Get summary of predicted

predicted_melt_testing <- melt(predicted_testing, id.vars = c("MDS", "repeated", "type"), 
                               variable.name = "sample",
                               value.name = "predicted")

predicted_summary_testing <- predicted_melt_testing[, .(predicted = mean(predicted),
                                                        predicted_sd = sd(predicted)),
                                                    by = c("MDS", "sample")]
predicted_summary_testing$sample <- as.character(predicted_summary_testing$sample)

ordination_testing <- ordination[!split, ]
ordination_testing$sample <- paste0("V", 3:(nrow(ordination_testing)+2))
ordination_testing_melt <- melt(ordination_testing, id.vars = c("sample"), 
                                variable.name = "MDS",
                                value.name = "observed")

#Merge
frame <- merge(ordination_testing_melt, predicted_summary_testing, by = c("sample", "MDS"))

#------------------------------------------
# Predict slopes and intercepts

values <- seq(-0.6, 0.6, by = 0.001)

#------------------------------------------
# Do plot

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th <- theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(frame[MDS == "MDS1"]) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_errorbar(aes(xmin = predicted - predicted_sd, xmax = predicted + predicted_sd, y = observed), colour = "grey") +
  geom_point(aes(x = predicted, y = observed), size = 2, shape = 21, fill = col[1], colour = "white") +
  #coord_cartesian(xlim = c(-0.55, 0.35), ylim = c(-0.55, 0.35)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted (MDS1)") + ylab("Observed (MDS1)") +
  th

ggplot(frame[MDS == "MDS2"]) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_errorbar(aes(xmin = predicted - predicted_sd, xmax = predicted + predicted_sd, y = observed), colour = "grey") +
  geom_point(aes(x = predicted, y = observed), size = 2, shape = 21, fill = col[2], colour = "white") +
  #coord_cartesian(xlim = c(-0.30, 0.30), ylim = c(-0.30, 0.30)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted (MDS2)") + ylab("Observed (MDS2)") +
  th

ggplot(frame[MDS == "MDS3"]) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_errorbar(aes(xmin = predicted - predicted_sd, xmax = predicted + predicted_sd, y = observed), colour = "grey") +
  geom_point(aes(x = predicted, y = observed), size = 2, shape = 21, fill = col[3], colour = "white") +
  #coord_cartesian(xlim = c(-0.35, 0.30), ylim = c(-0.35, 0.30)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted (MDS3)") + ylab("Observed (MDS3)") +
  th
  

new_ordination <- dcast(frame, sample ~ MDS, value.var = "predicted")
old_ordination <- dcast(frame, sample ~ MDS, value.var = "observed")
sd_ordination <- dcast(frame, sample ~ MDS, value.var = "predicted_sd")
sd_ordination <- rowMeans(sd_ordination[, 2:4])
plot(new_ordination$MDS2, new_ordination$MDS3)
plot(ordination_testing$MDS2, ordination_testing$MDS3)

plot3d(new_ordination[, 2:4])
plot3d(old_ordination[, 2:4])

nrow(new_ordination)

for(i in 1:nrow(new_ordination)) {
  arrow3d(old_ordination[i, 2:4], new_ordination[i, 2:4], 
          type = "rotation", 
          col = "grey", 
          width = 0.25,
          thickness = 0.2,
          n = 10,
          s = 0.25,
          theta = pi/128)
}
