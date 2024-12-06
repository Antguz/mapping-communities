#'------------------------------------------------------------------------------
#' @title Visualize observed and predicted occurrence
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(ggplot2)
library(ggpubr)
library(bit64)
library(ggdark)
library(RcppArmadillo)
library(Rcpp)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/probability_rcpp.cpp")

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/MMDB"
type <- "phylo"

# ------------------------------------------------------------------------------
# Read CWM from observed and predicted models

# EMIT ---------------------------------------
sensor <- "EMIT"
training <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                         sensor, "_", type, "_predicted_", "training", ".csv"))

testing <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                        sensor, "_", type, "_predicted_", "testing", ".csv"))
testing$sample <- testing$sample * 10000
dataset <- rbind(training, testing)

# Reshape predicted
# Select columns
obs <- dataset[, c("Axes", "sample", "type", "observed")]
pred <- dataset[, c("Axes", "sample", "type", "predicted")]

# Dcast columns
EMIT_obs <- dcast(obs, sample + type ~ Axes)
EMIT_pred <- dcast(pred, sample + type ~ Axes)

# DESIS ---------------------------------------
sensor <- "DESIS"
training <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                         sensor, "_", type, "_predicted_", "training", ".csv"))

testing <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                        sensor, "_", type, "_predicted_", "testing", ".csv"))
testing$sample <- testing$sample * 10000
dataset <- rbind(training, testing)

# Reshape predicted
# Select columns
obs <- dataset[, c("Axes", "sample", "type", "observed")]
pred <- dataset[, c("Axes", "sample", "type", "predicted")]

# Dcast columns
DESIS_obs <- dcast(obs, sample + type ~ Axes)
DESIS_pred <- dcast(pred, sample + type ~ Axes)

# ------------------------------------------------------------------------------
# Get and reshape coefficients
# Get glm coefficients
coefficients <- fread(paste0(root_path, "/Models/03_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

# Reshape coefficients
models <- unique(coefficients$plant)
coef_array <- array(dim = c(100, 3, length(models)))
intercept_mat <- matrix(nrow = 100, ncol = length(models))

for(i in 1:length(models)) {
  
  sub <- coefficients[plant == models[i]]
  coef_array[,,i] <- as.matrix(sub[, 3:5])
  intercept_mat[,i] <- sub$`(Intercept)`
  
}

# ------------------------------------------------------------------------------
# Predict

# EMIT
EMIT_traits_observed <- summary_prob(as.matrix(EMIT_obs[, 3:5]),
                                     coef_array,
                                     intercept_mat)[[1]]
EMIT_traits_observed <- as.data.table(EMIT_traits_observed)
colnames(EMIT_traits_observed) <- models

EMIT_traits_predicted <- summary_prob(as.matrix(EMIT_pred[, 3:5]),
                                      coef_array,
                                      intercept_mat)[[1]]
EMIT_traits_predicted <- as.data.table(EMIT_traits_predicted)
colnames(EMIT_traits_predicted) <- models

# DESIS
DESIS_traits_observed <- summary_prob(as.matrix(DESIS_obs[, 3:5]),
                                      coef_array,
                                      intercept_mat)[[1]]
DESIS_traits_observed <- as.data.table(DESIS_traits_observed)
colnames(DESIS_traits_observed) <- models

DESIS_traits_predicted <- summary_prob(as.matrix(DESIS_pred[, 3:5]),
                                       coef_array,
                                       intercept_mat)[[1]]
DESIS_traits_predicted <- as.data.table(DESIS_traits_predicted)
colnames(DESIS_traits_predicted) <- models

# plot(EMIT_traits_observed$Quercus ~ EMIT_traits_predicted$Quercus)
# plot(DESIS_traits_observed$drought ~ DESIS_traits_predicted$drought)
# 
# 
# traits_predicted <- EMIT_traits_predicted$`White oak / red oak / hickory`
# traits_observed <- EMIT_traits_observed$`White oak / red oak / hickory`
# type <- EMIT_obs$type
# fills <- c("#5e3c99", "transparent")
# fills <- c("#e66101", "transparent")

# ------------------------------------------------------------------------------
# Plot

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position.inside = c(0.90, 0.90), 
        legend.position = "top",
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #legend.position.inside = c(0.90, 0.90), 
        legend.position = "top",
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth=0),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"))

th <- th_trans
line_col <- "white"

th <- th_black
line_col <- "black"


plot_trait <- function(traits_predicted,
                       traits_observed,
                       type,
                       fills,
                       shape = c(21, 21),
                       lim = c(0, 1)) {
  
  frame <- data.table(traits_predicted = traits_predicted/10000, 
                      traits_observed = traits_observed/10000, 
                      type = type)
  
  frame$type <- factor(frame$type, levels = c("Training", "Testing"))
  
  p <- ggplot(frame) +
    geom_point(aes(x = traits_predicted,
                   y = traits_observed,
                   shape = type,
                   color = type,
                   fill = type),
               alpha = 0.75,
               size = 0.1) +
    geom_smooth(aes(x = traits_predicted,
                    y = traits_observed,
                    linetype = type),
                method = lm, 
                se = FALSE,
                linewidth = 0.4,
                color = line_col) +
    scale_shape_manual("Dataset", values = shape) +
    scale_fill_manual("Dataset", values = fills) + 
    scale_color_manual("Dataset", values = c(fills[1], fills[1])) +
    scale_linetype_manual("Dataset", values = c("dotted", "dashed")) +
    coord_cartesian(xlim = lim, ylim = lim) +
    scale_x_continuous(n.breaks = 4) +
    scale_y_continuous(n.breaks = 4) +
    xlab("Predicted (from satellite)") +
    ylab("Observed (from FIA)") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    th
  
  return(p)
  
}

EMIT <- plot_trait(traits_predicted = EMIT_traits_predicted$Quercus,
                   traits_observed = EMIT_traits_observed$Quercus,
                   type = EMIT_obs$type,
                   fills = c("#5e3c99", "transparent"),
                   shape = c(21, 21),
                   lim = c(0, 1))

DESIS <- plot_trait(traits_predicted = DESIS_traits_predicted$Quercus,
                    traits_observed = DESIS_traits_observed$Quercus,
                    type = DESIS_obs$type,
                    fills = c("#e66101", "transparent"),
                    shape = c(24, 24),
                    lim = c(0, 1))

rel <- ggarrange(DESIS, EMIT, 
                 ncol=2, nrow=1,
                 labels = c("", ""),
                 font.label =  list(size = 9, face = "plain", color = line_col),
                 label.x = -0.05,
                 label.y = 0.95)

# Export JPEG
jpeg(paste0(root_path, "/figures/Manuscript_black/", "occurrence-model_observed-predicted.jpeg"), 
     width = 150, 
     height = 60, 
     units = "mm", 
     res = 300, 
     bg = "white")

rel

dev.off()

# PNG

png(paste0(root_path, "/figures/Manuscript_white/", "occurrence-model_observed-predicted.png"), 
    width = 110, 
    height = 60, 
    units = "mm", 
    res = 300, 
    bg = "transparent")

rel

dev.off()

