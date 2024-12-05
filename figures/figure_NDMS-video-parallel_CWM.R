#'------------------------------------------------------------------------------
#' @title Visualize the NDMS as a video
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(gg3D)
library(stringr)
library(bit64)
library(plotly)
library(plot3D)
library(data.table)
library(ggpubr)
library(Rcpp)
library(RcppArmadillo)
library(pbmcapply)
library(viridis)
library(scales)
library(vegan)

#'------------------------------------------------------------------------------
#' @Source_code
#-------------------------------------------------------------------------------

sourceCpp("src/cwm_rcpp.cpp")
source("R/auxiliary/CWM.R")

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"

#'------------------------------------------------------------------------------
#' Func ordination

type <- "func"

# Communities
func_comm <- fread(paste0(root_path, "/Diversity/02_NMDS/NMDS_scores_comm_", type, "_total.csv"))

# ------------------------------------------------------------------------------
# Func observed

# Meta information
fia_meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))

# FIA data for aggregation
inventory <- readRDS(paste0(root_path, "/Diversity/01_beta/master_matched_comm_", type, "_VOL.rds"))
sp_meta <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))
match <- match(colnames(inventory), sp_meta$Taxon2)
sp_meta <- sp_meta[match, ]
aggregation <- CWM_batch(sp_meta[,c("Taxon2", "HTq", "SLq", "WD", "shade", "drought", "waterlogging", "LMA", "N")], 
                         decostand(inventory, "total"))

# Test for match
all.equal(aggregation$PLT_CN, func_comm$PLT_CN)

# ------------------------------------------------------------------------------
# Func modelled

# Get glm coefficients
coefficients <- fread(paste0(root_path, "/Models/03_model_CWM/", 
                             type, "_model-CWM.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 11:14)]

# Reshape coefficients
traits <- unique(coefficients$plant)
coef_array <- array(dim = c(100, 3, length(traits)))
intercept_mat <- matrix(nrow = 100, ncol = length(traits))

for(i in 1:length(traits)) {
  
  sub <- coefficients[plant == traits[i]]
  coef_array[,,i] <- as.matrix(sub[, 3:5])
  intercept_mat[,i] <- sub$`(Intercept)`
  
}

func_modelled <- summary_cwm(as.matrix(func_comm[, c(-1, -5)]),
                             coef_array, 
                             intercept_mat)[[1]]
func_modelled <- as.data.table(func_modelled)
colnames(func_modelled) <- traits

#'------------------------------------------------------------------------------
#' @Plot-CWM
#-------------------------------------------------------------------------------

# Plot features
theta <- c(seq(0, 180, length.out = 180), seq(-180, 0, length.out = 180))
phi <- c(seq(0, 0, length.out = length(theta)))
h_adjust_axis1 <- c(seq(1, 0, length.out = 180), seq(0, 1, length.out = 180))
h_adjust_axis2 <- c(seq(1, 1, length.out = 180), #0-180
                    seq(1, 0, length.out = 45), #180-225
                    seq(0, 0, length.out = 45), #225-270
                    seq(0, 0, length.out = 60),
                    seq(0, 1, length.out = 30)) #270-360

th <- theme_void(base_size = 15) +
  theme(panel.background = element_rect(fill = "black",
                                        colour = "black",
                                        linewidth = 0.5, 
                                        linetype = "solid"),
        plot.background = element_rect(fill = "black"),
        legend.title = element_text(colour="white", size = 20),
        legend.text = element_text(color = "white", size = 15),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0), 'pt'),
        legend.title.align=0.5)

# Vizual plot

application <- function(X,
                        func_comm, 
                        func_observed_trait,
                        func_modelled_trait,
                        legend, 
                        pallete,
                        lim = c(1, 5),
                        export, 
                        class,
                        theta, 
                        phi, 
                        h_adjust_axis1,
                        h_adjust_axis2,
                        th) {
  
  i <- X
  
  func_2 <- ggplot(func_comm, aes(x = Axis1, 
                                  y = Axis2, 
                                  z = Axis3,
                                  color = func_observed_trait)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_colour_gradientn(legend, colours = pallete, limits = lim, oob = scales::squish) +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    guides(colour = guide_colourbar(barwidth = 20, barheight = 0.7, title.position = "bottom")) +
    th 
  
  func_3 <- ggplot(func_comm, aes(x = Axis1, 
                                  y = Axis2, 
                                  z = Axis3,
                                  color = func_modelled_trait)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_colour_gradientn(legend, colours = pallete, limits = lim, oob = scales::squish) +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    guides(colour = guide_colourbar(barwidth = 20, barheight = 0.7, title.position = "bottom")) +
    th 
  
  name <- str_pad(i, 4, pad = "0")
  file_name = paste0(export, "/", class, "_", name, ".png", sep="")
  png(file_name,
      width = 400, 
      height = 210, 
      units = "mm", 
      pointsize = 32,
      bg = NA,  
      res = 300)
  print(ggarrange(func_2, func_3, ncol = 2, common.legend = TRUE,
                  legend = "bottom"))
  dev.off()
  
  return(NULL)
  
}

# CWM of HTq tolerance index
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           func_comm = func_comm, 
           func_observed_trait = aggregation$HTq,
           func_modelled_trait = func_modelled$HTq,
           legend = "CWM of maximun tree hight", 
           pallete = viridis(1000),
           lim = c(7, 35),
           export = "/media/antonio/Work/MMDB/figures/video/CWM/HTq", 
           class = "HTq",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 25,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

# CWM of SLq tolerance index
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           func_comm = func_comm, 
           func_observed_trait = aggregation$SLq,
           func_modelled_trait = func_modelled$SLq,
           legend = "CWM of tree slenderness", 
           pallete = viridis(1000),
           lim = c(60, 110),
           export = "/media/antonio/Work/MMDB/figures/video/CWM/SLq", 
           class = "SLq",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 25,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

# CWM of WD tolerance index
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           func_comm = func_comm, 
           func_observed_trait = aggregation$WD,
           func_modelled_trait = func_modelled$WD,
           legend = "CWM of wood density", 
           pallete = viridis(1000),
           lim = c(0.3, 0.7),
           export = "/media/antonio/Work/MMDB/figures/video/CWM/WD", 
           class = "WD",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 25,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

# CWM of shade tolerance index
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           func_comm = func_comm, 
           func_observed_trait = aggregation$shade,
           func_modelled_trait = func_modelled$shade,
           legend = "CWM of shade tolerance index", 
           pallete = viridis(1000),
           lim = c(1, 5),
           export = "/media/antonio/Work/MMDB/figures/video/CWM/shade", 
           class = "shade",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 25,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

# CWM of drought tolerance index
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           func_comm = func_comm, 
           func_observed_trait = aggregation$drought,
           func_modelled_trait = func_modelled$drought,
           legend = "CWM of drought tolerance index", 
           pallete = viridis(1000),
           lim = c(1, 5),
           export = "/media/antonio/Work/MMDB/figures/video/CWM/drought", 
           class = "drought",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 25,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

# CWM of waterlogging tolerance index
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           func_comm = func_comm, 
           func_observed_trait = aggregation$waterlogging,
           func_modelled_trait = func_modelled$waterlogging,
           legend = "CWM of Waterlogging tolerance index", 
           pallete = viridis(1000),
           lim = c(1, 5),
           export = "/media/antonio/Work/MMDB/figures/video/CWM/waterlogging", 
           class = "waterlogging",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 25,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

# CWM of LMA
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           func_comm = func_comm, 
           func_observed_trait = aggregation$LMA,
           func_modelled_trait = func_modelled$LMA,
           legend = "CWM of leaf mass area", 
           pallete = viridis(1000),
           lim = c(50, 300),
           export = "/media/antonio/Work/MMDB/figures/video/CWM/LMA", 
           class = "LMA",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 25,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

# CWM of LMA
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           func_comm = func_comm, 
           func_observed_trait = aggregation$N,
           func_modelled_trait = func_modelled$N,
           legend = "CWM of nitrogen", 
           pallete = viridis(1000),
           lim = c(5, 30),
           export = "/media/antonio/Work/MMDB/figures/video/CWM/N", 
           class = "N",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 25,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)

#-------------------------------------------------------------------------------
# CWM models

type <- "func"

# ------------------------------------------------------------------------------
# Read CWM from observed and predicted models
sensor <- "DESIS"
training <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                         sensor, "_", type, "_predicted_", "training", ".csv"))

testing <- fread(paste0(root_path, "/Models/01_model_beta/", sensor, "/",
                        sensor, "_", type, "_predicted_", "testing", ".csv"))
testing$sample <- testing$sample * 10000
dataset <- rbind(training, testing)

# Reshape predicted
# Select columns
obs <- dataset[, c("Axes", "sample", "observed")]
pred <- dataset[, c("Axes", "sample", "predicted")]

# Dcast columns
EMIT_obs <- dcast(obs, sample ~ Axes)
EMIT_pred <- dcast(pred, sample ~ Axes)

# ------------------------------------------------------------------------------
# Get and reshape coefficients
# Get glm coefficients
coefficients <- fread(paste0(root_path, "/Models/03_model_CWM/",
                             type, "_model-CWM.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 11:14)]

# Reshape coefficients
traits <- unique(coefficients$plant)
coef_array <- array(dim = c(100, 3, length(traits)))
intercept_mat <- matrix(nrow = 100, ncol = length(traits))

for(i in 1:length(traits)) {

 sub <- coefficients[plant == traits[i]]
 coef_array[,,i] <- as.matrix(sub[, 3:5])
 intercept_mat[,i] <- sub$`(Intercept)`

}

traits_observed <- summary_cwm(as.matrix(EMIT_obs[, -1]),
                               coef_array,
                               intercept_mat)[[1]]
traits_observed <- as.data.table(traits_observed)
colnames(traits_observed) <- traits

traits_predicted <- summary_cwm(as.matrix(EMIT_pred[, -1]),
                                coef_array,
                                intercept_mat)[[1]]
traits_predicted <- as.data.table(traits_predicted)
colnames(traits_predicted) <- traits


# CWM of drought tolerance index
pbmclapply(X = 1:length(theta),  #length(theta)
           FUN = application, 
           func_comm = EMIT_pred[, -1], 
           func_observed_trait = traits_observed$drought,
           func_modelled_trait = traits_predicted$drought,
           legend = "CWM of drought tolerance index", 
           pallete = viridis(1000),
           lim = c(1, 5),
           export = "/media/antonio/Work/MMDB/figures/video/CWM/drought_test", 
           class = "drought_test",
           theta = theta, 
           phi = phi, 
           h_adjust_axis1 = h_adjust_axis1,
           h_adjust_axis2 = h_adjust_axis2,
           th = th,
           mc.preschedule = TRUE, 
           mc.set.seed = FALSE,
           mc.silent = FALSE, 
           mc.cores = 25,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE, 
           affinity.list = NULL)
