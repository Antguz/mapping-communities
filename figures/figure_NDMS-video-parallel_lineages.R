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

sourceCpp("src/probability_rcpp.cpp")
source("R/auxiliary/aggregate_plants.R")

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"

#'------------------------------------------------------------------------------
#' Func ordination

type <- "phylo"

# Communities
phylo_comm <- fread(paste0(root_path, "/Diversity/02_NMDS/NMDS_scores_comm_", type, "_total.csv"))

# ------------------------------------------------------------------------------
# Phylo observed

# Meta information
fia_meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))

# FIA data for aggregation
inventory <- readRDS(paste0(root_path, "/Diversity/01_beta/master_matched_comm_", type, "_VOL.rds"))
sp_meta <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))
aggregation <- aggregate_batch(sp_meta[,c("Taxon2", "Family", "Genus", "growth_form", 
                                          "mycorrhizal_association", "leaf_shape", 
                                          "leaf_type", "leaf_habit")], 
                               decostand(inventory, "total"))
#aggregation <- aggregation[, 1:2]

# Test for match
all.equal(aggregation$PLT_CN, phylo_comm$PLT_CN)
fia_meta_match <- match.integer64(phylo_comm$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match,]
all.equal(aggregation$PLT_CN, fia_meta$PLT_CN)

# Select the most abundant groups
aggregation <- aggregation[, names(aggregation)[-1] := lapply(.SD, function(x) as.integer(x!=0)), 
                           .SDcols = 2:ncol(aggregation)]
presence <- colSums(aggregation[,-1])
rpresence <- round(presence/nrow(aggregation)*100, 4)
rpresence <- rpresence > 10.0
rpresence <- c(TRUE, rpresence)
aggregation <- aggregation[ , .SD, .SDcols = c(1:ncol(aggregation))[rpresence]]

# ------------------------------------------------------------------------------
# Phylo modelled

# Get glm coefficients
coefficients <- fread(paste0(root_path, "/Models/03_model_occurrence/", 
                             type, "_model-occurrence.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

# Reshape coefficients
plant <- unique(coefficients$plant)
coef_array <- array(dim = c(100, 3, length(plant)))
intercept_mat <- matrix(nrow = 100, ncol = length(plant))

for(i in 1:length(plant)) {
  
  sub <- coefficients[plant == plant[i]]
  coef_array[,,i] <- as.matrix(sub[, 3:5])
  intercept_mat[,i] <- sub$`(Intercept)`
  
}

phylo_modelled <- summary_prob(as.matrix(phylo_comm[, c(-1, -5)]),
                               coef_array, 
                               intercept_mat)[[1]]
phylo_modelled <- as.data.table(phylo_modelled)/10000
colnames(phylo_modelled) <- plant

#'------------------------------------------------------------------------------
#' @Plot-Lineages
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
                        phylo_comm, 
                        phylo_observed_plant,
                        phylo_modelled_plant,
                        legend, 
                        pallete,
                        export, 
                        class,
                        theta, 
                        phi, 
                        h_adjust_axis1,
                        h_adjust_axis2,
                        th) {
  
  i <- X
  
  Occurrence <- ifelse(phylo_observed_plant == 1, "Presence", "Absence")
  Occurrence <- factor(Occurrence, levels = c("Presence", "Absence"))
  
  phylo_2 <- ggplot(phylo_comm, aes(x = Axis1, 
                                  y = Axis2, 
                                  z = Axis3,
                                  color = Occurrence)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_color_manual(values = c("red", "grey")) +
    labs_3D(theta = theta[i], 
            phi = phi[i], 
            labs = c("Axis 1", "Axis 2", "Axis 3"),
            hjust = c(h_adjust_axis1[i], h_adjust_axis2[i], 1), 
            vjust = c(-0.5, -0.5,-0.5), 
            angle = c(0, 0, 90),
            color = "white", 
            size = 6) +
    guides(color = guide_legend(override.aes = list(size = 3), title.position = "bottom")) +
    th
  
  phylo_3 <- ggplot(phylo_comm, aes(x = Axis1, 
                                  y = Axis2, 
                                  z = Axis3,
                                  color = phylo_modelled_plant)) +
    axes_3D(theta = theta[i], 
            phi = phi[i], 
            linewidth = 0.5,
            colour = "white") +
    stat_3D(theta = theta[i], 
            phi = phi[i], 
            shape = 19,
            size = 0.6,
            alpha = 0.6) +
    scale_colour_gradientn(legend, colours = pallete, limits = c(0, 1), oob = scales::squish) +
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
  print(ggarrange(phylo_2, phylo_3, ncol = 2, common.legend = FALSE,
                  legend = "bottom"))
  dev.off()
  
  return(NULL)
  
}

# Acer / length(theta)
pbmclapply(X = 1:length(theta),  
           FUN = application, 
           phylo_comm = phylo_comm, 
           phylo_observed_plant = aggregation$Acer,
           phylo_modelled_plant = phylo_modelled$Acer,
           legend = "Occurrence probability", 
           pallete = viridis(1000),
           export = paste0(root_path, "/figures/video/occurrence/Acer"), 
           class = "Acer",
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
           mc.allow.recursive = TRUE)

# Quercus
pbmclapply(X = 1:length(theta), #length(theta),  
           FUN = application, 
           phylo_comm = phylo_comm, 
           phylo_observed_plant = aggregation$Quercus,
           phylo_modelled_plant = phylo_modelled$Quercus,
           legend = "Occurrence probability", 
           pallete = viridis(1000),
           export = paste0(root_path, "/figures/video/occurrence/Quercus"), 
           class = "Quercus",
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
           mc.allow.recursive = TRUE)

# Pinus
pbmclapply(X = 1:length(theta), #length(theta),  
           FUN = application, 
           phylo_comm = phylo_comm, 
           phylo_observed_plant = aggregation$Pinus,
           phylo_modelled_plant = phylo_modelled$Pinus,
           legend = "Occurrence probability", 
           pallete = viridis(1000),
           export = paste0(root_path, "/figures/video/occurrence/Pinus"), 
           class = "Pinus",
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
           mc.allow.recursive = TRUE)

