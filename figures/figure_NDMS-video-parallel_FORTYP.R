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

#'------------------------------------------------------------------------------
#' @Read-Information
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"
root_path <- "G:/Projects/MMDB"

#'------------------------------------------------------------------------------
#' Taxo ordination

type <- "taxo"

# Communities
taxo_comm <- fread(paste0(root_path, "/Diversity/02_NMDS/NMDS_scores_comm_", type, "_total.csv"))

# ------------------------------------------------------------------------------
# Taxo observed

# Meta information
fia_meta <- fread(paste0(root_path, "/FIA/00_FIA/plots_eastern_20230427.csv"))

# Test for match
fia_meta_match <- match.integer64(taxo_comm$PLT_CN, fia_meta$PLT_CN)
fia_meta <- fia_meta[fia_meta_match,]
all.equal(taxo_comm$PLT_CN, fia_meta$PLT_CN)

# FIA forest type
forest_type <- fread(paste0(root_path, "/FIA/00_FIA/ENTIRE_COND.csv"),  fill=TRUE)
forest_type <- forest_type[, c("PLT_CN", "FORTYPCD")]
FORTYPCD <- fread(paste0(root_path, "/FIA/00_FIA/ENTIRE_REF_FOREST_TYPE.csv"))
FORTYPCD <- FORTYPCD[, c("VALUE", "MEANING")]
forest_type <- merge(forest_type, FORTYPCD, 
                     by.x = "FORTYPCD", by.y = "VALUE",
                     all.x = TRUE, all.y = TRUE)
colnames(forest_type)[3] <- "FORTYP"

# Match plot order
FORTYP_match <- match.integer64(fia_meta$PLT_CN, forest_type$PLT_CN)
forest_type <- forest_type[FORTYP_match, ]
all.equal(fia_meta$PLT_CN, forest_type$PLT_CN)

# Clean forest type
to_remove <- c(181, 183, 223, 380, 382, 383, 803, 807, 925, 926, 
               932, 950, 951, 952, 953, 954, 955, 981, 999)
forest_type <- forest_type[!forest_type$FORTYPCD %in% to_remove, ]
forest_type <- forest_type[!is.na(FORTYP),]
forest_type <- forest_type[, c("PLT_CN", "FORTYP")]

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

phylo_modelled <- summary_cwm(as.matrix(phylo_comm[, c(-1, -5)]),
                              coef_array, 
                              intercept_mat)[[1]]
phylo_modelled <- as.data.table(phylo_modelled)/10000
colnames(phylo_modelled) <- plant

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
pbmclapply(X = 1:10,  
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
           #mc.silent = FALSE, 
           mc.cores = 1,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE)

# Quercus
pbmclapply(X = 1:2, #length(theta),  
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
           #mc.silent = FALSE, 
           mc.cores = 1,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE)

# Pinus
pbmclapply(X = 1:10, #length(theta),  
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
           #mc.silent = FALSE, 
           mc.cores = 1,
           mc.cleanup = TRUE, 
           mc.allow.recursive = TRUE)

