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
forest_type <- fread(paste0(root_path, "/FIA/00_FIA/ENTIRE_COND.csv"))
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
forest_type <- forest_type[forest_type$FORTYPCD %in% to_remove, FORTYP := NA]
forest_type <- forest_type[, c("PLT_CN", "FORTYP")]

# ------------------------------------------------------------------------------
# Taxo modelled

# Get glm coefficients
coefficients <- fread(paste0(root_path, "/Models/03_models_fortyp/", 
                             type, "_model-FORTYP.csv"))
coefficients <- subset(coefficients, dataset == "training")
coefficients <- coefficients[, .SD, .SDcols = c(1, 14:17)]

# Reshape coefficients
models <- unique(coefficients$forty)
coef_array <- array(dim = c(100, 3, length(models)))
intercept_mat <- matrix(nrow = 100, ncol = length(models))

for(i in 1:length(models)) {
  
  sub <- coefficients[forty == models[i]]
  coef_array[,,i] <- as.matrix(sub[, 3:5])
  intercept_mat[,i] <- sub$`(Intercept)`
  
}

fortyp_modelled <- summary_prob(as.matrix(taxo_comm[, c(-1, -5)]),
                                coef_array, 
                                intercept_mat)[[1]]
fortyp_modelled <- as.data.table(fortyp_modelled)/10000
colnames(fortyp_modelled) <- models

#'------------------------------------------------------------------------------
#' @Plot-FORTYP
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
                        taxo_comm, 
                        forest_type,
                        taxo_modelled_plant,
                        legend, 
                        pallete,
                        export, 
                        class,
                        exp_name,
                        theta, 
                        phi, 
                        h_adjust_axis1,
                        h_adjust_axis2,
                        th) {
  
  i <- X
  
  Occurrence <- forest_type
  Occurrence[FORTYP == class, Occurrence := "Presence"]
  Occurrence[FORTYP != class, Occurrence := "Absence"]
  Occurrence[is.na(FORTYP), Occurrence := "Absence"]
  Occurrence <- factor(Occurrence$Occurrence, levels = c("Presence", "Absence"))
  
  taxo_2 <- ggplot(taxo_comm, aes(x = Axis1, 
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
  
  taxo_3 <- ggplot(taxo_comm, aes(x = Axis1, 
                                    y = Axis2, 
                                    z = Axis3,
                                    color = taxo_modelled_plant)) +
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
  file_name = paste0(export, "/", exp_name, "_", name, ".png", sep="")
  png(file_name,
      width = 400, 
      height = 210, 
      units = "mm", 
      pointsize = 32,
      bg = NA,  
      res = 300)
  print(ggarrange(taxo_2, taxo_3, ncol = 2, common.legend = FALSE,
                  legend = "bottom"))
  dev.off()
  
  return(NULL)
  
}

# Red maple / lowland
pbmclapply(X = 1:length(theta),  
           FUN = application, 
           taxo_comm = taxo_comm, 
           forest_type = forest_type,
           taxo_modelled_plant = fortyp_modelled$`Red maple / lowland`,
           legend = "Occurrence probability", 
           pallete = viridis(1000),
           export = paste0(root_path, "/figures/video/FORTYP/Red-maple_lowland"), 
           class = "Red maple / lowland",
           exp_name = "Red-maple-lowland",
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

# White oak / red oak / hickory
pbmclapply(X = 1:length(theta),  
           FUN = application, 
           taxo_comm = taxo_comm, 
           forest_type = forest_type,
           taxo_modelled_plant = fortyp_modelled$`White oak / red oak / hickory`,
           legend = "Occurrence probability", 
           pallete = viridis(1000),
           export = paste0(root_path, "/figures/video/FORTYP/White-oak_red-oak_hickory"), 
           class = "White oak / red oak / hickory",
           exp_name = "White-oak_red-oak_hickory",
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

