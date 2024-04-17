#'------------------------------------------------------------------------------
#' @title Correlation among ordinations
#'------------------------------------------------------------------------------

#' @description It correlates and plots NMDS axes among multiple levels of 
#' diversity.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(ggcorrplot)
library(ggstatsplot)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#out_path <- "/media/antonio/Work/Projects/DESIS_mapping/07_model_occurrence"

#'------------------------------------------------------------------------------
#' @Read, reshaped, and test correlation
#-------------------------------------------------------------------------------

# Read data
taxo <- fread(paste0(root_path, "/05_NMDS/NMDS_scores_comm_", "taxo", ".csv"))
phylo <- fread(paste0(root_path, "/05_NMDS/NMDS_scores_comm_", "phylo-log", ".csv"))
func <- fread(paste0(root_path, "/05_NMDS/NMDS_scores_comm_", "func", ".csv"))

# Test order
all.equal(taxo$PLT_CN, phylo$PLT_CN)
all.equal(taxo$PLT_CN, func$PLT_CN)
all.equal(phylo$PLT_CN, func$PLT_CN)

# Compute correlation
axes <- cbind(taxo[, -1], phylo[, -1], func[, -1])
colnames(axes) <- c("T Axis 1", "T Axis 2", "T Axis 3",
                    "P Axis 1", "P Axis 2", "P Axis 3",
                    "F Axis 1", "F Axis 2", "F Axis 3")
axes <- axes[, 9:1]

#'------------------------------------------------------------------------------
#' @Plot and export
#-------------------------------------------------------------------------------

p <- ggstatsplot::ggcorrmat(
  data = axes,
  type = "bayes", # parametric for Pearson, nonparametric for Spearman's correlation
  colors = c("darkred", "white", "steelblue"), # change default colors
  cor.vars.names = c("Axis 1", "Axis 2", "Axis 3",
                     "Axis 1", "Axis 2", "Axis 3",
                     "Axis 1", "Axis 2", "Axis 3"),
  output = "ci",
  bf.prior = "ultrawide"
)

jpeg(paste0(root_path, "/figures/", "NMDS_correlation.jpeg"),
    width = 210, 
    height = 150, 
    units = "mm", 
    pointsize = 12,
    bg = "white",  
    res = 600)

p

dev.off()
