#'------------------------------------------------------------------------------
#' @title Correlation among ordinations (Fig S5)
#'------------------------------------------------------------------------------

#' @description It correlates and plots MDS axes among multiple dimensions of 
#' diversity.

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(ggcorrplot)
library(ggstatsplot)
library(RColorBrewer)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB"

#'------------------------------------------------------------------------------
#' @Read, reshaped, and test correlation
#-------------------------------------------------------------------------------

# Read data
taxo <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_scores_comm_", "taxo", "_total.csv"))
phylo <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_scores_comm_", "phylo", "_total.csv"))
func <- fread(paste0(root_path, "/Diversity/03_MDS-post/MDS_scores_comm_", "func", "_total.csv"))

# Test order
all.equal(taxo$PLT_CN, phylo$PLT_CN)
all.equal(taxo$PLT_CN, func$PLT_CN)
all.equal(phylo$PLT_CN, func$PLT_CN)

# Compute correlation
axes <- cbind(taxo[, c(-1)], phylo[, c(-1)], func[, c(-1)])
colnames(axes) <- c("T Axis 1", "T Axis 2", "T Axis 3",
                    "P Axis 1", "P Axis 2", "P Axis 3",
                    "F Axis 1", "F Axis 2", "F Axis 3")
axes <- axes[, 9:1]

#'------------------------------------------------------------------------------
#' @Plot and export
#-------------------------------------------------------------------------------

p <- ggstatsplot::ggcorrmat(
  data = axes,
  matrix.type = "upper",
  type = "parametric",
  colors = c("darkred", "white", "steelblue")
)

jpeg(paste0(root_path, "/figures/Manuscript_black/", "MDS_correlation.jpeg"),
    width = 210, 
    height = 150, 
    units = "mm", 
    pointsize = 12,
    bg = "white",  
    res = 600)

p

dev.off()
