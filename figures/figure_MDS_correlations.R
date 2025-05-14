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
library(ggplot2)
library(ggforce)
library(cowplot)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/MMDB/Diversity/04_MDS-corr/"

#'------------------------------------------------------------------------------
#' @Function
#-------------------------------------------------------------------------------

plot_correlation <- function(frame, title = "Taxonomic", y_lab = "Community-weighted mean of plant traits") {
  
  tframe <- transform(
    frame,
    angle = Values * (pi / 4),
    a = 0.45,
    b = (1 - abs(Values)) * 0.45 + 0.05
  )
  
  tframe[Significance <= 0.05, pvalue := "<= 0.05"]
  tframe[Significance > 0.05, pvalue := "> 0.05"]
  tframe$Values <- round(tframe$Values, 2)
  tframe$label <- sprintf("%.2f", tframe$Values)
  
  plot <- ggplot(tframe, aes(x = Axes, y = name)) +
    # geom_tile(aes(fill = Values), color = "white") +
    geom_ellipse(aes(x0 = as.numeric(Axes), 
                     y0 = as.numeric(name),
                     a = a, 
                     b = b, 
                     angle = angle,
                     color = pvalue,
                     fill = Values,
                     linetype = pvalue),
                 size = 0.5) +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    scale_fill_distiller(
      palette = "RdBu",  # or "PuOr", "BrBG", etc.
      direction = -1,     # or -1 to reverse
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    geom_text(aes(label = label, color = pvalue), size = 4) +
    scale_linetype_manual("p-value", values = c("solid", "dotted")) +
    scale_color_manual("p-value", values = c("grey25", "grey")) +
    labs(title = title) +
    ylab(y_lab)
  
  return(plot)

}

#'------------------------------------------------------------------------------
#' @Read, reshaped, and test correlation
#-------------------------------------------------------------------------------

# Plant traits -----------------------------------------------------------------

# Read data
taxo <- fread(paste0(root_path, "/MDS_traits-cor_", "taxo", "_total.csv"))
phylo <- fread(paste0(root_path, "/MDS_traits-cor_", "phylo", "_total.csv"))
func <- fread(paste0(root_path, "/MDS_traits-cor_", "func", "_total.csv"))

taxo <- melt(taxo,
             id.vars = c("plant"),
             measure.vars = list(
               Values = c("Axis1", "Axis2", "Axis3"),
               Significance = c("sig_Axis1", "sig_Axis2", "sig_Axis3")
             ),
             variable.name = "Axes")

taxo$Axes <- paste0("Axes ", taxo$Axes)

taxo[plant == "HTq", name := "Tree height (max)"]
taxo[plant == "SLq", name := "Slenderness"]
taxo[plant == "WD", name := "Wood density"]
taxo[plant == "shade", name := "Shade tolerance index"]
taxo[plant == "drought", name := "Drought tolerance index"]
taxo[plant == "waterlogging", name := "Waterlogging tolerance index"]
taxo[plant == "LMA", name := "Leaf mass per area"]
taxo[plant == "N", name := "Nitrogen"]

taxo$name <- as.factor(taxo$name)
taxo$name <- factor(taxo$name, 
                    levels = rev(c("Tree height (max)", 
                                   "Slenderness", 
                                   "Shade tolerance index", 
                                   "Drought tolerance index",
                                   "Waterlogging tolerance index",
                                   "Leaf mass per area",
                                   "Nitrogen",
                                   "Wood density")))

taxo$Axes <- as.factor(taxo$Axes)
taxo$Axes <- factor(taxo$Axes, levels = c("Axes 1", "Axes 2", "Axes 3"))

phylo <- melt(phylo,
             id.vars = c("plant"),
             measure.vars = list(
               Values = c("Axis1", "Axis2", "Axis3"),
               Significance = c("sig_Axis1", "sig_Axis2", "sig_Axis3")
             ),
             variable.name = "Axes")

phylo$Axes <- paste0("Axes ", phylo$Axes)

phylo[plant == "HTq", name := "Tree height (max)"]
phylo[plant == "SLq", name := "Slenderness"]
phylo[plant == "WD", name := "Wood density"]
phylo[plant == "shade", name := "Shade tolerance index"]
phylo[plant == "drought", name := "Drought tolerance index"]
phylo[plant == "waterlogging", name := "Waterlogging tolerance index"]
phylo[plant == "LMA", name := "Leaf mass per area"]
phylo[plant == "N", name := "Nitrogen"]

phylo$name <- as.factor(phylo$name)
phylo$name <- factor(phylo$name, 
                    levels = rev(c("Tree height (max)", 
                                   "Slenderness", 
                                   "Shade tolerance index", 
                                   "Drought tolerance index",
                                   "Waterlogging tolerance index",
                                   "Leaf mass per area",
                                   "Nitrogen",
                                   "Wood density")))

phylo$Axes <- as.factor(phylo$Axes)
phylo$Axes <- factor(phylo$Axes, levels = c("Axes 1", "Axes 2", "Axes 3"))

func <- melt(func,
             id.vars = c("plant"),
             measure.vars = list(
               Values = c("Axis1", "Axis2", "Axis3"),
               Significance = c("sig_Axis1", "sig_Axis2", "sig_Axis3")
             ),
             variable.name = "Axes")

func$Axes <- paste0("Axes ", func$Axes)

func[plant == "HTq", name := "Tree height (max)"]
func[plant == "SLq", name := "Slenderness"]
func[plant == "WD", name := "Wood density"]
func[plant == "shade", name := "Shade tolerance index"]
func[plant == "drought", name := "Drought tolerance index"]
func[plant == "waterlogging", name := "Waterlogging tolerance index"]
func[plant == "LMA", name := "Leaf mass per area"]
func[plant == "N", name := "Nitrogen"]

func$name <- as.factor(func$name)
func$name <- factor(func$name, 
                    levels = rev(c("Tree height (max)", 
                                   "Slenderness", 
                                   "Shade tolerance index", 
                                   "Drought tolerance index",
                                   "Waterlogging tolerance index",
                                   "Leaf mass per area",
                                   "Nitrogen",
                                   "Wood density")))

func$Axes <- as.factor(func$Axes)
func$Axes <- factor(func$Axes, levels = c("Axes 1", "Axes 2", "Axes 3"))

# Plots
taxo_plot <- plot_correlation(taxo, title = "Taxonomic", y_lab = "Community-weighted mean of plant traits")
phylo_plot <- plot_correlation(phylo, title = "Phylogenetic", y_lab = "Community-weighted mean of plant traits")
func_plot <- plot_correlation(func, title = "Functional", y_lab = "Community-weighted mean of plant traits")

trait_plot <- plot_grid(taxo_plot, phylo_plot, func_plot, 
                        ncol = 3, align = "v", labels = c("A", "B", "C"),
                        label_x = 0.02,
                        label_y = 0.90)

# Export
jpeg(paste0(root_path, "MDS_traits_cor.jpeg"), width = 75, height = 90, units = "mm", res = 600, bg = "white")

trait_plot

dev.off()
