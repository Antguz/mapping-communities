#'------------------------------------------------------------------------------
#' @title Codes to visualize the stressplot from NMDS
#'------------------------------------------------------------------------------

#' @description Export a stressplot from NMDS. Our NMDS are huge, so this task 
#' most be completed in MSI.

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(bit64)
library(ggdark)
library(viridis)
library(scales)
library(ggpubr)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/MMDB"

#'------------------------------------------------------------------------------
#' Read information

taxo_to_export <- readRDS(paste0(root_path, "/Diversity/02_NMDS/data_stressplot_", "taxo", ".rds"))
phylo_to_export <- readRDS(paste0(root_path, "/Diversity/02_NMDS/data_stressplot_", "phylo", ".rds"))
func_to_export <- readRDS(paste0(root_path, "/Diversity/02_NMDS/data_stressplot_", "func", ".rds"))

#-------------------------------------------------------------------------------
# Plot visuals

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(10, 10, 1, 1), "pt"),
        legend.position= c("top"), 
        strip.background = element_rect(color="black", 
                                        fill="black", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "white"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(legend.position = "top",
        legend.title.align = 0.5,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = unit(c(10, 10, 1, 1), "pt"),
        #legend.position = "top", 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth = 0))

th <- th_black
col <- "grey20"
text <- "black"
line_col <- "black"

th <- th_trans
col <- "white"
text <- "white"
line_col <- "white"

pal <- "plasma"
pal <- "inferno"

#'------------------------------------------------------------------------------
#' Plot

max_values <- c(taxo_to_export$max, phylo_to_export$max, func_to_export$max)
max_value <- max_values[which.max(max_values)]

# Plot taxo
p_taxo <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = line_col) +
  geom_raster(data = taxo_to_export$vox, aes(x = X, y = Y, fill = N)) +
  geom_line(aes(x = taxo_to_export$line$x_line, 
                y = taxo_to_export$line$y_line), 
            colour = "grey75") +
  geom_abline(slope = 1, 
              intercept = 0, 
              linetype = "dotted", 
              linewidth = 0.5,
              colour = line_col) +
  xlab(expression(paste("Observed T", beta[total]))) +
  ylab("Ordination distance") +
  coord_cartesian(xlim = c(0, 1), 
                  ylim = c(0, 1), 
                  expand = FALSE) +
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_viridis_c(name = "Count", 
                       option = pal,
                       trans = "log10", 
                       limits = c(1, max_value),  
                       breaks = c(1, 1000, 1000000),
                       labels = trans_format("log10", math_format(10^.x))) +
  guides(fill = guide_colourbar(barwidth = 15, 
                                barheight = 0.7, 
                                title.position = "bottom")) +
  geom_text(aes(label = paste0("Non-metric fit, R2 = ", taxo_to_export$Rst, "\n",
                               " Linear fit, R2 = ", taxo_to_export$Ral),
                x = -Inf, 
                y = Inf, 
                hjust = -0.1, 
                vjust = 1.50),
            colour = col,
            size = 2,
            family = "Helvetica") +
  th +
  theme(panel.grid = element_blank())

# Plot phylo
p_phylo <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = line_col) +
  geom_raster(data = phylo_to_export$vox, aes(x = X, y = Y, fill = N)) +
  geom_line(aes(x = phylo_to_export$line$x_line, 
                y = phylo_to_export$line$y_line), 
            colour = "grey75") +
  geom_abline(slope = 1, 
              intercept = 0, 
              linetype = "dotted", 
              linewidth = 0.5,
              colour = line_col) +
  xlab(expression(paste("Observed P", beta[total]))) +
  ylab("Ordination distance") +
  coord_cartesian(xlim = c(0, 1), 
                  ylim = c(0, 1), 
                  expand = FALSE) +
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_viridis_c(name = "Count", 
                       option = pal,
                       trans = "log10", 
                       limits = c(1, max_value),  
                       breaks = c(1, 1000, 1000000),
                       labels = trans_format("log10", math_format(10^.x))) +
  guides(fill = guide_colourbar(barwidth = 15, 
                                barheight = 0.7, 
                                title.position = "bottom")) +
  geom_text(aes(label = paste0("Non-metric fit, R2 = ", phylo_to_export$Rst, "\n",
                               " Linear fit, R2 = ", phylo_to_export$Ral),
                x = -Inf,
                y = Inf,
                hjust = -0.1,
                vjust = 1.50),
            colour = col,
            size = 2,
            family = "Helvetica") +
  th +
  theme(panel.grid = element_blank())

# Func
p_func <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = line_col) +
  geom_raster(data = func_to_export$vox, aes(x = X, y = Y, fill = N)) +
  geom_line(aes(x = func_to_export$line$x_line, 
                y = func_to_export$line$y_line), 
            colour = "grey75") +
  geom_abline(slope = 1, 
              intercept = 0, 
              linetype = "dotted", 
              linewidth = 0.5,
              colour = line_col) +
  xlab(expression(paste("Observed F", beta[total]))) +
  ylab("Ordination distance") +
  coord_cartesian(xlim = c(0, 1), 
                  ylim = c(0, 1), 
                  expand = FALSE) +
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 3) +
  scale_fill_viridis_c(name = "Count", 
                       option = pal,
                       trans = "log10", 
                       limits = c(1, max_value),  
                       breaks = c(1, 1000, 1000000),
                       labels = trans_format("log10", math_format(10^.x))) +
  guides(fill = guide_colourbar(barwidth = 15, 
                                barheight = 0.7, 
                                title.position = "bottom")) +
  geom_text(aes(label = paste0("Non-metric fit, R2 = ", func_to_export$Rst, "\n",
                               " Linear fit, R2 = ", func_to_export$Ral),
                x = -Inf, 
                y = Inf, 
                hjust = -0.1, 
                vjust = 1.50),
            colour = col,
            size = 2,
            family = "Helvetica") +
  th +
  theme(panel.grid = element_blank())


# Export

panels <- ggarrange(p_taxo, p_phylo, p_func, 
                    ncol=3, nrow=1,
                    labels = c("A", "B", "C"),
                    font.label =  list(size = 9, face = "plain", color = col),
                    label.x = 0.00,
                    label.y = 1.00,
                    common.legend = TRUE,
                    legend = "bottom")


# Export
jpeg(paste0(root_path, "/figures/Manuscript_black/", "figure_stressplot.jpeg"), width = 210, height = 80, units = "mm", res = 600, bg = "white")

panels

dev.off()

png(paste0(root_path, "/figures/Manuscript_white/", "figure_stressplot.png"), width = 210, height = 80, units = "mm", res = 600, bg = "transparent")

panels

dev.off()
