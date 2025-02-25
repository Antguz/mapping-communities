#'------------------------------------------------------------------------------
#' @title Codes to visualize the spectral and community dissimilarity
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(bit64)
library(ggdark)
library(parallel)
library(viridis)
library(smatr)
library(scales)
library(ggpubr)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/MMDB"
root_path <- "E:/MMDB"

#-------------------------------------------------------------------------------
#' @Plot-details

th_black <- theme_bw() + 
  theme(legend.position = "bottom",
        legend.title.align = 0.5,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        plot.margin = unit(c(10, 10, 1, 1), "pt"),
        rect = element_rect(fill = "transparent"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.margin = unit(c(10, 10, 1, 1), "pt"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth = 0))

th <- th_trans
line <- "white"
text <- "white"

th <- th_black
line <- "grey25"
text <- "black"

#'------------------------------------------------------------------------------
#' @Read-files

taxo <- readRDS(paste0(root_path, "/Diversity/03_MDS-post/data_goodness_of_fit_taxo.rds"))
phylo <- readRDS(paste0(root_path, "/Diversity/03_MDS-post/data_goodness_of_fit_phylo.rds"))
func <- readRDS(paste0(root_path, "/Diversity/03_MDS-post/data_goodness_of_fit_func.rds"))

taxo_st <- fread(paste0(root_path, "/Diversity/03_MDS-post/goodness_of_fit_taxo_total.csv"))
phylo_st <- fread(paste0(root_path, "/Diversity/03_MDS-post/goodness_of_fit_phylo_total.csv"))
func_st <- fread(paste0(root_path, "/Diversity/03_MDS-post/goodness_of_fit_func_total.csv"))

#'------------------------------------------------------------------------------
#' @Plot

palette <- "viridis"

ptaxo <- ggplot() +
  geom_raster(data = taxo$vox, 
              aes(x = (X-5)/10000, y = (Y-5)/10000, fill = N)) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_abline(data = taxo$linear, 
              aes(intercept = mean(elevation-5)/10000, 
                  slope = mean(slope)),
              colour = line) +
  scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.119), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_fill_viridis_c(name = "Count",
                       option = palette,
                       trans = "log10",
                       limits = c(1, 10000),
                       breaks = c(1, 100, 10000),
                       begin = 0.0,
                       end = 1.0,
                       labels = trans_format("log10", math_format(10^.x)),
                       oob = scales::squish) +
  xlab(expression(paste("T", beta[total]))) +
  ylab("Ordination distances") + 
  guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7, title.position = "top")) +
  geom_text(mapping = aes(label = paste0("Intercept = ", round(mean((taxo$linear$elevation-5)/10000), 2), " ± ", round(sd((taxo$linear$elevation-5)/10000), 2), "\n",
                                         "Slope = ", round(mean(taxo$linear$slope), 2), " ± ", round(sd(taxo$linear$slope), 2), "\n",
                                         "Goodness-of-fit = ", round(taxo_st$rsquared, 2)),
                x = 0.03, #0.03
                y = 0.90, #0.26
                hjust = 0,
                vjust = 0),
            colour = text,
            size = 2.5) +
  th

ptaxo
  
pphylo <- ggplot() +
  geom_raster(data = phylo$vox, 
              aes(x = (X-5)/10000, y = (Y-5)/10000, fill = N)) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_abline(data = phylo$linear, 
              aes(intercept = mean(elevation-5)/10000, 
                  slope = mean(slope)),
              colour = line) +
  scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.119), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_fill_viridis_c(name = "Count",
                       option = palette,
                       trans = "log10",
                       limits = c(1, 10000),
                       breaks = c(1, 100, 10000),
                       begin = 0.0,
                       end = 1.0,
                       labels = trans_format("log10", math_format(10^.x)),
                       oob = scales::squish) +
  xlab(expression(paste("P", beta[total]))) +
  ylab(" ") + 
  guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7, title.position = "top")) +
  geom_text(mapping = aes(label = paste0("Intercept = ", round(mean((phylo$linear$elevation-5)/10000), 2), " ± ", round(sd((phylo$linear$elevation-5)/10000), 2), "\n",
                                         "Slope = ", round(mean(phylo$linear$slope), 2), " ± ", round(sd(phylo$linear$slope), 2), "\n",
                                         "Goodness-of-fit = ", round(phylo_st$rsquared, 2)),
                          x = 0.03, #0.03
                          y = 0.90, #0.26
                          hjust = 0,
                          vjust = 0),
            colour = text,
            size = 2.5) +
  th

pphylo

pfunc <- ggplot() +
  geom_raster(data = func$vox, 
              aes(x = (X-5)/10000, y = (Y-5)/10000, fill = N)) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_abline(data = func$linear, 
              aes(intercept = mean(elevation-5)/10000, 
                  slope = mean(slope)),
              colour = line) +
  scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.119), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_fill_viridis_c(name = "Count",
                       option = palette,
                       trans = "log10",
                       limits = c(1, 10000),
                       breaks = c(1, 100, 10000),
                       begin = 0.0,
                       end = 1.0,
                       labels = trans_format("log10", math_format(10^.x)),
                       oob = scales::squish) +
  xlab(expression(paste("F", beta[total]))) +
  ylab(" ") + 
  guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7, title.position = "top")) +
  geom_text(mapping = aes(label = paste0("Intercept = ", round(mean((func$linear$elevation-5)/10000), 2), " ± ", round(sd((func$linear$elevation-5)/10000), 2), "\n",
                                         "Slope = ", round(mean(func$linear$slope), 2), " ± ", round(sd(func$linear$slope), 2), "\n",
                                         "Goodness-of-fit = ", round(func_st$rsquared, 2)),
                          x = 0.03, #0.03
                          y = 0.90, #0.26
                          hjust = 0,
                          vjust = 0),
            colour = text,
            size = 2.5) +
  th

pfunc

# Merge plots
row1 <- ggarrange(ptaxo, pphylo, pfunc, 
                  ncol=3, nrow=1,
                  labels = c("A", "B", "C"),
                  font.label =  list(size = 9, face = "plain", color = text),
                  label.x = 0.02,
                  label.y = 1.00,
                  common.legend = TRUE,
                  legend = "bottom")

jpeg(paste0(root_path, "/figures/Manuscript_black/", "figure_goodness-of-fit.jpeg"), 
     width = 210, 
     height = 85, 
     units = "mm", 
     res = 600, 
     bg = "white")

row1

dev.off()
