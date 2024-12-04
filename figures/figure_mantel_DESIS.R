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
library(ggpointdensity)
library(viridis)
library(rcartocolor)
library(smatr)
library(scales)
library(ggpubr)
library(rTLS)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/MMDB"

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
line <- "grey75"
text <- "black"

#'------------------------------------------------------------------------------
#' @Function

mantel_plot <- function(mantel_results, 
                        exp = expression(paste("T", beta[total], ")")),
                        y_lab = "Spectral dissimilarity",
                        palette = "viridis") {
  
  # Get slopes
  dissimilarities <- mantel_results$frame
  
  #-----------------------------------------------------------------------------
  # Slope estimation
  sma_application <- function(X, plot_id, frame) {
    
    # Subset
    plot <- frame[PLOT_CN == plot_id[X], ]
    
    # Test
    test <- sma(spectral ~ Beta, data = plot)
    #test <- sma(spectral ~ Beta, data = plot)
    elevation <- test$coef[[1]][1,1]
    slope <- test$coef[[1]][2,1]
    
    # Results
    return(data.table(PLOT_CN = X,
                      elevation = elevation,
                      slope = slope,
                      rsquar = test$r2[[1]]))
    
  }
  
  unique_plots <- unique(dissimilarities$PLOT_CN)
  sma <- mclapply(X = 1:length(unique_plots), 
                  FUN = sma_application,
                  plot_id = unique_plots,
                  frame = dissimilarities,
                  mc.cores = 28)
  sma <- do.call(rbind.data.frame, sma)
  
  pre <- matrix(rep(seq(0.00, 1.00, by = 0.1), nrow(sma)), 
                nrow = nrow(sma), byrow = TRUE)
  pre <- pre*sma$slope
  pre_mean <- apply(pre, 2, mean)
  pre_sd <- apply(pre, 2, sd)
  x <- seq(0.00, 1.00, by = 0.1)
  
  #-----------------------------------------------------------------------------
  # Count estimations
  
  vox <- voxels(cbind(dissimilarities[, 4:3], z = 1), 
                edge_length = c(0.001, 0.001, 0.001), 
                threads = 16, 
                obj.voxels = FALSE)
  
  vox[X >= 1.0, X := 0.9995]
  
  print(max(vox$N))
  
  #-----------------------------------------------------------------------------
  # Plot
  
  p <- ggplot() +
    geom_raster(data = vox, 
                aes(x = X, y = Y, fill = N)) +
    geom_ribbon(aes(x = x,
                    ymin = pre_mean - pre_sd,
                    ymax = pre_mean + pre_sd),
                colour = line,
                fill = NA,
                linetype = "dashed",
                alpha = 0.5) +
    geom_abline(data = sma, 
                aes(intercept = mean(elevation), 
                    slope = mean(slope)),
                colour = line) +
    scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.05, 0.1, 0.15), expand = c(0, 0)) +
    scale_fill_viridis_c(name = "Count",
                         option = palette,
                         trans = "log10",
                         limits = c(1, 1000000),
                         breaks = c(1, 1000, 1000000),
                         begin = 0.1,
                         labels = trans_format("log10", math_format(10^.x))) +
    xlab(exp) +
    ylab(y_lab) + 
    guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7, title.position = "top")) +
    geom_text(aes(label = paste0("Intercept = ", round(mean(sma$elevation), 2), " ± ", round(sd(sma$elevation), 2), "\n",
                                 " Slope = ", round(mean(sma$slope), 2), " ± ", round(sd(sma$slope), 2), "\n",
                                 "Mantel test = ", round(mantel_results$mantel, 2), " (p = ", round(mantel_results$mantel_significance, 3), ")"),
                  x = -Inf, 
                  y = Inf, 
                  hjust = -0.1, 
                  vjust = 1.50),
              colour = text,
              size = 2.5) +
    th
  #p
  
  return(p)
  
}

#-------------------------------------------------------------------------------
#' @figures

pal <- "inferno"
sensor <- "DESIS"

# total ------------------------------------------------

type <- "taxo"
diversity <- "total"
DESIS_taxo_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_taxo_total_fig <- mantel_plot(mantel_results = DESIS_taxo_total,
                                    exp = expression(paste("T", beta[total])),
                                    palette = pal) #1030468
# rm(DESIS_taxo_total)
gc()

type <- "phylo"
diversity <- "total"
DESIS_phylo_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_phylo_total_fig <- mantel_plot(mantel_results = DESIS_phylo_total,
                                     exp = expression(paste("P", beta[total])),
                                     palette = pal) #33640
# rm(DESIS_phylo_total)
gc()

type <- "func"
diversity <- "total"
DESIS_func_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_func_total_fig <- mantel_plot(mantel_results = DESIS_func_total,
                                    exp = expression(paste("F", beta[total])),
                                    palette = pal) #37090
# rm(DESIS_func_total)
gc()

# replacement ------------------------------------------------

type <- "taxo"
diversity <- "replacement"
DESIS_taxo_repl <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_taxo_repl_fig <- mantel_plot(mantel_results = DESIS_taxo_repl,
                                    exp = expression(paste("T", beta[repl/total])),
                                    palette = pal) #1030468
# rm(DESIS_taxo_total)
gc()

type <- "phylo"
diversity <- "replacement"
DESIS_phylo_repl <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_phylo_repl_fig <- mantel_plot(mantel_results = DESIS_phylo_repl,
                                     exp = expression(paste("P", beta[repl/total])),
                                     palette = pal) #33640
# rm(DESIS_phylo_total)
gc()

type <- "func"
diversity <- "replacement"
DESIS_func_repl <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_func_repl_fig <- mantel_plot(mantel_results = DESIS_func_repl,
                                    exp = expression(paste("F", beta[repl/total])),
                                    palette = pal) #37090
# rm(DESIS_func_total)
gc()

# richness ------------------------------------------------

type <- "taxo"
diversity <- "richness"
DESIS_taxo_rich <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_taxo_rich_fig <- mantel_plot(mantel_results = DESIS_taxo_rich,
                                   exp = expression(paste("T", beta[rich/total])),
                                   palette = pal) #1030468
# rm(DESIS_taxo_rich)
gc()

type <- "phylo"
diversity <- "richness"
DESIS_phylo_rich <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_phylo_rich_fig <- mantel_plot(mantel_results = DESIS_phylo_rich,
                                    exp = expression(paste("P", beta[rich/total])),
                                    palette = pal) #33640
# rm(DESIS_phylo_rich)
gc()

type <- "func"
diversity <- "richness"
DESIS_func_rich <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_func_rich_fig <- mantel_plot(mantel_results = DESIS_func_rich,
                                   exp = expression(paste("F", beta[rich/total])),
                                   palette = pal) #37090
# rm(DESIS_func_rich)
gc()

# Merge plots
row1 <- ggarrange(DESIS_taxo_total_fig, DESIS_phylo_total_fig, DESIS_func_total_fig, 
                  ncol=3, nrow=1,
                  labels = c("A", "B", "C"),
                  font.label =  list(size = 9, face = "plain", color = text),
                  label.x = 0.02,
                  label.y = 1.05,
                  common.legend = TRUE,
                  legend = "top")

row2 <- ggarrange(DESIS_taxo_repl_fig, DESIS_phylo_repl_fig, DESIS_func_repl_fig, 
                  ncol=3, nrow=1,
                  labels = c("D", "E", "F"),
                  font.label =  list(size = 9, face = "plain", color = text),
                  label.x = 0.02,
                  label.y = 1.05,
                  common.legend = TRUE,
                  legend = "top")

row3 <- ggarrange(DESIS_taxo_rich_fig, DESIS_phylo_rich_fig, DESIS_func_rich_fig, 
                  ncol=3, nrow=1,
                  labels = c("G", "H", "I"),
                  font.label =  list(size = 9, face = "plain", color = text),
                  label.x = 0.02,
                  label.y = 1.05,
                  common.legend = TRUE,
                  legend = "top")

# Export
jpeg(paste0(root_path, "/figures/Manuscript_black/", "figure_mantel_DESIS.jpeg"), width = 230, height = 240, units = "mm", res = 600, bg = "white")

annotate_figure(ggarrange(row1, row2, row3, ncol= 1, nrow=3, common.legend = TRUE),
                top = text_grob("Taxonomic                                  Phylogenetic                                  Functional", face = "bold", color = "black"),
                right = text_grob("Total Replacement/Total Richness/Total", face = "bold", color = text, rot = -90))

dev.off()

png(paste0(root_path, "/figures/Manuscript_white/", "figure_mantel_DESIS.png"), width = 230, height = 240, units = "mm", res = 600, bg = "transparent")

annotate_figure(ggarrange(row1, row2, row3, ncol= 1, nrow=3, common.legend = TRUE),
                top = text_grob("Taxonomic                                  Phylogenetic                                  Functional", face = "bold", color = text),
                right = text_grob("Total Replacement/Total Richness/Total", face = "bold", color = text, rot = -90))

dev.off()
