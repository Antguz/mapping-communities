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
library(smatr)
library(scales)
library(ggpubr)
options(scipen = 100000)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/MMDB"

#-------------------------------------------------------------------------------
#' @Plot-details

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th_black <- theme_bw() + 
  theme(legend.position = "bottom",
        legend.title.align = 0.5,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #legend.position = c(0.90, 0.90), 
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

line <- "grey75"
text <- "white"

#'------------------------------------------------------------------------------
#' @Function

mantel_plot <- function(mantel_results, 
                        exp = expression(paste("Taxonomic dissimilarity (", beta[total], ")")),
                        y_lab = "Spectral dissimilarity") {
  
  # Get slopes
  dissimilarities <- mantel_results$frame
  
  #-----------------------------------------------------------------------------
  # Slope estimation
  sma_application <- function(X, plot_id, frame) {
    
    # Subset
    plot <- frame[PLOT_CN == plot_id[X], ]
    
    # Test
    test <- sma(spectral ~ 0 + Beta, data = plot)
    elevation <- test$coef[[1]][1,1]
    slope <- test$coef[[1]][2,1]
    
    # Results
    return(data.table(PLOT_CN = X,
                      elevation = elevation,
                      slope = slope))
    
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
  # Plot
  
  p <- ggplot() +
    geom_hex(data = dissimilarities,#[1:10000,],
             aes(x = Beta, 
                 y = spectral),
             bins = 100) +
    geom_abline(aes(intercept = 0, 
                    slope = 1),
                colour = line,
                linetype = "dotted") +
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
    scale_x_continuous(limits = c(-0.01, 1.01), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.01, 1.01), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_fill_viridis_c(name = "Count", 
                         option = "inferno",
                         trans = "log10", 
                         limits = c(1, 10000000),  
                         breaks = c(1, 3162.278, 10000000),
                         labels = trans_format("log10", math_format(10^.x))) +
    xlab(exp) + 
    ylab(y_lab) + 
    guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7, title.position = "top")) +
    geom_text(aes(label = paste0("Slope = ", round(mean(sma$slope), 2), " ± ", round(sd(sma$slope), 2), "\n",
                                "r = ", round(mantel_results$mantel$statistic, 2), " (p < ", round(mantel_results$mantel$signif, 3), ")"),
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

sensor <- "DESIS"
type <- "taxo"

diversity <- "total"
DESIS_taxo_total <- readRDS(paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_taxo_total_fig <- mantel_plot(mantel_results = DESIS_taxo_total,
                                    exp = expression(paste("T", beta[total]))) #expression(paste("T", beta[total], ")")))

diversity <- "replacement"
DESIS_taxo_replacement <- readRDS(paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_taxo_replacement_fig <- mantel_plot(mantel_results = DESIS_taxo_replacement,
                                    exp = expression(paste("T", beta[replacement]))) #expression(paste("T", beta[replacement], ")")))

diversity <- "richness"
DESIS_taxo_richness <- readRDS(paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_taxo_richness_fig <- mantel_plot(mantel_results = DESIS_taxo_richness,
                                       exp = expression(paste("T", beta[richness]))) #expression(paste("T", beta[richness], ")")))

sensor <- "DESIS"
type <- "phylo-log"

diversity <- "total"
DESIS_phylo <- readRDS(paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_phylo_fig <- mantel_plot(mantel_results = DESIS_phylo,
                               exp = " ",#expression(paste("Phylogenetic dissimilarity (", beta[total], ")")))
                               y_lab = " ")

sensor <- "DESIS"
type <- "func"

diversity <- "total"
DESIS_func <- readRDS(paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_func_fig <- mantel_plot(mantel_results = DESIS_func,
                              exp = " ",#expression(paste("Functional dissimilarity (", beta[total], ")")))
                              y_lab = " ")

sensor <- "EMIT"
type <- "taxo"
diversity <- "total"

EMIT_taxo <- readRDS(paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
EMIT_taxo_fig <- mantel_plot(mantel_results = EMIT_taxo,
                              exp = expression(paste("Taxonomic dissimilarity (", beta[total], ")")))

sensor <- "EMIT"
type <- "phylo-log"
diversity <- "total"

EMIT_phylo <- readRDS(paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
EMIT_phylo_fig <- mantel_plot(mantel_results = EMIT_phylo,
                              exp = expression(paste("Phylogenetic dissimilarity (", beta[total], ")")),
                              y_lab = " ")

sensor <- "EMIT"
type <- "func"
diversity <- "total"
EMIT_func <- readRDS(paste0(root_path, "/06_mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
EMIT_func_fig <- mantel_plot(mantel_results = EMIT_func,
                              exp = expression(paste("Functional dissimilarity (", beta[total], ")")),
                             y_lab = " ")

# Merge plots
row1 <- ggarrange(DESIS_taxo_fig, DESIS_phylo_fig, DESIS_func_fig, 
                  ncol=3, nrow=1,
                  labels = c("A", "B", "C"),
                  font.label =  list(size = 9, face = "plain", color = text),
                  label.x = 0.02,
                  label.y = 1.05,
                  common.legend = TRUE,
                  legend = "bottom")

row2 <- ggarrange(EMIT_taxo_fig, EMIT_phylo_fig, EMIT_func_fig, 
                  ncol=3, nrow=1,
                  labels = c("D", "E", "F"),
                  font.label =  list(size = 9, face = "plain", color = text),
                  label.x = 0.02,
                  label.y = 1.05,
                  common.legend = TRUE,
                  legend = "bottom")

# Export
jpeg(paste0(root_path, "/figures/", "figure_mantel.jpeg"), width = 210, height = 150, units = "mm", res = 600, bg = "white")

annotate_figure(ggarrange(row1, row2, ncol= 1, nrow=2, common.legend = TRUE),
                top = text_grob("Taxonomic                                  Phylogenetic                                  Functional", face = "bold", color = "black"),
                right = text_grob("DESIS                                  EMIT", face = "bold", color = text, rot = -90))

dev.off()

png(paste0(root_path, "/figures/", "figure_mantel_black.png"), width = 210, height = 150, units = "mm", res = 600, bg = "transparent")

annotate_figure(ggarrange(row1, row2, ncol= 1, nrow=2, common.legend = TRUE),
                top = text_grob("Taxonomic                                  Phylogenetic                                  Functional", face = "bold", color = text),
                right = text_grob("DESIS                                  EMIT", face = "bold", color = text, rot = -90))

dev.off()
