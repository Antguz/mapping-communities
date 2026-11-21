#'------------------------------------------------------------------------------
#' @title Codes to visualize the spectral and community dissimilarity (Fig 2)
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(bit64)
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

root_path <- "/media/guzman/Work/MMDB"

#-------------------------------------------------------------------------------
#' @Plot-details

th_black <- theme_bw(base_size = 12) + 
  theme(legend.position = "bottom",
        legend.title.align = 0.5,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        plot.margin = unit(c(18, 10, 1, 1), "pt"),
        rect = element_rect(fill = "transparent"))

th_trans <- theme_bw(base_size = 11) +
  theme(text = element_text(colour = "white"),
        axis.text = element_text(colour = "grey70"),
        axis.ticks = element_line(colour = "grey80"),
        panel.border = element_rect(colour = "grey80", fill = NA),
        panel.grid = element_line(colour = "grey8"),
        panel.background = element_rect(colour = NA, fill = "black"),
        plot.background = element_rect(colour = "black", fill = "black"),
        legend.background = element_rect(colour = NA, fill = "black"),
        legend.key = element_rect(colour = NA, fill = "black"),
        strip.background = element_rect(colour = "grey80", fill = "grey15"),
        strip.text = element_text(colour = "grey90")) +
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
#' @Function

mantel_plot <- function(mantel_results, 
                        msr_results,
                        msr_type,
                        exp = expression(paste("T", beta[total], ")")),
                        y_lab = "Spectral dissimilarity",
                        palette = "viridis") {
  
  # Get slopes
  dissimilarities <- mantel_results$frame
  
  # MSR Mantel results
  msr <- as.data.frame(msr_results)
  msr <- msr[msr$beta_type == msr_type, ]
  
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
  # sma <- mclapply(X = 1:length(unique_plots), 
  #                 FUN = sma_application,
  #                 plot_id = unique_plots,
  #                 frame = dissimilarities,
  #                 mc.cores = 28)
  # sma <- do.call(rbind.data.frame, sma)
  
  # pre <- matrix(rep(seq(0.00, 1.00, by = 0.1), nrow(sma)), 
  #               nrow = nrow(sma), byrow = TRUE)
  # pre <- pre*sma$slope
  # pre_mean <- apply(pre, 2, mean)
  # pre_sd <- apply(pre, 2, sd)
  x <- seq(0.00, 1.00, by = 0.1)
  
  #-----------------------------------------------------------------------------
  # Count estimations
  
  # Stored SAM is angle/180; x2 gives angle/90 (0-1) as in src/sam.cpp
  vox <- voxels(cbind(dissimilarities[, 4], spectral = dissimilarities$spectral, z = 1),
                edge_length = c(0.001, 0.002, 0.001),
                threads = 16, 
                obj.voxels = FALSE)
  
  vox[X >= 1.0, X := 0.9995]
  
  print(max(vox$N))
  
  #-----------------------------------------------------------------------------
  # Plot
  
  p <- ggplot() +
    geom_raster(data = vox, 
                aes(x = X, y = Y, fill = N)) +
    # geom_ribbon(aes(x = x,
    #                 ymin = pre_mean - pre_sd,
    #                 ymax = pre_mean + pre_sd),
    #             colour = line,
    #             fill = NA,
    #             linetype = "dashed",
    #             alpha = 0.5) +
    # geom_abline(data = sma, 
    #             aes(intercept = mean(elevation), 
    #                 slope = mean(slope)),
    #             colour = line) +
    scale_x_continuous(limits = c(0, 1.0), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 0.32), breaks = c(0, 0.16, 0.32), expand = c(0, 0)) +
    scale_fill_viridis_c(name = "Count",
                         option = palette, #palette
                         trans = "log10",
                         limits = c(1, 10^6.5), 
                         breaks = c(1, 10^2, 10^4, 10^6),
                         begin = 0.0,
                         end = 1.0,
                         labels = trans_format("log10", math_format(10^.x)),
                         oob = scales::squish) +
    xlab(exp) +
    ylab(y_lab) + 
    guides(fill = guide_colourbar(barwidth = 15, barheight = 0.7, title.position = "top")) +
    geom_text(aes(label = paste("Mantel test = ", round(msr$mantel, 2),
                                " (p = ", round(msr$mantel_significance, 3), ")", sep = ""),
                  x = 0.03, #0.03
                  y = 0.304, #0.26
                  hjust = 0,
                  vjust = 1),
              colour = text,
              size = 4) +
    th
  # p
   
  return(p)
  
}

#-------------------------------------------------------------------------------
#' @figures

# pal <- "cividis"
pal <- "inferno"
val <- c(2663174, 48882, 25456, 1638026, 27080, 13526, 2404258, 40086, 19638)

# DESIS ------------------------------------------------
sensor <- "DESIS"

type <- "taxo"
diversity <- "total"
DESIS_taxo_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_taxo_total_msr <- readRDS(paste0(root_path, "/Diversity/03_Mantel_MSR/", sensor, "_", type, "_mantel_msr.rds"))
DESIS_taxo_total_fig <- mantel_plot(mantel_results = DESIS_taxo_total,
                                    msr_results = DESIS_taxo_total_msr,
                                    msr_type = diversity,
                                    exp = expression(paste("T", beta[total])),
                                    palette = pal) #2663174
# rm(DESIS_taxo_total)
gc()

type <- "phylo"
diversity <- "total"
DESIS_phylo_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_phylo_total_msr <- readRDS(paste0(root_path, "/Diversity/03_Mantel_MSR/", sensor, "_", type, "_mantel_msr.rds"))
DESIS_phylo_total_fig <- mantel_plot(mantel_results = DESIS_phylo_total,
                                    msr_results = DESIS_phylo_total_msr,
                                    msr_type = diversity,
                                    exp = expression(paste("P", beta[total])),
                                    y_lab = " ",
                                    palette = pal) #48882
# rm(DESIS_phylo_total)
gc()

type <- "func"
diversity <- "total"
DESIS_func_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
DESIS_func_total_msr <- readRDS(paste0(root_path, "/Diversity/03_Mantel_MSR/", sensor, "_", type, "_mantel_msr.rds"))
DESIS_func_total_fig <- mantel_plot(mantel_results = DESIS_func_total,
                                     msr_results = DESIS_func_total_msr,
                                     msr_type = diversity,
                                     exp = expression(paste("F", beta[total])),
                                    y_lab = " ",
                                     palette = pal) #25456
# rm(DESIS_func_total)
gc()

# EMIT ------------------------------------------------
sensor <- "EMIT"

type <- "taxo"
diversity <- "total"
EMIT_taxo_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
EMIT_taxo_total_msr <- readRDS(paste0(root_path, "/Diversity/03_Mantel_MSR/", sensor, "_", type, "_mantel_msr.rds"))
EMIT_taxo_total_fig <- mantel_plot(mantel_results = EMIT_taxo_total,
                                    msr_results = EMIT_taxo_total_msr,
                                    msr_type = diversity,
                                    exp = expression(paste("T", beta[total])),
                                    palette = pal) #1017678
# rm(EMIT_taxo_total)
gc()

type <- "phylo"
diversity <- "total"
EMIT_phylo_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
EMIT_phylo_total_msr <- readRDS(paste0(root_path, "/Diversity/03_Mantel_MSR/", sensor, "_", type, "_mantel_msr.rds"))
EMIT_phylo_total_fig <- mantel_plot(mantel_results = EMIT_phylo_total,
                                     msr_results = EMIT_phylo_total_msr,
                                     msr_type = diversity,
                                     exp = expression(paste("P", beta[total])),
                                    y_lab = " ",
                                     palette = pal) #27080
# rm(EMIT_phylo_total)
gc()

type <- "func"
diversity <- "total"
EMIT_func_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
EMIT_func_total_msr <- readRDS(paste0(root_path, "/Diversity/03_Mantel_MSR/", sensor, "_", type, "_mantel_msr.rds"))
EMIT_func_total_fig <- mantel_plot(mantel_results = EMIT_func_total,
                                    msr_results = EMIT_func_total_msr,
                                    msr_type = diversity,
                                    exp = expression(paste("F", beta[total])),
                                   y_lab = " ",
                                    palette = pal) #13526
# rm(EMIT_func_total)
gc()

# EMITs ------------------------------------------------
sensor <- "EMITs"

type <- "taxo"
diversity <- "total"
EMITs_taxo_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
EMITs_taxo_total_msr <- readRDS(paste0(root_path, "/Diversity/03_Mantel_MSR/", sensor, "_", type, "_mantel_msr.rds"))
EMITs_taxo_total_fig <- mantel_plot(mantel_results = EMITs_taxo_total,
                                   msr_results = EMITs_taxo_total_msr,
                                   msr_type = diversity,
                                   exp = expression(paste("T", beta[total])),
                                   palette = pal) #2404258
# rm(EMITs_taxo_total)
gc()

type <- "phylo"
diversity <- "total"
EMITs_phylo_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
EMITs_phylo_total_msr <- readRDS(paste0(root_path, "/Diversity/03_Mantel_MSR/", sensor, "_", type, "_mantel_msr.rds"))
EMITs_phylo_total_fig <- mantel_plot(mantel_results = EMITs_phylo_total,
                                    msr_results = EMITs_phylo_total_msr,
                                    msr_type = diversity,
                                    exp = expression(paste("P", beta[total])),
                                    y_lab = " ",
                                    palette = pal) #40086
# rm(EMITs_phylo_total)
gc()

type <- "func"
diversity <- "total"
EMITs_func_total <- readRDS(paste0(root_path, "/Diversity/03_Mantel/", sensor, "_", type, "_", diversity, "_mantel.rds"))
EMITs_func_total_msr <- readRDS(paste0(root_path, "/Diversity/03_Mantel_MSR/", sensor, "_", type, "_mantel_msr.rds"))
EMITs_func_total_fig <- mantel_plot(mantel_results = EMITs_func_total,
                                   msr_results = EMITs_func_total_msr,
                                   msr_type = diversity,
                                   exp = expression(paste("F", beta[total])),
                                   y_lab = " ",
                                   palette = pal) #19638
# rm(EMITs_func_total)
gc()


# Merge plots
# Common legend
fig_legend <- get_legend(DESIS_taxo_total_fig + theme(legend.position = "top"))

# Sensor and spectral range (right side)
row_label <- function(sensor, from, to) {
  as_ggplot(grid::textGrob(bquote(atop(bold(.(sensor)), scriptstyle((lambda[.(from)] ~ "—" ~ lambda[.(to)])))),
                           rot = -90,
                           gp = grid::gpar(col = text, fontsize = 14)))
}

# Column headers
headers <- ggarrange(text_grob("Taxonomic", face = "bold", size = 13, color = text),
                     text_grob("Phylogenetic", face = "bold", size = 13, color = text),
                     text_grob("Functional", face = "bold", size = 13, color = text),
                     NULL,
                     ncol = 4, nrow = 1,
                     widths = c(1, 1, 1, 0.2))

row1 <- ggarrange(DESIS_taxo_total_fig + xlab(NULL) + ylab(NULL),
                  DESIS_phylo_total_fig + xlab(NULL) + ylab(NULL),
                  DESIS_func_total_fig + xlab(NULL) + ylab(NULL),
                  row_label("DESIS", 447, 895),
                  ncol = 4, nrow = 1,
                  widths = c(1, 1, 1, 0.2),
                  labels = c("A", "B", "C", ""),
                  font.label =  list(size = 12, face = "plain", color = text),
                  label.x = 0.02,
                  label.y = 1.0,
                  vjust = 1.2,
                  legend = "none")

row3 <- ggarrange(EMITs_taxo_total_fig + xlab(NULL) + ylab(NULL),
                  EMITs_phylo_total_fig + xlab(NULL) + ylab(NULL),
                  EMITs_func_total_fig + xlab(NULL) + ylab(NULL),
                  row_label("EMIT", 447, 894),
                  ncol = 4, nrow = 1,
                  widths = c(1, 1, 1, 0.2),
                  labels = c("D", "E", "F", ""),
                  font.label =  list(size = 12, face = "plain", color = text),
                  label.x = 0.02,
                  label.y = 1.0,
                  vjust = 1.2,
                  legend = "none")

row2 <- ggarrange(EMIT_taxo_total_fig + ylab(NULL),
                  EMIT_phylo_total_fig + ylab(NULL),
                  EMIT_func_total_fig + ylab(NULL),
                  row_label("EMIT", 447, 2426),
                  ncol = 4, nrow = 1,
                  widths = c(1, 1, 1, 0.2),
                  labels = c("G", "H", "I", ""),
                  font.label =  list(size = 12, face = "plain", color = text),
                  label.x = 0.02,
                  label.y = 1.0,
                  vjust = 1.2,
                  legend = "none")

# Export
jpeg(paste0(root_path, "/figures/Manuscript_black/", "figure_mantel_msr.jpeg"), width = 230, height = 230, units = "mm", res = 600, bg = "white")

annotate_figure(ggarrange(as_ggplot(fig_legend), headers, row1, row3, row2,
                          ncol = 1, nrow = 5,
                          heights = c(0.22, 0.07, 1*0.8, 1*0.8, 1.1*0.8)),
                left = text_grob("Spectral dissimilarity", rot = 90, color = text, size = 13),
                bottom = text_grob(expression(paste("Dimensions of ", beta, "-diversity")), color = text, size = 14))

dev.off()

# Transparent version (DESIS and EMIT); row1 has no x title here
png(paste0(root_path, "/figures/Manuscript_white/", "figure_mantel_msr.png"), width = 230, height = 170, units = "mm", res = 600, bg = "transparent")

annotate_figure(ggarrange(as_ggplot(fig_legend), headers, row1, row2,
                          ncol = 1, nrow = 4,
                          heights = c(0.22, 0.07, 1, 1.1)),
                left = text_grob("Spectral dissimilarity", rot = 90, color = text, size = 13),
                bottom = text_grob(expression(paste("Dimensions of ", beta, "-diversity")), color = text, size = 14))

dev.off()
