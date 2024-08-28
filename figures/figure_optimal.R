#'------------------------------------------------------------------------------
#' @title Codes to visualize optimal number of components of the models
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(ggplot2)
library(data.table)
library(bit64)
library(ggdark)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' @Read files
#-------------------------------------------------------------------------------

#' Working directory

root_path <- "/media/antonio/Work/DESIS_mapping"

# DESIS
sensor <- "DESIS"
out_path <- paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor)
DESIS_taxo <- fread(paste0(out_path, "/", sensor, "_", "taxo_optimal.csv"))
DESIS_taxo$dimension <- "Taxonomic"
DESIS_phylo <- fread(paste0(out_path, "/", sensor, "_", "phylo-log_optimal.csv"))
DESIS_phylo$dimension <- "Phylogenetic"
DESIS_func <- fread(paste0(out_path, "/", sensor, "_", "func_optimal.csv"))
DESIS_func$dimension <- "Functional"
DESIS <- rbind(DESIS_taxo, DESIS_phylo, DESIS_func)
DESIS$sensor <- "DESIS"

# EMIT
sensor <- "EMIT"
out_path <- paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor)
EMIT_taxo <- fread(paste0(out_path, "/", sensor, "_", "taxo_optimal.csv"))
EMIT_taxo$dimension <- "Taxonomic"
EMIT_phylo <- fread(paste0(out_path, "/", sensor, "_", "phylo-log_optimal.csv"))
EMIT_phylo$dimension <- "Phylogenetic"
EMIT_func <- fread(paste0(out_path, "/", sensor, "_", "func_optimal.csv"))
EMIT_func$dimension <- "Functional"
EMIT <- rbind(EMIT_taxo, EMIT_phylo, EMIT_func)
EMIT$sensor <- "EMIT"

# Merge sensors
optimal <- rbind(EMIT, DESIS)

# Melt data in columns and clean
optimal <- melt(optimal, id.vars = c("sensor", "dimension", "Axes", "type", "estimate", "iteration"),
                variable.name = "comp", value.name = "value")
optimal <- optimal[comp != "(Intercept)",]
optimal$comp <- as.character(optimal$comp)
optimal[, comp := strsplit(comp, "comps")[[1]][1], by = seq_along(1:nrow(optimal))]
optimal$comp <- as.numeric(optimal$comp)

# Get summary
summary_optimal <- optimal[, .(mean = mean(value), sd = sd(value)),
                           by = c( "sensor", "dimension", "Axes", "type", "estimate", "comp")]

summary_optimal$dimension <- as.factor(summary_optimal$dimension)
summary_optimal$dimension <- factor(summary_optimal$dimension, levels = c("Taxonomic",
                                                                          "Phylogenetic",
                                                                          "Functional"))

# Enter the optimal number of comp
opt <- data.table(sensor = c(rep("DESIS", 9), rep("EMIT", 9)),
                  dimension = c(rep("Taxonomic", 3), rep("Phylogenetic", 3), rep("Functional", 3)),
                  Axes = rep(c("Axis1", "Axis2", "Axis3"), 6),
                  opt = c(23, 20, 17,
                          26, 12, 24,
                          25, 25, 25,
                          10, 15, 12,
                          15, 16, 10,
                          12, 14, 9))

opt$dimension <- as.factor(opt$dimension)
opt$dimension <- factor(opt$dimension, levels = c("Taxonomic",
                                                  "Phylogenetic",
                                                  "Functional"))

#-------------------------------------------------------------------------------
# Plot details

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.5, 2.0, 0.5, 0.5), "pt"),
        legend.position= c("top"), 
        strip.background = element_rect(color="black", 
                                        fill="black", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "white"))

th_trans <- dark_theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position= c("top"),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth=0),
        plot.margin = unit(c(0.5, 2.0, 0.5, 0.5), "pt"),
        strip.background = element_rect(color="white", 
                                        fill="white", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "black"))

#th <- th_trans
th <- th_black
line <- "black"
text <- "black"

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

#-------------------------------------------------------------------------------
#' @Plot 

frame <- summary_optimal[type == "RMSEP" & estimate == "train"]

plot_optimal <- ggplot(frame) +
  geom_ribbon(aes(x = comp, ymin = mean - sd, ymax = mean + sd, fill = Axes), alpha = 0.25) +
  geom_line(aes(x = comp, y = mean, colour = Axes), linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, 30)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0.11, 0.23), n.breaks = 4) +
  #geom_vline(xintercept = EMIT, colour = col, linetype = "dotted", linewidth = 1.0) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name= "", values = col) +
  xlab("Number of components") +
  ylab("RMSEP") +
  facet_grid(sensor ~ dimension) +
  geom_vline(data = opt,
             aes(xintercept = opt,
                 colour = Axes),
             linetype = "dotted",
             linewidth = 0.5) +
  th

plot_optimal


jpeg(paste0(root_path, "/figures/optimal-components.jpg"), 
     width = 210, 
     height = 150, 
     units = "mm", 
     res = 600, 
     #bg = "transparent",
     bg = "white",
     pointsize = 14)
plot_optimal
dev.off()

png(paste0(root_path, "/figures/sensors-model_CWM.png"), 
    width = 210, 
    height = 200, 
    units = "mm", 
    res = 600, 
    #bg = "transparent",
    bg = "white",
    pointsize = 14)
plot_optimal
dev.off()