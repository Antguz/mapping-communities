#'------------------------------------------------------------------------------
#' @title Codes to visualize optimal number of components of the models
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(bit64)
library(ggdark)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "EMIT"
out_path <- paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor)
type <- "taxo"

#'------------------------------------------------------------------------------
#' Read information

# Get RMSEP
optimal <- fread(paste0(out_path, "/", sensor, "_", type, "_optimal.csv"))

rmse <- optimal[type == "RMSEP" & estimate == "train"]
rmse <- rmse[, .SD, .SDcols = c(4, 6:ncol(rmse))]
colnames(rmse) <- c("Axes", 1:(ncol(rmse)-1))
rmse <- melt(rmse, id.vars = "Axes", variable.name = "comp", value.name = "RMSEP")
rmse$comp <- as.numeric(rmse$comp)

rmse_cv <- optimal[type == "RMSEP" & estimate == "CV"]
rmse_cv <- rmse_cv[, .SD, .SDcols = c(4, 6:ncol(rmse_cv))]
colnames(rmse_cv) <- c("Axes", 1:(ncol(rmse_cv)-1))
rmse_cv <- melt(rmse_cv, id.vars = "Axes", variable.name = "comp", value.name = "RMSEP")
rmse_cv$comp <- as.numeric(rmse_cv$comp)

# Summary by MDS
summary_rmsep <- rmse[, .(mean = mean(RMSEP), sd = sd(RMSEP)), 
                       by = c("Axes", "comp")]

summary_rmsep_cv <- rmse_cv[, .(mean = mean(RMSEP), sd = sd(RMSEP)), 
                      by = c("Axes", "comp")]

min_rmse <- optimal[type == "RMSEP" & estimate == "CV",]
opt_min_rmse <- min_rmse[, .SD, .SDcols = c(6:ncol(min_rmse))]
opt_min_rmse <- opt_min_rmse[, .(opt = which.min(.SD)), by = .I]
opt_min_rmse$Axes <- min_rmse$Axes

mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

opt <- opt_min_rmse[, .(mode_opt = mode(opt), 
                        min_opt = min(opt), 
                        max_opt = max(opt),
                        meadian_opt = median(opt),
                        mean_opt = mean(opt)), by = Axes] 

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Plot details

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #legend.position = c(0.93, 0.90),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #legend.position = c(0.08, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth = 0))

#-------------------------------------------------------------------------------
# Get plot

DESIS <- c(10.5, 3, 6)
EMIT <- c(10.5, 3, 6)

plot_optimal <- ggplot(summary_rmsep) +
  geom_ribbon(aes(x = comp, ymin = mean - sd, ymax = mean + sd, fill = Axes), alpha = 0.25) +
  geom_line(aes(x = comp, y = mean, colour = Axes), linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, max(summary_rmsep$comp))) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0.11, 0.23), n.breaks = 4) +
  geom_vline(xintercept = DESIS, colour = col, linetype = "dotted", linewidth = 1.0) +
  #geom_vline(xintercept = EMIT, colour = col, linetype = "dotted", linewidth = 1.0) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name= "", values = col) +
  xlab("Number of components") +
  ylab("RMSEP") +  
  #th_trans + 
  theme(legend.position = c(0.93, 0.90)) + theme_bw()

plot_optimal_cv <- ggplot(summary_rmsep_cv) +
  geom_ribbon(aes(x = comp, ymin = mean - sd, ymax = mean + sd, fill = Axes), alpha = 0.25) +
  geom_line(aes(x = comp, y = mean, colour = Axes), linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(1, max(summary_rmsep$comp))) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0.11, 0.23), n.breaks = 4) +
  geom_vline(xintercept = DESIS, colour = col, linetype = "dotted", linewidth = 1.0) +
  #geom_vline(xintercept = EMIT, colour = col, linetype = "dotted", linewidth = 1.0) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name= "", values = col) +
  xlab("Number of components") +
  ylab("RMSEP") +  
  #th_trans + 
  theme(legend.position = c(0.93, 0.90)) + theme_bw()

png(paste0(root_path, "/figures/", sensor, "_optimal_transparent.png"), 
    width = 150, 
    height = 90, 
    pointsize = 12, 
    units = "mm", 
    res = 300, 
    bg = "transparent")

plot_optimal

dev.off()

# Taxonomic 18
# Phylogenetic 15
# Phylogenetic-log 15
# Func 24