#'------------------------------------------------------------------------------
#' @title Codes to visualize the MDS_predicted_mean values
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(ggplot2)
library(ggpubr)
library(bit64)
library(ggdark)
library(ggpointdensity)
library(rcartocolor)
library(viridis)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
sensor <- "EMIT"
folder <- paste0(root_path, "/06_model/SW_beta_PLS/", sensor)

#'------------------------------------------------------------------------------
#' Read information

performance <- fread(paste0(folder, "/", sensor, "_performance.csv"))
predicted_training <- fread(paste0(folder, "/", sensor, "_predicted_training.csv"))
predicted_testing <- fread(paste0(folder, "/", sensor, "_predicted_testing.csv"))

#'------------------------------------------------------------------------------
#' Data reshaping

work <- "Training"
#work <- "Testing"

testing <- performance[type == work]
testing <- testing[, .(B0 = mean(B0),
                           B1 = mean(B1)),
                       by = c("MDS")]

# Predict slopes and intercepts
values <- seq(-0.6, 0.6, by = 0.001)

#-------------------------------------------------------------------------------
# Plot details

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.93, 0.90),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.90, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth=0))

th <- th_trans
#th <- th_black
line_col <- "white"

point_size <- 2
point_shape <- 20
point_alpha <- 0.9

#-------------------------------------------------------------------------------
# Plot

frame <- predicted_training
#frame <- predicted_testing

range <- frame[, .(min_obs = min(MDS_observed),
          min_pred = min(MDS_predicted_mean),
          max_obs = max(MDS_observed),
          max_pred = max(MDS_predicted_mean)),
      by = "MDS"]

a1 <- ggplot(frame[MDS == "MDS1"], 
             aes(x = MDS_predicted_mean, 
                 y = MDS_observed)) +
  geom_abline(intercept = 0, 
              slope = 1, 
              colour = line_col, 
              linetype = "dotted") +
  geom_pointdensity(#adjust = 1, 
                    size = point_size,
                    shape = point_shape,
                    alpha = point_alpha) +
  #geom_point(size = point_size,
  #           shape = 21,
  #           alpha = 0.1,
  #           fill = NA,
  #           colour = "white") +
  #scale_colour_gradient2(low = "#ffc5aeff",
  #                       high = col[1],
  #                       guide = "none") +
  scale_colour_carto_c(type = "sequential", 
                       palette = "Peach", 
                       direction = 1,
                       guide = "none") +
  geom_abline(intercept = as.numeric(testing[1,2]), 
              slope = as.numeric(testing[1,3]), 
              colour = line_col, 
              linetype = "solid") +
  coord_cartesian(xlim = c(-0.58, 0.58), ylim = c(-0.58, 0.58)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(" ") + ylab("Observed") +
  th +
  labs(subtitle = "Axis 1")

a2 <- ggplot(frame[MDS == "MDS2"], 
             aes(x = MDS_predicted_mean, y = MDS_observed)) +
  geom_abline(intercept = 0, 
              slope = 1, 
              colour = line_col, 
              linetype = "dotted") +
  geom_pointdensity(#adjust = 3, 
                    size = point_size,
                    shape = point_shape,
                    alpha = point_alpha) +
  #geom_point(size = point_size,
  #           shape = 21,
  #           alpha = 0.1,
  #           fill = NA,
  #           colour = "white") +
  scale_colour_carto_c(type = "sequential", 
                       palette = "Mint", 
                       direction = 1,
                       guide = "none") +
  geom_abline(intercept = as.numeric(testing[1,2]), 
              slope = as.numeric(testing[1,3]), 
              colour = line_col, 
              linetype = "solid") +
  coord_cartesian(xlim = c(-0.58, 0.58), ylim = c(-0.58, 0.58)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted") + ylab(" ") +
  th +
  labs(subtitle = "Axis 2")

a3 <- ggplot(frame[MDS == "MDS3"], 
             aes(x = MDS_predicted_mean, y = MDS_observed)) +
  geom_abline(intercept = 0, 
              slope = 1, 
              colour = line_col, 
              linetype = "dotted") +
  geom_pointdensity(adjust = 2, 
                    size = point_size,
                    shape = point_shape,
                    alpha = point_alpha) +
  scale_colour_carto_c(type = "sequential", 
                       palette = "Teal", 
                       direction = 1,
                       guide = "none") +
  #geom_point(size = point_size,
  #           shape = 21,
  #           alpha = 0.1,
  #           fill = NA,
  #           colour = "white") +
  geom_abline(intercept = as.numeric(testing[1,2]), 
              slope = as.numeric(testing[1,3]), 
              colour = line_col, 
              linetype = "solid") +
  coord_cartesian(xlim = c(-0.48, 0.48), ylim = c(-0.48, 0.48)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(" ") + ylab(" ") +
  th +
  labs(subtitle = "Axis 3")

# Merge plots
plot <- ggarrange(a1, a2, a3, ncol=3, nrow=1)

# Export
png(paste0(root_path, "/figures/", sensor, "_predicted_training.png"), width = 210, height = 70, units = "mm", res = 300, bg = "transparent")

plot

dev.off()

