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
#' Read information

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
performance <- fread(paste0(root_path, "/06_model/SW_beta_PLS/performance.csv"))
predicted_training <- fread(paste0(root_path, "/06_model/SW_beta_PLS/predicted_training.csv"))
predicted_testing <- fread(paste0(root_path, "/06_model/SW_beta_PLS/predicted_testing.csv"))

#------------------------------------------
# Get summary of performance

testing <- performance[type == "Testing"]

testing <- performance[, .(B0 = mean(B0),
                           B1 = mean(B1)),
                       by = c("MDS")]

#------------------------------------------
# Predict slopes and intercepts

values <- seq(-0.6, 0.6, by = 0.001)

#------------------------------------------
# Do plot

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th <- theme_bw() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill= "transparent"),
            panel.background = element_rect(fill = "transparent"),
            rect = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent"))

frame <- predicted_testing

a1 <- ggplot(frame[MDS == "MDS1"]) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_errorbar(aes(xmin = MDS_predicted_mean - MDS_predicted_sd, xmax = MDS_predicted_mean + MDS_predicted_sd, y = MDS_observed), colour = "grey") +
  geom_point(aes(x = MDS_predicted_mean, y = MDS_observed), size = 2, shape = 21, fill = col[1], colour = "white") +
  geom_abline(intercept = as.numeric(testing[1,2]), slope = as.numeric(testing[1,3]), colour = "black", linetype = "solid") +
  coord_cartesian(xlim = c(-1.2, 1.31), ylim = c(-1.2, 1.31)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted (Axis 1)") + ylab("Observed (Axis 1)") +
  th

a2 <- ggplot(frame[MDS == "MDS2"]) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_errorbar(aes(xmin = MDS_predicted_mean - MDS_predicted_sd, xmax = MDS_predicted_mean + MDS_predicted_sd, y = MDS_observed), colour = "grey") +
  geom_point(aes(x = MDS_predicted_mean, y = MDS_observed), size = 2, shape = 21, fill = col[2], colour = "white") +
  geom_abline(intercept = as.numeric(testing[2,2]), slope = as.numeric(testing[2,3]), colour = "black", linetype = "solid") +
  coord_cartesian(xlim = c(-1.1, 0.8), ylim = c(-1.1, 0.8)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted (Axis 2)") + ylab("Observed (Axis 2)") +
  th

a3 <- ggplot(frame[MDS == "MDS3"]) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dotted") +
  geom_errorbar(aes(xmin = MDS_predicted_mean - MDS_predicted_sd, xmax = MDS_predicted_mean + MDS_predicted_sd, y = MDS_observed), colour = "grey") +
  geom_point(aes(x = MDS_predicted_mean, y = MDS_observed), size = 2, shape = 21, fill = col[3], colour = "white") +
  geom_abline(intercept = as.numeric(testing[3,2]), slope = as.numeric(testing[3,3]), colour = "black", linetype = "solid") +
  coord_cartesian(xlim = c(-1.17, 1.32), ylim = c(-1.17, 1.32)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted (Axis 3)") + ylab("Observed (Axis 3)") +
  th

p <- ggarrange(a1, a2, a3, ncol=3, nrow=1)

jpeg("figures/MDS_observed-MDS_predicted_mean.jpg", width = 210, height = 210, units = "mm", res = 300)
p
dev.off()

png("figures/MDS_observed-MDS_predicted_mean.png", width = 210, height = 80, units = "mm", res = 300, bg = "transparent")
p
dev.off()
  

new_ordination <- dcast(frame, sample ~ MDS, value.var = "MDS_predicted_mean")
old_ordination <- dcast(frame, sample ~ MDS, value.var = "MDS_observed")
sd_ordination <- dcast(frame, sample ~ MDS, value.var = "MDS_predicted_sd")
sd_ordination <- rowMeans(sd_ordination[, 2:4])
plot(new_ordination$MDS2, new_ordination$MDS3)
plot(ordination_testing$MDS2, ordination_testing$MDS3)

plot3d(new_ordination[, 2:4])
plot3d(old_ordination[, 2:4])

nrow(new_ordination)

for(i in 1:nrow(new_ordination)) {
  arrow3d(old_ordination[i, 2:4], new_ordination[i, 2:4], 
          type = "rotation", 
          col = "grey", 
          width = 0.25,
          thickness = 0.2,
          n = 10,
          s = 0.25,
          theta = pi/128)
}

# Black background -------------------------------------------------------------
frame <- predicted_testing

a1 <- ggplot(frame[MDS == "MDS1"]) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dotted") +
  geom_errorbar(aes(xmin = MDS_predicted_mean - MDS_predicted_sd, xmax = MDS_predicted_mean + MDS_predicted_sd, y = MDS_observed), colour = "grey30") +
  geom_point(aes(x = MDS_predicted_mean, y = MDS_observed), size = 0.5, shape = 19, colour = col[1]) +
  geom_abline(intercept = as.numeric(testing[1,2]), slope = as.numeric(testing[1,3]), colour = "white", linetype = "solid", linewidth = 0.5) +
  #coord_cartesian(xlim = c(-0.50, 0.50), ylim = c(-0.50, 0.50)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted_(Axis 1)") + ylab("MDS_observed (Axis 1)") +
  dark_theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.90, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white", size = 12),
        legend.text = element_text(color = "white"),
        axis.text=element_text(color="white"), 
        plot.background = element_rect(colour = "transparent", 
                                       fill = "transparent", 
                                       linewidth = 0))

a2 <- ggplot(frame[MDS == "MDS2"]) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dotted") +
  geom_errorbar(aes(xmin = MDS_predicted_mean - MDS_predicted_sd, xmax = MDS_predicted_mean + MDS_predicted_sd, y = MDS_observed), colour = "grey30") +
  geom_point(aes(x = MDS_predicted_mean, y = MDS_observed), size = 0.5, shape = 19, colour = col[2]) +
  geom_abline(intercept = as.numeric(testing[1,2]), slope = as.numeric(testing[1,3]), colour = "white", linetype = "solid", linewidth = 0.5) +
  #coord_cartesian(xlim = c(-0.41, 0.32), ylim = c(-0.41, 0.32)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted (Axis 2)") + ylab("MDS_observed (Axis 2)") +
  dark_theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.90, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white", size = 12),
        legend.text = element_text(color = "white"),
        axis.text=element_text(color="white"), 
        plot.background = element_rect(colour = "transparent", 
                                       fill = "transparent", 
                                       linewidth = 0))

a3 <- ggplot(frame[MDS == "MDS3"]) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dotted") +
  geom_errorbar(aes(xmin = MDS_predicted_mean - MDS_predicted_sd, xmax = MDS_predicted_mean + MDS_predicted_sd, y = MDS_observed), colour = "grey30") +
  geom_point(aes(x = MDS_predicted_mean, y = MDS_observed), size = 0.5, shape = 19, colour = col[3]) +
  geom_abline(intercept = as.numeric(testing[1,2]), slope = as.numeric(testing[1,3]), colour = "white", linetype = "solid", linewidth = 0.5) +
  #coord_cartesian(xlim = c(-0.52, 0.25), ylim = c(-0.52, 0.25)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("MDS_predicted_mean (Axis 3)") + ylab("MDS_observed (Axis 3)") +
  dark_theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.90, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white", size = 12),
        legend.text = element_text(color = "white"),
        axis.text=element_text(color="white"), 
        plot.background = element_rect(colour = "transparent", 
                                       fill = "transparent", 
                                       linewidth = 0))

p <- ggarrange(a1, a2, a3, ncol=3, nrow=1)


png("figures/MDS_observed-MDS_predicted_mean_black.png", width = 210, height = 70, units = "mm", res = 600, bg = "transparent")
p
dev.off()

size_point <- 0.2

a1 <- ggplot(frame[MDS == "MDS1"]) +
  geom_errorbar(aes(xmin = MDS_predicted_mean - MDS_predicted_sd, xmax = MDS_predicted_mean + MDS_predicted_sd, y = MDS_observed), colour = "grey30") +
  #geom_point(aes(x = MDS_predicted_mean, y = MDS_observed), size = 0.5, shape = 19, colour = col[1]) +
  geom_pointdensity(aes(x = MDS_predicted_mean, y = MDS_observed), size = size_point) +
  #scale_color_viridis(option = "G") +
  scale_color_carto_c(type = "quantitative", palette = "Mint", direction = -1) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dotted") +
  geom_abline(intercept = as.numeric(testing[1,2]), slope = as.numeric(testing[1,3]), colour = "white", linetype = "solid", linewidth = 0.5) +
  coord_cartesian(xlim = c(-1.1, 1.30), ylim = c(-1.1, 1.30)) +
  scale_x_continuous(expand = c(0, 0), n.breaks = 3) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 3) +
  xlab(" ") + ylab("Observed") +
  dark_theme_bw(base_size = 12) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white", size = 12),
        legend.text = element_text(color = "white"),
        axis.text=element_text(color="white"), 
        plot.background = element_rect(colour = "transparent", 
                                       fill = "transparent", 
                                       linewidth = 0)) +
  labs(subtitle = "Axis 1")

a2 <- ggplot(frame[MDS == "MDS2"]) +
  geom_errorbar(aes(xmin = MDS_predicted_mean - MDS_predicted_sd, xmax = MDS_predicted_mean + MDS_predicted_sd, y = MDS_observed), colour = "grey30") +
  #geom_point(aes(x = MDS_predicted_mean, y = MDS_observed), size = 0.5, shape = 19, colour = col[1]) +
  geom_pointdensity(aes(x = MDS_predicted_mean, y = MDS_observed), size = size_point) +
  #scale_color_viridis(option = "G") +
  scale_color_carto_c(type = "quantitative", palette = "Mint", direction = -1) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dotted") +
  geom_abline(intercept = as.numeric(testing[1,2]), slope = as.numeric(testing[1,3]), colour = "white", linetype = "solid", linewidth = 0.5) +
  coord_cartesian(xlim = c(-1.1, 0.93), ylim = c(-1.1, 0.93)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Predicted") + ylab(" ") +
  dark_theme_bw(base_size = 12) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white", size = 12),
        legend.text = element_text(color = "white"),
        axis.text=element_text(color="white"), 
        plot.background = element_rect(colour = "transparent", 
                                       fill = "transparent", 
                                       linewidth = 0)) +
  labs(subtitle = "Axis 2")

a3 <- ggplot(frame[MDS == "MDS3"]) +
  geom_errorbar(aes(xmin = MDS_predicted_mean - MDS_predicted_sd, xmax = MDS_predicted_mean + MDS_predicted_sd, y = MDS_observed), colour = "grey30") +
  #geom_point(aes(x = MDS_predicted_mean, y = MDS_observed), size = 0.5, shape = 19, colour = col[1]) +
  geom_pointdensity(aes(x = MDS_predicted_mean, y = MDS_observed), size = size_point) +
  #scale_color_viridis(option = "G") +
  scale_color_carto_c(type = "quantitative", palette = "Mint", direction = -1) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dotted") +
  geom_abline(intercept = as.numeric(testing[1,2]), slope = as.numeric(testing[1,3]), colour = "white", linetype = "solid", linewidth = 0.5) +
  coord_cartesian(xlim = c(-1.11, 1.18), ylim = c(-1.11, 1.18)) +
  scale_x_continuous(expand = c(0, 0), n.breaks = 3) +
  scale_y_continuous(expand = c(0, 0), n.breaks = 3) +
  xlab(" ") + ylab(" ") +
  dark_theme_bw(base_size = 12) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white", size = 12),
        legend.text = element_text(color = "white"),
        axis.text=element_text(color="white"), 
        plot.background = element_rect(colour = "transparent", 
                                       fill = "transparent", 
                                       linewidth = 0)) +
  labs(subtitle = "Axis 3")


p <- ggarrange(a1, a2, a3, ncol=3, nrow=1)


png("figures/observed-predicted_black.png", width = 210, height = 70, units = "mm", res = 600, bg = "transparent")
p
dev.off()

