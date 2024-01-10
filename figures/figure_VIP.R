#'------------------------------------------------------------------------------
#' @title Codes to visualize the VIP from the models
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
sensor <- "DESIS"
out_path <- paste0(root_path, "/06_model/SW_beta_PLS/", sensor)

#'------------------------------------------------------------------------------
#' Read information

vip <- coefficients <- fread(paste0(out_path, "/", sensor, "_vip.csv"))

#'------------------------------------------------------------------------------
#' Data reshaping

# Get bands names
bands <- colnames(vip)
bands <- bands[3:length(bands)]
bands <- as.numeric(bands)

# Get summary of vip
vip_melt <- melt(vip, id.vars = c("MDS", "iteration"), 
                 variable.name = "wavelength",
                 value.name = "VIP")


vip_summary <- vip_melt[, .(mean = mean(VIP), sd = sd(VIP)), 
                        by = c("MDS", "wavelength")]

vip_summary$wavelength <- as.numeric(as.character(vip_summary$wavelength))
vip_summary$MDS <- paste0("Axis ", vip_summary$MDS)

if(sensor == "EMIT") {
  
  add_na <- data.table(MDS = rep(c("Axis 1", "Axis 2", "Axis 3"), 2),
                       wavelength = c(1300, 1300, 1300, 1800, 1800, 1800),
                       mean = rep(c(NA, NA, NA), 2),
                       sd = rep(c(NA, NA, NA), 2))
  
  vip_summary <- rbind(vip_summary, add_na)
  vip_summary <- vip_summary[order(MDS, wavelength)]
  
}

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
                  legend.position = c(0.08, 0.90), 
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

plot <- ggplot(vip_summary) +
  geom_ribbon(aes(x = wavelength, ymin = mean - sd, ymax = mean + sd, fill = MDS), alpha = 0.25) +
  geom_line(aes(x = wavelength, y = mean, colour = MDS), linewidth = 0.5) +
  coord_cartesian(xlim = c(450, 900), ylim = c(0, 3.5), expand = FALSE) +
  #coord_cartesian(xlim = c(447, 2425), ylim = c(0, 3.5), expand = FALSE) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name = "", values = col) +
  xlab("Wavelength (nm)") + 
  ylab("VIP") +
  labs(subtitle = "Variable of Importance in Projection") + 
  th_trans

png(paste0(root_path, "/figures/", sensor, "_vip_transparent.png"), 
    width = 180, 
    height = 90, 
    pointsize = 12, 
    units = "mm", 
    res = 300, 
    bg = "transparent")

plot

dev.off()

