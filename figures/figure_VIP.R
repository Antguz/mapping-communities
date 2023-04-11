#'------------------------------------------------------------------------------
#' @title Codes to visualize the VIP from the models
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)

options(scipen = 100)

#'------------------------------------------------------------------------------
#' Read information

#' Folder root
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
out_path <- paste0(root_path, "/06_model/SW_beta_PLS")

# Get VIP
vip <- fread(paste0(out_path, "/vip.csv"))

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
vip_summary$MDS <- paste0("MDS", vip_summary$MDS)

#-------------------------------------------------------------------------------
# Get plot


col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th <- theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot <- ggplot(vip_summary) +
  geom_ribbon(aes(x = wavelength, ymin = mean - sd, ymax = mean + sd, fill = MDS), alpha = 0.25) +
  geom_line(aes(x = wavelength, y = mean, colour = MDS), size = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(460, 895)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3.1)) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name = "", values = col) +
  xlab("Wavelength (nm)") + ylab("Variable Importance in Projection") +
  th + theme(legend.position = c(0.93, 0.90), 
             legend.background = element_rect(fill= 'transparent'))

jpeg("figures/VIP.jpg", width = 180, height = 90, units = "mm", res = 300)
plot
dev.off()
