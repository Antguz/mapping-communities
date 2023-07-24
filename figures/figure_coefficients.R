#'------------------------------------------------------------------------------
#' @title Codes to visualize the coefficients from the models
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(bit64)
library(ggdark)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Read information

# Folder root
root_path <- "/media/antonio/Work/Projects/DESIS_mapping/06_model/SW_beta_PLS"

# Get VIP
coefficients <- fread(paste0(root_path, "/coefficients.csv"))

# Get bands names
bands <- colnames(coefficients)
bands <- bands[4:length(bands)]
bands <- as.numeric(bands)

# Get summary of coefficients
coef_melt <- melt(coefficients, id.vars = c("MDS", "iteration"), 
                  variable.name = "wavelength",
                  value.name = "coefficients")

coef_melt <- coef_melt[wavelength != "Intercept"]
coef_melt$coefficients <- coef_melt$coefficients*10000

coef_summary <- coef_melt[, .(mean = mean(coefficients), sd = sd(coefficients)), 
                        by = c("MDS", "wavelength")]

coef_summary$wavelength <- as.numeric(as.character(coef_summary$wavelength))
coef_summary$MDS <- paste0("Axis ", coef_summary$MDS)

coef_summary[MDS == "Axis 1", mean := (mean-min(mean))/(max(mean)-min(mean)) * 2 - 1]
coef_summary[MDS == "Axis 2", mean := (mean-min(mean))/(max(mean)-min(mean)) * 2 - 1]
coef_summary[MDS == "Axis 3", mean := (mean-min(mean))/(max(mean)-min(mean)) * 2 - 1]
coef_summary[MDS == "Axis 1", sd := (sd-min(sd))/(max(sd)-min(sd)) * 2 - 1]
coef_summary[MDS == "Axis 2", sd := (sd-min(sd))/(max(sd)-min(sd)) * 2 - 1]
coef_summary[MDS == "Axis 3", sd := (sd-min(sd))/(max(sd)-min(sd)) * 2 - 1]

#-------------------------------------------------------------------------------
# Get plot

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th <- theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_coeffic <- ggplot(coef_summary) +
  geom_ribbon(aes(x = wavelength, ymin = mean - sd, ymax = mean + sd, fill = MDS), alpha = 0.25) +
  geom_line(aes(x = wavelength, y = mean, colour = MDS), linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(460, 895)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-3.7, 3.7)) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name= "", values = col) +
  xlab("Wavelength (nm)") +
  ylab(expression(paste("PLS coefficients (x10"^-6, ")", sep = ""))) + 
  th + theme(legend.position = c(0.93, 0.90), 
             legend.background = element_rect(fill= "transparent"),
             panel.background = element_rect(fill = "transparent"),
             rect = element_rect(fill = "transparent"))


jpeg("figures/pls_coefficients.jpg", width = 180, height = 90, units = "mm", res = 300)
plot
dev.off()

png("figures/pls_coefficients.png", width = 180, height = 90, units = "mm", res = 300, bg = "transparent")
plot_coeffic
dev.off()



# Black background
plot_coeffic <- ggplot(coef_summary) +
  geom_ribbon(aes(x = wavelength, ymin = mean - sd, ymax = mean + sd, fill = MDS), alpha = 0.1) +
  geom_hline(yintercept = c(0), colour = "grey25", linetype = "solid", linewidth = 0.5) +
  geom_line(aes(x = wavelength, y = mean, colour = MDS), linewidth = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(467, 918)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name= "", values = col) +
  xlab("Wavelength (nm)") +
  ylab(expression(paste("PLS coefficients (Normalized)", sep = ""))) + 
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
                                       size=0)) +
  labs(subtitle = "PLS coefficients")

png("figures/pls_coefficients_black.png", width = 140, height = 90, units = "mm", res = 600, bg = "transparent", pointsize = 12)
plot_coeffic
dev.off()


