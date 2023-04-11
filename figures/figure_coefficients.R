#'------------------------------------------------------------------------------
#' @title Codes to visualize the coefficients from the models
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
coefficients <- fread(paste0(out_path, "/coefficients.csv"))

# Get bands names
bands <- colnames(coefficients)
bands <- bands[4:length(bands)]
bands <- as.numeric(bands)

# Get summary of coefficients
coef_melt <- melt(coefficients, id.vars = c("MDS", "iteration"), 
                  variable.name = "wavelength",
                  value.name = "coefficients")

coef_melt <- coef_melt[wavelength != "intercept"]

coef_summary <- coef_melt[, .(mean = mean(coefficients), sd = sd(coefficients)), 
                        by = c("MDS", "wavelength")]

coef_summary$wavelength <- as.numeric(as.character(coef_summary$wavelength))
coef_summary$MDS <- paste0("MDS", coef_summary$MDS)

#-------------------------------------------------------------------------------
# Get plot

col <- c("#fc8d62", "#66c2a5", "#8da0cb")

th <- theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot <- ggplot(coef_summary) +
  geom_ribbon(aes(x = wavelength, ymin = mean - sd, ymax = mean + sd, fill = MDS), alpha = 0.25) +
  geom_line(aes(x = wavelength, y = mean, colour = MDS), size = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(460, 895)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.00027, 0.00027)) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name= "", values = col) +
  xlab("Wavelength (nm)") + ylab("PLS coefficients") +
  th + theme(legend.position = c(0.93, 0.90), 
             legend.background = element_rect(fill= 'transparent'))


jpeg("figures/pls_coefficients.jpg", width = 180, height = 90, units = "mm", res = 300)
plot
dev.off()
