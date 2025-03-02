#'------------------------------------------------------------------------------
#' @title Codes to visualize the coefficients from the models
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(bit64)
library(ggdark)
library(ggpubr)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/MMDB"

#'------------------------------------------------------------------------------
#' Read information

# DESIS
sensor <- "DESIS"
out_path <- paste0(root_path, "/Models/01_model_beta/", sensor, "/", sensor, "_")
DESIS_taxo_coefficients <- fread(paste0(out_path, "taxo", "_coefficients.csv"))
DESIS_taxo_coefficients$type <- "Taxonomic"
DESIS_phylo_coefficients <- fread(paste0(out_path, "phylo", "_coefficients.csv"))
DESIS_phylo_coefficients$type <- "Phylogenetic"
DESIS_func_coefficients <- fread(paste0(out_path, "func", "_coefficients.csv"))
DESIS_func_coefficients$type <- "Functional"

DESIS <- rbind(DESIS_taxo_coefficients, 
               DESIS_phylo_coefficients,
               DESIS_func_coefficients)
DESIS$sensor <- sensor

DESIS <- melt(DESIS, id.vars = c("sensor", "type", "iteration", "Axes"),
              variable.name = "wavelength", value.name = "Coefficient")

DESIS <- DESIS[wavelength != "(Intercept)",]
DESIS$wavelength <- as.numeric(as.character(DESIS$wavelength))
DESIS$Coefficient <- DESIS$Coefficient*10000

# EMIT
sensor <- "EMIT"
out_path <- paste0(root_path, "/Models/01_model_beta/", sensor, "/", sensor, "_")
EMIT_taxo_coefficients <- fread(paste0(out_path, "taxo", "_coefficients.csv"))
EMIT_taxo_coefficients$type <- "Taxonomic"
EMIT_phylo_coefficients <- fread(paste0(out_path, "phylo", "_coefficients.csv"))
EMIT_phylo_coefficients$type <- "Phylogenetic"
EMIT_func_coefficients <- fread(paste0(out_path, "func", "_coefficients.csv"))
EMIT_func_coefficients$type <- "Functional"

EMIT <- rbind(EMIT_taxo_coefficients, 
              EMIT_phylo_coefficients,
              EMIT_func_coefficients)
EMIT$sensor <- sensor

EMIT <- melt(EMIT, id.vars = c("sensor", "type", "iteration", "Axes"),
              variable.name = "wavelength", value.name = "Coefficient")

EMIT <- EMIT[wavelength != "(Intercept)",]
EMIT$wavelength <- as.numeric(as.character(EMIT$wavelength))

#EMIT$Axis1 <- EMIT$Axis1/10000
#EMIT$Axis2 <- EMIT$Axis2/10000
#EMIT$Axis3 <- EMIT$Axis3/10000

# Merge
coefficients <- rbind(DESIS, EMIT)

# Summary 
coef_summary <- coefficients[, .(Coefficients = mean(Coefficient), 
                                 Coefficients_sd = sd(Coefficient)),
                             by = c("wavelength", "type", "sensor", "Axes")]
coef_summary$wavelength <- as.numeric(coef_summary$wavelength)
add_na <- data.table(wavelength = c(1300, 1300, 1300, 1800, 1800, 1800,
                               1300, 1300, 1300, 1800, 1800, 1800,
                               1300, 1300, 1300, 1800, 1800, 1800),
                     type = c(rep("Taxonomic", 6),
                              rep("Phylogenetic", 6), 
                              rep("Functional", 6)),
                     sensor = rep("EMIT", 18),
                     Axes = c("Axis1", "Axis2", "Axis3",
                              "Axis1", "Axis2", "Axis3",
                              "Axis1", "Axis2", "Axis3"),
                     Coefficients = rep(NA, 18),
                     Coefficients_sd = rep(NA, 18))

coef_summary <- rbind(coef_summary, add_na)
coef_summary <- coef_summary[order(wavelength)]

coef_summary$Axes <- factor(coef_summary$Axes, labels = c("Axis 1",
                                                          "Axis 2",
                                                          "Axis 3"))

coef_summary$type <- as.factor(coef_summary$type)
coef_summary$type <- factor(coef_summary$type, levels = c("Taxonomic",
                                                          "Phylogenetic",
                                                          "Functional"))

#-------------------------------------------------------------------------------
# Plot details

col <- c("#fc8d62", "#66c2a5", "#8da0cb")


th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"),
        legend.position= c("top"), 
        strip.background = element_rect(color="black", 
                                        fill="black", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "white"))

th_trans <- dark_theme_bw(base_size = 11) +
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
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "pt"),
        strip.background = element_rect(color="white", 
                                        fill="white", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "black"))

th <- th_trans
line_col <- "white"

th <- th_black
line_col <- "black"

#-------------------------------------------------------------------------------
# Get plot

DESIS_plot <- ggplot(coef_summary[sensor == "DESIS",]) +
  geom_ribbon(aes(x = wavelength, 
                  ymin = Coefficients - Coefficients_sd, 
                  ymax = Coefficients + Coefficients_sd, 
                  fill = Axes), 
              alpha = 0.25) +
  geom_line(aes(x = wavelength, 
                y = Coefficients, 
                colour = Axes), 
            linewidth = 0.5) +
  coord_cartesian(xlim = c(440, 910), ylim = c(-4, 4), expand = FALSE) +
  scale_y_continuous(breaks = c(-3, 0, 3)) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name = "", values = col) +
  xlab(" ") +
  ylab(" ") +
  #labs(subtitle = "PLS coefficients") +
  facet_grid(type ~ sensor) +
  th

DESIS_plot

EMIT_plot <- ggplot(coef_summary[sensor == "EMIT",]) +
  geom_ribbon(aes(x = wavelength, 
                  ymin = Coefficients - Coefficients_sd, 
                  ymax = Coefficients + Coefficients_sd, 
                  fill = Axes), 
              alpha = 0.25) +
  geom_line(aes(x = wavelength, 
                y = Coefficients, 
                colour = Axes), 
            linewidth = 0.5) +
  coord_cartesian(xlim = c(440, 2428), ylim = c(-35, 35), expand = FALSE) +
  scale_y_continuous(breaks = c(-20, 0, 20)) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name = "", values = col) +
  xlab(" ") +
  ylab(" ") +
  #labs(subtitle = "PLS coefficients") +
  facet_grid(type ~ sensor) +
  th

EMIT_plot

# Export JPEG

jpeg(paste0(root_path, "/figures/Manuscript_black/PLSR_coefficients.jpeg"), 
     width = 220, 
     height = 110, 
     pointsize = 12, 
     units = "mm", 
     res = 300, 
     bg = "white")

annotate_figure(ggarrange(DESIS_plot, EMIT_plot, ncol= 2, nrow=1, widths = c(0.3, 0.7)),
                left = text_grob("PLS coefficients", face = "plain", color = line_col, rot = 90),
                bottom = text_grob("Wavelength (nm)", face = "plain", color = line_col))

dev.off()


# Export PNG

png(paste0(root_path, "/figures/Manuscript_white/PLSR_coefficients.png"), 
     width = 220, 
     height = 110, 
     pointsize = 12, 
     units = "mm", 
     res = 300, 
     bg = "transparent")

annotate_figure(ggarrange(DESIS_plot, EMIT_plot, ncol= 2, nrow=1, widths = c(0.3, 0.7)),
                left = text_grob("PLS coefficients", face = "plain", color = line_col, rot = 90),
                bottom = text_grob("Wavelength (nm)", face = "plain", color = line_col))

dev.off()



