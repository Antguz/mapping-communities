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
#' Working directory

root_path <- "E:/Projects/DESIS_mapping"

#'------------------------------------------------------------------------------
#' Read information

# DESIS
sensor <- "DESIS"
out_path <- paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/", sensor, "_")
DESIS_taxo_coefficients <- fread(paste0(out_path, "taxo", "_coefficients.csv"))
DESIS_taxo_coefficients$type <- "Taxonomic"
DESIS_phylo_coefficients <- fread(paste0(out_path, "phylo-log", "_coefficients.csv"))
DESIS_phylo_coefficients$type <- "Phylogenetic"
DESIS_func_coefficients <- fread(paste0(out_path, "func", "_coefficients.csv"))
DESIS_func_coefficients$type <- "Functional"

DESIS <- rbind(DESIS_taxo_coefficients, 
               DESIS_phylo_coefficients,
               DESIS_func_coefficients)
DESIS$sensor <- sensor

DESIS$Axis1 <- DESIS$Axis1*10000
DESIS$Axis2 <- DESIS$Axis2*10000
DESIS$Axis3 <- DESIS$Axis3*10000

# EMIT
sensor <- "EMIT"
out_path <- paste0(root_path, "/06_model_beta/SW_beta_PLS/", sensor, "/", sensor, "_")
EMIT_taxo_coefficients <- fread(paste0(out_path, "taxo", "_coefficients.csv"))
EMIT_taxo_coefficients$type <- "Taxonomic"
EMIT_phylo_coefficients <- fread(paste0(out_path, "phylo-log", "_coefficients.csv"))
EMIT_phylo_coefficients$type <- "Phylogenetic"
EMIT_func_coefficients <- fread(paste0(out_path, "func", "_coefficients.csv"))
EMIT_func_coefficients$type <- "Functional"

EMIT <- rbind(EMIT_taxo_coefficients, 
              EMIT_phylo_coefficients,
              EMIT_func_coefficients)
EMIT$sensor <- sensor

#EMIT$Axis1 <- EMIT$Axis1/10000
#EMIT$Axis2 <- EMIT$Axis2/10000
#EMIT$Axis3 <- EMIT$Axis3/10000

# Merge
coefficients <- rbind(DESIS, EMIT)
coefficients <- coefficients[bands != "(Intercept)",]
coefficients <- coefficients[,-1]
coefficients <- melt(coefficients, 
                     id.vars = c("bands", "type", "sensor"),
                     measure.vars = c("Axis1", "Axis2", "Axis3"),
                     variable.name = "Axes",
                     value.name = "Coefficients")

# Summary 
coef_summary <- coefficients[, .(Coefficients = mean(Coefficients), 
                                 Coefficients_sd = sd(Coefficients)),
                             by = c("bands", "type", "sensor", "Axes")]
coef_summary$bands <- as.numeric(coef_summary$bands)
add_na <- data.table(bands = c(1300, 1300, 1300, 1800, 1800, 1800,
                               1300, 1300, 1300, 1800, 1800, 1800,
                               1300, 1300, 1300, 1800, 1800, 1800),
                     sensor = rep("EMIT", 18),
                     type = c(rep("Taxonomic", 6),
                              rep("Phylogenetic", 6), 
                              rep("Functional", 6)),
                     Axes = c("Axis1", "Axis2", "Axis3",
                              "Axis1", "Axis2", "Axis3",
                              "Axis1", "Axis2", "Axis3"),
                     Coefficients = rep(NA, 18),
                     Coefficients_sd = rep(NA, 18))

coef_summary <- rbind(coef_summary, add_na)
coef_summary <- coef_summary[order(bands)]

coef_summary$Axes <- factor(coef_summary$Axes, labels = c("Axis 1",
                                                          "Axis 2",
                                                          "Axis 3"))

coef_summary$type <- as.factor(coef_summary$type)
coef_summary$type <- factor(coef_summary$type, levels = c("Taxonomic",
                                                          "Phylogenetic",
                                                          "Functional"))

#-------------------------------------------------------------------------------
# Plot details

col <- c("#1b7837", "#762a83")

th_black <- theme_bw() + 
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.margin = margin(4, 4, 0, 1, "pt"),
        legend.position= c("top"), 
        legend.direction = "horizontal", 
        legend.background = element_rect(fill = "transparent"), 
        legend.box.background = element_blank(),
        strip.background = element_rect(color="black", 
                                        fill="black", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "white"))

th_trans <- dark_theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #legend.position = c(0.90, 0.90), 
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(colour="white"),
        legend.text = element_text(color = "white"),
        axis.text = element_text(color="white"),
        plot.background = element_rect(colour = "transparent",
                                       fill = "transparent",
                                       linewidth = 0))

th <- th_black


#-------------------------------------------------------------------------------
# Get plot

plot <- ggplot(coef_summary) +
  geom_ribbon(aes(x = bands, 
                  ymin = Coefficients - Coefficients_sd, 
                  ymax = Coefficients + Coefficients_sd, 
                  fill = sensor), 
              alpha = 0.25) +
  geom_line(aes(x = bands, 
                y = Coefficients, 
                colour = sensor), 
            linewidth = 0.5) +
  coord_cartesian(xlim = c(440, 2425), ylim = c(-19, 19), expand = FALSE) +
  scale_colour_manual(name = "", values = col) +
  scale_fill_manual(name = "", values = col) +
  xlab("Wavelength (nm)") +
  ylab("PLS coefficients") +
  #labs(subtitle = "PLS coefficients") +
  facet_grid(Axes ~ type) +
  th

plot

jpeg(paste0(root_path, "/figures/PLSR_coefficients.jpeg"), 
     width = 210, 
     height = 110, 
     pointsize = 12, 
     units = "mm", 
     res = 300, 
     bg = "transparent")

plot

dev.off()




