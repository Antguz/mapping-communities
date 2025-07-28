#'------------------------------------------------------------------------------
#' @title Codes to visualize the VIP from the models (Fig 4)
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
DESIS_taxo_vip <- fread(paste0(out_path, "taxo", "_vip.csv"))
DESIS_taxo_vip$type <- "Taxonomic"
DESIS_phylo_vip <- fread(paste0(out_path, "phylo", "_vip.csv"))
DESIS_phylo_vip$type <- "Phylogenetic"
DESIS_func_vip <- fread(paste0(out_path, "func", "_vip.csv"))
DESIS_func_vip$type <- "Functional"

DESIS <- rbind(DESIS_taxo_vip, 
               DESIS_phylo_vip,
               DESIS_func_vip)
DESIS$sensor <- sensor

DESIS <- melt(DESIS, id.vars = c("sensor", "type", "iteration", "Axes"),
              variable.name = "wavelength", value.name = "VIP")

DESIS <- DESIS[wavelength != "(Intercept)",]
DESIS$wavelength <- as.numeric(as.character(DESIS$wavelength))

# EMIT
sensor <- "EMIT"
out_path <- paste0(root_path, "/Models/01_model_beta/", sensor, "/", sensor, "_")
EMIT_taxo_vip <- fread(paste0(out_path, "taxo", "_vip.csv"))
EMIT_taxo_vip$type <- "Taxonomic"
EMIT_phylo_vip <- fread(paste0(out_path, "phylo", "_vip.csv"))
EMIT_phylo_vip$type <- "Phylogenetic"
EMIT_func_vip <- fread(paste0(out_path, "func", "_vip.csv"))
EMIT_func_vip$type <- "Functional"

EMIT <- rbind(EMIT_taxo_vip, 
              EMIT_phylo_vip,
              EMIT_func_vip)
EMIT$sensor <- sensor

EMIT <- melt(EMIT, id.vars = c("sensor", "type", "iteration", "Axes"),
             variable.name = "wavelength", value.name = "VIP")

EMIT <- EMIT[wavelength != "(Intercept)",]
EMIT$wavelength <- as.numeric(as.character(EMIT$wavelength))

# Merge
vip <- rbind(DESIS, EMIT)

# Summary 
vip_summary <- vip[, .(VIP = mean(VIP), 
                                 vip_sd = sd(VIP)),
                             by = c("wavelength", "type", "sensor", "Axes")]
vip_summary$wavelength <- as.numeric(vip_summary$wavelength)
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
                     VIP = rep(NA, 18),
                     vip_sd = rep(NA, 18))

vip_summary <- rbind(vip_summary, add_na)
vip_summary <- vip_summary[order(wavelength)]

vip_summary$Axes <- factor(vip_summary$Axes, labels = c("Axis 1",
                                                          "Axis 2",
                                                          "Axis 3"))

vip_summary$type <- as.factor(vip_summary$type)
vip_summary$type <- factor(vip_summary$type, levels = c("Taxonomic",
                                                          "Phylogenetic",
                                                          "Functional"))

vip_summary$sensor <- as.factor(vip_summary$sensor)
vip_summary$sensor <- factor(vip_summary$sensor, levels = c("DESIS",
                                                        "EMIT"))

#-------------------------------------------------------------------------------
# Plot details

#col <- c("#b35806", "#542788")
col <- c("#e08214", "#8073ac")


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

plot <- ggplot(vip_summary) +
  geom_ribbon(aes(x = wavelength, 
                  ymin = VIP - vip_sd, 
                  ymax = VIP + vip_sd, 
                  fill = sensor), 
              alpha = 0.25) +
  geom_line(aes(x = wavelength, 
                y = VIP, 
                colour = sensor), 
            linewidth = 0.5) +
  geom_hline(yintercept = 1, colour = "grey", linetype = "dotted") + 
  coord_cartesian(xlim = c(440, 2430), ylim = c(0, 4.5), expand = FALSE) +
  scale_y_continuous(breaks = c(0, 2, 4)) +
  scale_colour_manual(name = "Sensor", values = col) +
  scale_fill_manual(name = "Sensor", values = col) +
  xlab("Wavelength (nm)") +
  ylab("Variable of Importance in Projection") +
  #labs(subtitle = "PLS vip") +
  facet_grid(Axes ~ type) +
  th

plot

# Export JPEG

jpeg(paste0(root_path, "/figures/Manuscript_black/PLSR_vip.jpeg"), 
     width = 210, 
     height = 110, 
     pointsize = 12, 
     units = "mm", 
     res = 300, 
     bg = "white")

plot

dev.off()


# Export PNG

png(paste0(root_path, "/figures/Manuscript_white/PLSR_vip.png"), 
    width = 210, 
    height = 110, 
    pointsize = 12, 
    units = "mm", 
    res = 300, 
    bg = "transparent")

plot

dev.off()



