#'------------------------------------------------------------------------------
#' @title Visualize performance of occurrence of plant groups
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(bit64)
library(data.table)
library(rcartocolor)
library(viridis)
library(ggplot2)
library(ggdark)

#'------------------------------------------------------------------------------
#' @Read-performance files
#-------------------------------------------------------------------------------

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"

# Categorical 
categorical <- fread(paste0(root_path, "/08_occurence_performance/performance_summary_occurance2.csv"))
categorical$Model <- "Occurrence"
#cwm <- fread(paste0(root_path, "/08_CWM_performance/performance_summary_CWM.csv"))
#cwm$Model <- "CWM"

data <- rbind(categorical)
data <- data[dataset == "testing",]

data[type == "taxo", type := "Taxonomic"]
data[type == "phylo-log", type := "Phylogenetic"]
data[type == "func", type := "Functional"]
data$type <- as.factor(data$type)
data$type <- factor(data$type, levels = c("Taxonomic", "Phylogenetic", "Functional"))

data <- data[order(category, plants),]
data$plants <- as.factor(data$plants)
data$plants <- factor(data$plants, levels = rev(c("Altingiaceae", "Betulaceae", "Cupressaceae", "Fagaceae", "Juglandaceae",
                                             "Magnoliaceae", "Nyssaceae", "Oleaceae", "Pinaceae", "Rosaceae",
                                             "Salicaceae", "Sapindaceae", "Ulmaceae", "Abies", "Acer",
                                             "Betula", "Carya", "Fagus", "Fraxinus", "Liquidambar",
                                             "Liriodendron", "Nyssa", "Picea", "Pinus", "Populus",
                                             "Prunus", "Quercus", "Ulmus",
                                             "Angiosperm", "Gymnosperm",
                                             "Arbuscular mycorrhizal", "Ectomycorrhizal",
                                             "Broadleaf", "Evergreen",
                                             "Simple leaf", "Compound leaf", "Needle leaf")))

unique(data$plants)

#-------------------------------------------------------------------------------
# Plot details

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

th_trans <- dark_theme_bw(base_size = 14) +
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
#th <- th_black
line <- "black"
text <- "black"

#-------------------------------------------------------------------------------
#' @Plot 

unique(data$Model)

p <- ggplot(data[Model == "Occurrence"]) +
  geom_point(aes(x = perRMSE, 
                 y = plants, 
                 colour = type,
                 fill = type,
                 shape = type),
             size = 2) +
  coord_cartesian(xlim = c(0.04, 0.5), expand = TRUE) +
  scale_shape_manual("Dimension of biodiversity", values = c(21, 22, 23)) +
  scale_fill_manual("Dimension of biodiversity", values = c("#af8dc3", "#d8b365", "#7fbf7b")) +
  scale_colour_manual("Dimension of biodiversity", values = c("#762a83", "#8c510a", "#1b7837")) +
  facet_grid(. ~ sensor) +
  xlab("%RMSE") + ylab("Occurence") +
  th

png(paste0(root_path, "/figures/Occurence_performance.png"), 
    width = 210, 
    height = 200, 
    units = "mm", 
    res = 600, 
    bg = "transparent", 
    pointsize = 14)
p
dev.off()

