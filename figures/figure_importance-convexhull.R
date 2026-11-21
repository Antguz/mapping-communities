#'------------------------------------------------------------------------------
#' @title Relationship between relative importance and occurrance
#'------------------------------------------------------------------------------

#' @description It is used to evaluate the relationship between relative importance
#' and occurrence

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(bit64)
library(ape)
library(geometry)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(vegan)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' @Working_directory
#-------------------------------------------------------------------------------

root_path <- "/media/guzman/Work/MMDB"

#'------------------------------------------------------------------------------
#' @PCA_of_climate_and_topography
#-------------------------------------------------------------------------------

# Select file of climate and topography
lineages_path <- paste0(root_path, "/Models/03_model_occurrence")
variables <- readRDS(paste0(lineages_path, "/LINEAGES_variables.rds"))

# 
pca <- rda(variables[, 5:12], scale = TRUE)
scores <- as.data.table(scores(pca, choices = c(1, 2, 3), display = "sites"))
scores <- cbind(variables[, "PLT_CN"], scores)

#'------------------------------------------------------------------------------
#' @Convexhull_volume
#-------------------------------------------------------------------------------

estimate_hull_volume <- function(group, scores) {
  
  # Clean group
  observed <- as.data.table(copy(group))
  group_name <- colnames(observed)[2]
  colnames(observed)[2] <- "occurrance"
  observed <- subset(observed, occurrance == 1)
  
  # Merge group with pca
  observed <- merge(observed, scores, by = "PLT_CN", all.x = TRUE, all.y = FALSE)
  observed_scores <- observed[, 3:5]
  observed_scores <- as.matrix(na.exclude(observed_scores))
  
  # Compute convex hull
  hull <- tryCatch(
    convhulln(observed_scores, options = "FA"),
    error = function(e) NULL
  )
  
  return(data.table(group = group_name, hull_vol = hull$vol))
  
}

#'------------------------------------------------------------------------------
#' @Forest_type
#-------------------------------------------------------------------------------

fortyp_path <- paste0(root_path, "/Models/03_models_fortyp")

# Occurrance and variables
FORTYPCD_wide <- readRDS(paste0(fortyp_path, "/FORTYPCD_occurrance.rds"))
FORTYPCD_ordination <- readRDS(paste0(fortyp_path, "/FORTYPCD_variables.rds"))

# Estimate volume
group_cols <- setdiff(names(FORTYPCD_wide), "PLT_CN")

FORTYPCD_list <- mclapply(X = group_cols,
                          function(col) {
                            group_sub <- FORTYPCD_wide[, .(PLT_CN, tmp = get(col))]
                            setnames(group_sub, "tmp", col)
                            estimate_hull_volume(group_sub, scores)
                          },
                          mc.cores = 25)

FORTYPCD_list <- rbindlist(FORTYPCD_list, fill = TRUE)

#'------------------------------------------------------------------------------
#' @Plant_lineages
#-------------------------------------------------------------------------------

lineages_path <- paste0(root_path, "/Models/03_model_occurrence")

LINEAGES_wide <- readRDS(paste0(lineages_path, "/LINEAGES_occurrance.rds"))
LINEAGES_ordination <- readRDS(paste0(lineages_path, "/LINEAGES_variables.rds"))

# Estimate volume
group_cols <- setdiff(names(LINEAGES_wide), "PLT_CN")

LINEAGES_list <- mclapply(X = group_cols,
                          function(col) {
                            group_sub <- LINEAGES_wide[, .(PLT_CN, tmp = get(col))]
                            setnames(group_sub, "tmp", col)
                            estimate_hull_volume(group_sub, scores)
                          },
                          mc.cores = 25)

LINEAGES_list <- rbindlist(LINEAGES_list, fill = TRUE)

#'------------------------------------------------------------------------------
#' @Importance_FORTYPCD
#-------------------------------------------------------------------------------

# Categorical 
data <- fread(paste0(fortyp_path, "/taxo_altogether_model-FORTYP.csv"))
data <- data[, .SD, .SDcols = c(1:3, 27:37)]
colnames(data)[4:14] <- gsub("imp_","", colnames(data)[4:14])
colnames(data)[4:14] <- c("Axis 1", "Axis 2", "Axis 3", "MAT", "TS", "TAR", "MAP", "PS", "ELE", "SLO", "CTI")

# Estimate relative importance
total <- rowSums(data[,.SD, .SDcols = 4:ncol(data)])
data[, 4:ncol(data)] <- data[,.SD, .SDcols = 4:ncol(data)]/total
data <- subset(data, dataset == "training")
fortyp_axes <- rowSums(data[, 4:6])
fortyp_axes <- data.table(group = data$forty, ri_axes = fortyp_axes)

fortyp_axes <- fortyp_axes[, .(axes_mean = mean(ri_axes), 
                               axes_sd = sd(ri_axes)),
                           by = group]

fortyp_axes <- merge(fortyp_axes, FORTYPCD_list, by = "group", all.x = TRUE, all.y = FALSE)

#'------------------------------------------------------------------------------
#' @Importance_LINEAGES
#-------------------------------------------------------------------------------

# Categorical 
data <- fread(paste0(lineages_path, "/phylo_altogether_model-occurrence.csv"))
data <- data[, .SD, .SDcols = c(1:3, 27:37)]
colnames(data)[4:14] <- gsub("imp_","", colnames(data)[4:14])
colnames(data)[4:14] <- c("Axis 1", "Axis 2", "Axis 3", "MAT", "TS", "TAR", "MAP", "PS", "ELE", "SLO", "CTI")
data <- data[plant != "ECM",]
data <- data[plant != "AM",]

# Estimate relative importance
total <- rowSums(data[,.SD, .SDcols = 4:ncol(data)])
data[, 4:ncol(data)] <- data[,.SD, .SDcols = 4:ncol(data)]/total
data <- subset(data, dataset == "training")
lineages_axes <- rowSums(data[, 4:6])
lineages_axes <- data.table(group = data$plant, ri_axes = lineages_axes)

lineages_axes <- lineages_axes[, .(axes_mean = mean(ri_axes), 
                               axes_sd = sd(ri_axes)),
                           by = group]

lineages_axes <- merge(lineages_axes, LINEAGES_list, by = "group", all.x = TRUE, all.y = FALSE)
families <- c("Altingiaceae", "Betulaceae", "Cupressaceae", 
              "Fagaceae", "Juglandaceae", "Magnoliaceae", 
              "Nyssaceae", "Oleaceae", "Pinaceae", "Rosaceae",
              "Salicaceae", "Sapindaceae", "Ulmaceae")
families_axes <- lineages_axes[group %in% families,]

genus <- c("Abies", "Acer", "Betula", "Carya", "Fagus", 
           "Fraxinus", "Liquidambar", "Liriodendron", 
           "Nyssa", "Picea", "Pinus", "Populus",
           "Prunus", "Quercus", "Ulmus")
genus_axes <- lineages_axes[group %in% genus,]

#'------------------------------------------------------------------------------
#' @Plot
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Plot details

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(0.5, 4, 0.5, 0.5), "mm"),
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
        plot.margin = unit(c(0.5, 4, 0.5, 0.5), "mm"),
        strip.background = element_rect(color="white", 
                                        fill="white", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "black"))

th <- th_trans
line <- "white"
text <- "white"

th <- th_black
line <- "black"
text <- "black"

#--------------------------------------------
# Forest types

A <- ggplot(fortyp_axes) +
  geom_errorbar(aes(x = hull_vol,
                    ymin = axes_mean - axes_sd,
                    ymax = axes_mean + axes_sd),
                width = 0,
                colour = "darkred",
                alpha = 0.6) +
  geom_point(aes(x = hull_vol, y = axes_mean), 
             shape = 21, 
             size = 3,
             fill = "darkred",
             colour = "white") +
  geom_text_repel(aes(x = hull_vol, y = axes_mean, label = group),
                  size = 3,
                  colour = "grey25",
                  box.padding = 0.35,
                  point.padding = 0.2,
                  segment.color = "grey70",
                  max.overlaps = 10,
                  min.segment.length = 0) +
  geom_smooth(aes(x = hull_vol, y = axes_mean),
              method = "glm",
              se = TRUE,
              colour = "black") +
  stat_regline_equation(aes(x = hull_vol, y = axes_mean, label = paste(..rr.label..)),
                        label.x.npc = "right",
                        label.y.npc = "bottom") +
  coord_cartesian(ylim = c(0.15, 0.85), xlim = c(0, 0.20)) +
  th +
  labs(subtitle = "Forest types") +
  #lab("Climatic and topographic hypervolume") +
  ylab(" ") +
  xlab(" ")
  #ylab("Relative importance of MDS β-diversity axes")

B <- ggplot(families_axes) +
  geom_errorbar(aes(x = hull_vol,
                    ymin = axes_mean - axes_sd,
                    ymax = axes_mean + axes_sd),
                width = 0,
                colour = "darkred",
                alpha = 0.6) +
  geom_point(aes(x = hull_vol, y = axes_mean), 
             shape = 21, 
             size = 3,
             fill = "darkred",
             colour = "white") +
  geom_text_repel(aes(x = hull_vol, y = axes_mean, label = group),
                  size = 3,
                  colour = "grey25",
                  box.padding = 0.35,
                  point.padding = 0.2,
                  segment.color = "grey70",
                  max.overlaps = Inf,
                  min.segment.length = 0) +
  geom_smooth(aes(x = hull_vol, y = axes_mean),
              method = "glm",
              se = TRUE,
              colour = "black") +
  stat_regline_equation(aes(x = hull_vol, y = axes_mean, label = paste(..rr.label..)),
                        label.x.npc = "left",
                        label.y.npc = "top") +
  coord_cartesian(ylim = c(0.15, 0.85), xlim = c(0, 0.20)) +
  th +
  labs(subtitle = "Plant lineages: families") +
  #xlab("Climatic and topographic hypervolume") +
  xlab(" ") +
  ylab("Relative importance of MDS β-diversity axes")
  

C <- ggplot(genus_axes) +
  geom_errorbar(aes(x = hull_vol,
                    ymin = axes_mean - axes_sd,
                    ymax = axes_mean + axes_sd),
                width = 0,
                colour = "darkred",
                alpha = 0.6) +
  geom_point(aes(x = hull_vol, y = axes_mean), 
             shape = 21, 
             size = 3,
             fill = "darkred",
             colour = "white") +
  geom_text_repel(aes(x = hull_vol, y = axes_mean, label = group),
                  size = 3,
                  colour = "grey25",
                  box.padding = 0.35,
                  point.padding = 0.2,
                  segment.color = "grey70",
                  max.overlaps = Inf,
                  min.segment.length = 0) +
  geom_smooth(aes(x = hull_vol, y = axes_mean),
              method = "glm",
              se = TRUE,
              colour = "black") +
  stat_regline_equation(aes(x = hull_vol, y = axes_mean, label = paste(..rr.label..)),
                        label.x.npc = "left",
                        label.y.npc = "top") +
  coord_cartesian(ylim = c(0.15, 0.85), xlim = c(0, 0.20)) +
  th +
  labs(subtitle = "Plant lineages: genus") +
  xlab("Climatic and topographic hypervolume") +
  ylab(" ") 
  #ylab("Relative importance of MDS β-diversity axes")

plot <- ggarrange(A, B, C, ncol = 1)

jpeg(paste0(root_path, "/figures/Manuscript_black/importance_hypervolume_new.jpg"), 
     width = 150, 
     height = 230, 
     units = "mm", 
     res = 600, 
     bg = "white",
     #bg = "transparent", 
     pointsize = 14)
plot
dev.off()
