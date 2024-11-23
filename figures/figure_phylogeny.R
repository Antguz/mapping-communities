#'------------------------------------------------------------------------------
#' @title Codes to visualize the phylogenetic tree
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(ggplot2)
library(data.table)
library(picante)
library(Taxonstand)
library(phytools)
library(bit64)
library(ggdark)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(tidytree)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Read information

# Folder root
root_path <- "/media/antonio/Work/MMDB"
#root_path <- "/media/antonio/Extreme_Pro/Projects/DESIS_mapping"

# Get tree
phylo <- readRDS(paste0(root_path, "/03_diversity/master_matched_phylo.rds"))
features <- fread(paste0(root_path, "/00_FIA/SPCD_meta.csv"))
#features <- features[, c("Taxon2", "growth_form", "mycorrhizal_association", "leaf_habit")]

frame <- data.table(species = phylo$tip.label)
frame$ID <- 1:nrow(frame)
frame <- merge(frame, features, by.x= "species", by.y = "Taxon2", all.x = TRUE, all.y = FALSE)
frame <- frame[order(ID)]
frame <- unique(frame)
frame[growth_form == "gymnosperm", growth_form := "Gymnosperm"]
frame[growth_form == "angiosperm", growth_form := "Angiosperms"]
frame[mycorrhizal_association == "AM", mycorrhizal_association := "Arbuscular mycorrhizal"]
frame[mycorrhizal_association == "ECM", mycorrhizal_association := "Ectomycorrhizal"]
frame[mycorrhizal_association == "ERM", mycorrhizal_association := "Ectomycorrhizal"]

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
th <- th_black

line <- "black"
text <- "black"

#'------------------------------------------------------------------------------
#' Plot

d <- data.frame(node = nodeid(phylo, frame$species),
                growth_form = frame$growth_form,
                mycorrhizal_association = frame$mycorrhizal_association)

d$node <- as.numeric(d$node)
tree <- full_join(phylo, d, by = 'node')
tree <- groupClade(tree, .node = c(244, 452))

frame <- frame[order(-ID)]


p <- ggtree(tree, aes(color=group)) %<+% frame +
  #geom_tippoint(data = tree, aes(fill = mycorrhizal_association), shape = 21, size = 1) +  
  scale_color_manual("Plant group", values = c("#ca562c", "#3d5941")) +
  scale_fill_manual("Mycorrhizal association", values = c("grey75", "grey25")) +
  new_scale_color() +
  new_scale_fill() +
  geom_facet(panel = "Tree height", 
             data = frame, 
             geom = geom_col, 
             aes(x = frame$HTq, 
                 color = frame$mycorrhizal_association, 
                 fill = frame$growth_form),
             orientation = 'y') +
  geom_facet(panel = "Slenderness", 
             data = frame, 
             geom = geom_col, 
             aes(x = frame$SLq, 
                 color = frame$mycorrhizal_association, 
                 fill = frame$growth_form), 
             orientation = 'y') + 
  geom_facet(panel = "Wood density", 
             data = frame, 
             geom = geom_col, 
             aes(x = frame$WOOD_SPGR_GREENVOL_DRYWT, 
                 color = frame$mycorrhizal_association, 
                 fill = frame$growth_form), 
             orientation = 'y') +
  geom_facet(panel = "Moisture content",
            data = frame,
            geom = geom_col,
            aes(x = frame$MC_PCT_GREEN_WOOD,
                color = frame$mycorrhizal_association,
                fill = frame$growth_form),
            orientation = 'y') +
  geom_facet(panel = "Carbon fraction",
            data = frame,
            geom = geom_col,
            aes(x = frame$CARBON_RATIO_LIVE,
                color = frame$mycorrhizal_association,
                fill = frame$growth_form),
            orientation = 'y') +
  geom_facet(panel = "Shade", 
             data = frame, 
             geom = geom_col, 
             aes(x = frame$shade, 
                 color = frame$mycorrhizal_association, 
                 fill = frame$growth_form), 
             orientation = 'y') +
  geom_facet(panel = "Drought", 
             data = frame, 
             geom = geom_col, 
             aes(x = frame$drought, 
                 color = frame$mycorrhizal_association, 
                 fill = frame$growth_form), 
             orientation = 'y') +
  geom_facet(panel = "Waterlogging", 
             data = frame, 
             geom = geom_col, 
             aes(x = frame$waterlogging, 
                 color = frame$mycorrhizal_association, 
                 fill = frame$growth_form), 
             orientation = 'y') +
  geom_facet(panel = "Leaf area",
            data = frame,
            geom = geom_col,
            aes(x = log10(frame$leaf_area),
                color = frame$mycorrhizal_association,
                fill = frame$growth_form),
            orientation = 'y') +
  geom_facet(panel = "LMA", 
             data = frame, 
             geom = geom_col, 
             aes(x = frame$LMA, 
                 color = frame$mycorrhizal_association, 
                 fill = frame$growth_form), 
             orientation = 'y') +
  geom_facet(panel = "Nitrogen", 
             data = frame, 
             geom = geom_col, 
             aes(x = frame$N, 
                 color = frame$mycorrhizal_association, 
                 fill = frame$growth_form), 
             orientation = 'y') +
  scale_fill_manual("Plant group", values = c("#ca562c", "#3d5941")) +
  scale_color_manual("Mycorrhizal association", values = c("grey75", "grey25")) +
  theme_tree2(base_size = 12) +
  th

p

png(paste0(root_path, "/figures/phylogenies_trans.png"), 
    width = 300, 
    height = 200, 
    units = "mm", 
    res = 600, 
    bg = "transparent", 
    pointsize = 12)
p
dev.off()

png(paste0(root_path, "/figures/phylogenies.png"), 
    width = 300, 
    height = 200, 
    units = "mm", 
    res = 600, 
    bg = "white", 
    pointsize = 12)
p
dev.off()
