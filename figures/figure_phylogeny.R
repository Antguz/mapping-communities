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
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Read information

# Folder root
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"

# Get tree
phylo <- readRDS(paste0(root_path, "/03_diversity/master_matched_phylo-tree.rds"))
features <- fread(paste0(root_path, "//00_FIA/SPCD_features2.csv"))
features <- features[, c("Taxon2", "Evergreen", "Gymnosperm", "mycorrhizal_association")]

frame <- data.table(species = phylo$tip.label)
frame$ID <- 1:nrow(frame)
frame <- merge(frame, features, by.x= "species", by.y = "Taxon2", all.x = TRUE, all.y = FALSE)
frame <- frame[order(ID)]
frame <- unique(frame)
frame[Gymnosperm == "Y", Gymnosperm := "Gymnosperm"]
frame[Gymnosperm == "N", Gymnosperm := "Angiosperms"]
frame[mycorrhizal_association == "AM", mycorrhizal_association := "Arbuscular mycorrhizal"]
frame[mycorrhizal_association == "ECM", mycorrhizal_association := "Ectomycorrhizal"]

colnames(frame)[4] <- "Group"

p <- ggtree(phylo, layout="circular", branch.length = "branch.length", colour = "white")  

p1 <- p + geom_fruit(data = frame,
               geom = geom_tile,
               mapping = aes(y=species,
                             fill=Group),
               offset = 0.1,
               width=15) +
  scale_fill_manual("Plant group", values = c("#ca562c", "#3d5941"))

p2 <- p1 + 
  new_scale_fill() +
  geom_fruit(data = frame,
             geom = geom_tile,
             mapping = aes(y=species,
                           fill=mycorrhizal_association),
             offset = 0.1,
             width = 15) +
  scale_fill_manual("Mycorrhizal association", values = c("#009392", "#d0587e")) +
  theme(legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        legend.title = element_text(color = "white"),
        legend.text = element_text(color = "white"),
        axis.text=element_text(color = "transparent"), 
        plot.background = element_rect(colour = "transparent", 
                                       fill = "transparent", 
                                       size=0),
        panel.border = element_rect(colour = "black", fill=NA, size=5))

png("figures/phylogenies_black.png", width = 140, height = 90, units = "mm", res = 600, bg = "transparent", pointsize = 12)
p2
dev.off()


