#'------------------------------------------------------------------------------
#' @title PCA and functional tree of functional traits
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' @Libraries
#-------------------------------------------------------------------------------

library(data.table)
library(ape)
library(picante)
library(bit64)
library(ggplot2)

#'------------------------------------------------------------------------------
#' @Read files
#-------------------------------------------------------------------------------

#root_path <- "/media/antonio/Work/MMDB"
root_path <- "G:/Projects/MMDB"

fia <- readRDS(paste0(root_path, "/FIA/01_FIA-reshape/fia_reshape.rds"))
SPCD <- fread(paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv"))


#'------------------------------------------------------------------------------
#' @PCA analysis
#-------------------------------------------------------------------------------

# Run PCA similar to 01-beta_functional.R function

# Merge inventory with taxonomic information
fia_summary <- merge(fia, SPCD[, c("SPCD", "Taxon2")], 
                     by = "SPCD", all.x = TRUE, all.y = FALSE)

# Summary of species and families
species <- fia[, .N, by = "SPCD"]
species <- species[, 1]
species <- merge(species, 
                 SPCD[, c("SPCD", "Taxon2", 
                          "HTq", "SLq", "WD", 
                          "shade", "drought", "waterlogging", 
                          "LMA", "N")],
                 by = "SPCD",
                 all.x = TRUE,
                 all.y = FALSE)

species <- as.data.frame(species)
species_tips <- species$Taxon2

# Run PCA
pca <- rda(species[, c(3:10)], scale = TRUE)
traits <- as.data.table(scores(pca, choices = c(1, 2, 3), display = "species"))
traits <- cbind(Traits = c("Tree height", "Slenderness", "Wood density", "Shade", "Drought", "Waterlogging", "LMA", "Nitrogen"),
                traits)
val <- as.data.table(scores(pca, choices = c(1, 2, 3), display = "sites"))
val <- cbind(Taxon2 = species$Taxon2, val)
val <- merge(SPCD[, c("Taxon2", "Taxon", "growth_form", "mycorrhizal_association")],
             val,
             by = "Taxon2", all.x = FALSE, all.y = TRUE)

val[growth_form == "gymnosperm", growth_form := "Gymnosperm"]
val[growth_form == "angiosperm", growth_form := "Angiosperms"]
val[mycorrhizal_association == "AM", mycorrhizal_association := "Arbuscular mycorrhizal"]
val[mycorrhizal_association == "ECM", mycorrhizal_association := "Ectomycorrhizal"]
val[mycorrhizal_association == "ERM", mycorrhizal_association := "Ectomycorrhizal"]

# Proportion of variance explained
eig_vals <- pca$CA$eig
var_explained <- (eig_vals / sum(eig_vals))*100
cum_var_explained <- cumsum(var_explained)
var_explained <- data.frame(
  PC = factor(paste0("PC", seq_along(var_explained)), levels = paste0("PC", seq_along(var_explained))),
  Eigenvalue = eig_vals,
  Variance = var_explained,
  Cumulative = cum_var_explained
)

#'------------------------------------------------------------------------------
#' @Plot
#-------------------------------------------------------------------------------

th_black <- theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill= "transparent"),
        panel.background = element_rect(fill= "transparent"),
        rect = element_rect(fill = "transparent"),
        plot.margin = unit(c(4, 4, 0.5, 0.5), "mm"),
        legend.position= c("top"), 
        strip.background = element_rect(color="black", 
                                        fill="black", 
                                        linewidth=1.5, 
                                        linetype="solid"),
        strip.text = element_text(color = "white"))


pc12 <- ggplot(val, aes(x = PC1, y = PC2)) +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_point(aes(colour = mycorrhizal_association, 
                 fill = growth_form,
                 shape = mycorrhizal_association), size = 3, alpha = 0.75) +
  geom_segment(data = traits,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black") +
  geom_text(data = traits,
            aes(x = PC1, y = PC2, label = Traits),
            color = "black", vjust = 1.5) +
  scale_fill_manual("Plant group", values = c("#ca562c", "#3d5941")) +
  scale_colour_manual("Mycorrhizal association", values = c("grey25", "grey25")) +
  scale_shape_manual("Mycorrhizal association", values = c(21, 24)) +
  coord_fixed(xlim = c(-2.2, 2.2), ylim = c(-2, 2), expand = FALSE) +
  xlab("PC1 (24.28%)") + ylab("PC2 (17.88%)") +
  th_black

pc13 <- ggplot(val, aes(x = PC1, y = PC3)) +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_point(aes(colour = mycorrhizal_association, 
                 fill = growth_form,
                 shape = mycorrhizal_association), size = 3, alpha = 0.75) +
  geom_segment(data = traits,
               aes(x = 0, y = 0, xend = PC1, yend = PC3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black") +
  geom_text(data = traits,
            aes(x = PC1, y = PC3, label = Traits),
            color = "black", vjust = 1.5) +
  scale_fill_manual("Plant group", values = c("#ca562c", "#3d5941")) +
  scale_colour_manual("Mycorrhizal association", values = c("grey25", "grey25")) +
  scale_shape_manual("Mycorrhizal association", values = c(21, 24)) +
  coord_fixed(xlim = c(-2.2, 2.2), ylim = c(-1.65, 1.65), expand = FALSE) +
  xlab("PC1 (24.28%)") + ylab("PC3 (16.32%)") +
  th_black

pc23 <- ggplot(val, aes(x = PC2, y = PC3)) +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_point(aes(colour = mycorrhizal_association, 
                 fill = growth_form,
                 shape = mycorrhizal_association), size = 3, alpha = 0.75) +
  geom_segment(data = traits,
               aes(x = 0, y = 0, xend = PC2, yend = PC3),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black") +
  geom_text(data = traits,
            aes(x = PC2, y = PC3, label = Traits),
            color = "black", vjust = 1.5) +
  scale_fill_manual("Plant group", values = c("#ca562c", "#3d5941")) +
  scale_colour_manual("Mycorrhizal association", values = c("grey25", "grey25")) +
  scale_shape_manual("Mycorrhizal association", values = c(21, 24)) +
  coord_fixed(xlim = c(-2.2, 2.2), ylim = c(-1.65, 1.65), expand = FALSE) +
  xlab("PC2 (17.88%)") + ylab("PC3 (16.32%)") +
  th_black

variance <- ggplot(var_explained, aes(x = PC)) +
  geom_bar(aes(y = Variance), stat = "identity", fill = "#05445e") +
  geom_line(aes(y = Cumulative, group = 1), color = "#189AB4", size = 1) +
  geom_point(aes(y = Cumulative), color = "#189AB4", size = 2) +
  ylab("Proportion of Variance Explained") + 
  scale_y_continuous(lim = c(0, 100), expand = c(0, 0)) +
  th_black

