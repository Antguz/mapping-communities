#'------------------------------------------------------------------------------
#' @title Non-metric Multidimensional Scaling
#'------------------------------------------------------------------------------

#' @description It creates a new ordination of species abundance of plots
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(doSNOW)
library(foreach)
library(ecodist)
library(bigmds)
library(ape)
library(ggplot2)
library(bit64)
library(ggrepel)
#library(BiodiversityR)
options(scipen = 100)
library(viridis)
library(ggpointdensity)
library(ggpubr)

#'------------------------------------------------------------------------------
#' Arguments
#' @param stand
#' @param out_path

root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/media/antonio/Extreme Pro/Projects/DESIS_mapping"
#root_path <- "E:/Projects/DESIS_mapping"
stand <- readRDS(paste0(root_path, "/03_diversity/master_taxonomic_DRYBIO.rds"))
out_path <- paste0(root_path, "/05_ordination")

#'------------------------------------------------------------------------------
#' @analysis  

#-------------------------------------------------------------------------------
# Selecting the NMDS for taxonomic or phylogenetic ordinations

NMDS_taxonomic <- readRDS(paste0(root_path, "/05_ordination/taxonomic/NMDS_taxonomic_DRYBIO_k3.rds"))
NMDS_phylogenetic <- readRDS(paste0(root_path, "/05_ordination/beta/NMDS_beta_DRYBIO_k3.rds"))

# Plot for inspection
plot(NMDS_taxonomic, choices = c(1, 2))
plot(NMDS_taxonomic, choices = c(1, 3))
plot(NMDS_taxonomic, choices = c(2, 3))

# Plot for inspection
plot(NMDS_phylogenetic, choices = c(1, 2))
plot(NMDS_phylogenetic, choices = c(1, 3))
plot(NMDS_phylogenetic, choices = c(2, 3))

#-------------------------------------------------------------------------------
# Get stress plot
stressplot(NMDS_taxonomic)
stressplot(NMDS_phylogenetic)

#-------------------------------------------------------------------------------
# Get weighted averages scores of species
species_scores_taxonomic <- wascores(NMDS_taxonomic$points, stand, expand = FALSE)
species_names_taxonomic <- rownames(species_scores_taxonomic)
species_scores_taxonomic <- cbind(Taxon2 = species_names_taxonomic, as.data.table(species_scores_taxonomic))
saveRDS(species_scores_taxonomic, paste0(out_path, "/species_scores_taxonomic_DRYBIO.rds"))

species_scores_phylogenetic <- wascores(NMDS_phylogenetic$points, stand, expand = FALSE)
species_names_phylogenetic <- rownames(species_scores_phylogenetic)
species_scores_phylogenetic <- cbind(Taxon2 = species_names_phylogenetic, as.data.table(species_scores_phylogenetic))
saveRDS(species_scores_phylogenetic, paste0(out_path, "/species_scores_beta_DRYBIO.rds"))

#-------------------------------------------------------------------------------
# Get dimensions
plot_names_taxonomic <- as.integer64(rownames(NMDS_taxonomic$points))
dimensions_plots_taxonomic <- as.data.table(NMDS_taxonomic$points)
dimensions_plots_taxonomic$PLT_CN <- plot_names_taxonomic
dimensions_plots_taxonomic$goodness <- goodness(NMDS_taxonomic)
saveRDS(dimensions_plots_taxonomic, paste0(out_path, "/dimensions_plots_taxonomic_DRYBIO.rds"))

plot_names_phylogenetic <- as.integer64(rownames(NMDS_phylogenetic$points))
dimensions_plots_phylogenetic <- as.data.table(NMDS_phylogenetic$points)
dimensions_plots_phylogenetic$PLT_CN <- plot_names_phylogenetic
dimensions_plots_phylogenetic$goodness <- goodness(NMDS_phylogenetic)
saveRDS(dimensions_plots_phylogenetic, paste0(out_path, "/dimensions_plots_beta_DRYBIO.rds"))
