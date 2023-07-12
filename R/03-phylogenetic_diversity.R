#'------------------------------------------------------------------------------
#' @title Beta phylogenetic diversity between communities 
#'------------------------------------------------------------------------------

#' @description It estimates beta phylogenetic diversity between communities.
#' @return Return a distance matrices as .rds files.

################################################################################
#' @warming If you have to many communities this need to be run using HPC
################################################################################

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(picante)
library(Taxonstand)
library(V.PhyloMaker)
library(dplyr)
library(tidyr)
library(phytools)
library(geiger)
library(bit64)

#Set working directory

#'------------------------------------------------------------------------------
#' Source code
setwd("/panfs/jay/groups/17/cavender/guzman/DESIS_mapping/mapping-species")
source("R/beta_rcpp.R")

#'------------------------------------------------------------------------------
#' Arguments
#' @param fia
#' @param SPCD
#' @param domain
#' @param domains_of_interest
#' @param out_path
#' @param threads
#' @param progress

#'------------------------------------------------------------------------------
#' Function
phylogenetic_diversity <- function(fia, 
                                   SPCD, 
                                   domain, 
                                   domains_of_interest = NULL, 
                                   outpath, 
                                   threads, 
                                   progress = TRUE) {
  
  ### Prepare plots (Cleaning) -------------------------------------------------
  
  # Read fia
  fia <- fia[order(PLT_CN, SUBP)]
  
  # Inventories from 2018 until now
  fia <- fia[MEASYEAR >= 2018]
  
  # Merge species code
  fia <- merge(fia, 
               SPCD[, c("SPCD", "PLANTS-Code-new", 
                        "Taxon2", "Genus", "Family", 
                        "Remove")], 
               by = "SPCD", 
               all.x = TRUE, 
               all.y = TRUE)
  
  # Remove dead trees
  fia <- fia[STATUSCD == 1]
  
  # Get forest metrics by TPA
  fia$BA <- pi*(fia$DIA/2)^2
  fia$BA_TPA <- fia$BA*fia$TPA_UNADJ
  fia$DRYBIO_AG_TPA <- fia$DRYBIO_AG*fia$TPA_UNADJ
  
  # Get plots with complete
  complete_plots <- fia[, sum(is.na(DRYBIO_AG_TPA)), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  # Select plots with complete information
  fia[is.na(Remove), Remove := 1]
  
  complete_plots <- fia[, sum(Remove), by = "PLT_CN"]
  complete_plots <- complete_plots[V1 == 0]
  
  fia <- merge(complete_plots[, 1], 
               fia, 
               by = "PLT_CN", 
               all.x = TRUE, 
               all.y = FALSE)
  
  # Summarize
  fia_summary <- fia[, .(ntrees = .N,
                         BA_TPA = sum(BA, na.rm = TRUE), 
                         DRYBIO_AG_TPA = sum(DRYBIO_AG_TPA, na.rm = TRUE)),
                     by = c("PLT_CN", "PLANTS-Code-new", 
                            "Taxon2", "Genus", "Family")]
  fia_summary <- fia_summary[order(PLT_CN)]
  
  # Select plots with proper boundaries within domains
  colnames(domain)[1] <- "Domain"
  
  if(is.null(domains_of_interest)) {
    fia_summary <- merge(domain, fia_summary, by = "PLT_CN", all.x = FALSE, all.y = FALSE)
  } else {
    domain <- subset(domain, Domain == domains_of_interest)
    fia_summary <- merge(domain, fia_summary, by = "PLT_CN", all.x = FALSE, all.y = FALSE)
  }
  
  
  ### Prepare phylo  -----------------------------------------------------------
  # Summary of species and families
  species <- unique(fia_summary[, c("Taxon2", "Genus", "Family")])
  colnames(species) <- c("species", "genus", "family")
  species <- na.exclude(species)
  hypotheses <- phylo.maker(species, scenarios = "S3")
  
  # Resolve Multichotomies
  phylo <- multi2di(hypotheses$scenario.3)
  is.binary.phylo(phylo) #Test for Binary Tree
  is.ultrametric(phylo) #Test if a Tree is Ultrametric
  
  # Transform phylo
  tol = 1e-9
  phylo$edge.length[phylo$edge.length <= 0] <- tol
  is.ultrametric(phylo)
  
  # Scale phylo
  rtrns <- rescale(phylo, "nrate")
  phylo <- rtrns(time=c(0.2, 0.4, 0.6, 0.8), rate=c(2, 4, 8, 16))
  is.ultrametric(phylo)
  
  ### Community matrix ---------------------------------------------------------
  
  # dcast
  master_matrix <- sample2matrix(fia_summary[, c("PLT_CN", "DRYBIO_AG_TPA", "Taxon2")])
  
  # Remove plots without species
  nsp <- apply(master_matrix, 1, function(x) sum(x > 0))
  master_matrix <- master_matrix[nsp != 0, ]
  
  # Standardization of comunities
  master_matrix <- decostand(master_matrix, method = "hellinger") 
  
  ### Get phylogenetic beta diversity ------------------------------------------
  
  #Get match
  matched <- picante::match.phylo.comm(phy = phylo, comm = master_matrix)
  
  cat(paste0("Number of communities:", nrow(master_matrix)))
  
  # Export
  if(is.null(domains_of_interest)) {
    saveRDS(matched$comm, paste0(out_path, "/master_matched_comm-DRYBIO.rds"))
    saveRDS(matched$phy, paste0(out_path, "/master_matched_phylo-tree.rds"))
    
  } else {
    saveRDS(matched$comm, paste0(out_path, "/master_matched_comm-DRYBIO_",domains_of_interest,".rds"))
    saveRDS(matched$phy, paste0(out_path, "/master_matched_phylo-tree_",domains_of_interest,".rds"))
  }
  
  beta <- beta_rcpp(comm = as.matrix(matched$comm), 
                    tree = matched$phy, 
                    threads = threads,
                    progress = progress)
  
  if(is.null(domains_of_interest)) {
    saveRDS(beta, paste0(out_path, "/beta_DRYBIO.rds"))
    
  } else {
    saveRDS(beta, paste0(out_path, "/beta_DRYBIO_",domains_of_interest,".rds"))
  }
  
}

#'------------------------------------------------------------------------------
#' @example 

#root_path <- "F:/Projects/DESIS_mapping"
root_path <- "/media/antonio/Work/Projects/DESIS_mapping"
#root_path <- "/panfs/jay/groups/17/cavender/guzman/DESIS_mapping"
fia <- fread(paste0(root_path, "/00_FIA/trees_eastern_20230427.csv"))
SPCD <- fread(paste0(root_path, "/00_FIA/SPCD.csv"), na.strings="NA")
domain <- fread(paste0(root_path, "/00_FIA/plot_domain.csv"))
domain <- domain[, c("DomainName", "PLT_CN")]
domains_of_interest <- NULL
out_path <- paste0(root_path, "/03_diversity")
threads <- 28
progress = TRUE

domain <- domain[DomainName != "Ozarks Complex",]
domain <- domain[DomainName != "Atlantic Neotropical",]
domain <- domain[DomainName != "Southeast",]

phylogenetic_diversity(fia, 
                       SPCD, 
                       domain, 
                       domains_of_interest = NULL,
                       outpath, 
                       threads, 
                       progress)

domain <- domain[DomainName == "Great Lakes", DomainName := "great_lakes"]
domain <- domain[DomainName == "Northeast", DomainName := "northeast"]
domain <- domain[DomainName == "Appalachians / Cumberland Plateau", DomainName := "appalachians"]
domain <- domain[DomainName == "Mid Atlantic", DomainName := "mid_atlantic"]

domain <- domain[DomainName == "Ozarks Complex", DomainName := "ozarks_complex"]
domain <- domain[DomainName == "Atlantic Neotropical", DomainName := "atlantic_neotropical"]
domain <- domain[DomainName == "Southeast", DomainName := "southeast"]


phylogenetic_diversity(fia, 
                       SPCD, 
                       domain, 
                       domains_of_interest = "northeast",
                       outpath, 
                       threads, 
                       progress)

phylogenetic_diversity(fia, 
                       SPCD, 
                       domain, 
                       domains_of_interest = "great_lakes",
                       outpath, 
                       threads, 
                       progress)

phylogenetic_diversity(fia, 
                       SPCD, 
                       domain, 
                       domains_of_interest = "ozarks_complex",
                       outpath, 
                       threads, 
                       progress)

phylogenetic_diversity(fia, 
                       SPCD, 
                       domain, 
                       domains_of_interest = "appalachians",
                       outpath, 
                       threads, 
                       progress)

phylogenetic_diversity(fia, 
                       SPCD, 
                       domain, 
                       domains_of_interest = "atlantic_neotropical",
                       outpath, 
                       threads, 
                       progress)

phylogenetic_diversity(fia, 
                       SPCD, 
                       domain, 
                       domains_of_interest = "southeast",
                       outpath, 
                       threads, 
                       progress)

phylogenetic_diversity(fia, 
                       SPCD, 
                       domain, 
                       domains_of_interest = "mid_atlantic",
                       outpath, 
                       threads, 
                       progress)
