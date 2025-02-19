#'------------------------------------------------------------------------------
#' @title A script for post processing the MDS ordinations
#'------------------------------------------------------------------------------

#' @description It rotates ordinations according to environmental gradients to 
#' then extract information associated with the projection of species, traits, and
#' environments. It also extract the goodness of fit the communities in the ordinations 
#' and the data for the stressplot.

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
library(bit64)
library(Rcpp)
library(RcppArmadillo)
library(smatr)
library(dynutils)
library(pbmcapply)

#'------------------------------------------------------------------------------
#' Source code

sourceCpp("src/euclidean_distance.cpp")
source("R/auxiliary/MDS_rotate.R")
source("R/auxiliary/aggregate_plants.R")
source("R/auxiliary/voxels.R")

#'------------------------------------------------------------------------------
#' Function

postprocessing <- function(ordination, 
                           beta, 
                           fia_environment, 
                           sp_meta, 
                           species, 
                           fia_traits, 
                           threads, 
                           out_path, 
                           type) {
  
  # ----------------------------------------------------------------------------
  # Read ordination
  ord <- readRDS(ordination)
  comm_names <- as.integer64(rownames(ord))
  
  # ----------------------------------------------------------------------------
  # Read climate and rotate ordination
  
  # Climate
  env <- fread(fia_environment)
  vector_match <- match(comm_names, env$PLT_CN)
  env <- env[vector_match,]
  
  # Rotate
  ord_new <- MDS_rotate(ord, env[, c("MAT", "MAP")], na.rm = TRUE)
  rm(list = c("ord"))
  gc()
  
  #-----------------------------------------------------------------------------
  # Get community ordinations and goodness of fit
  
  # Get community ordinations
  comm <- ord_new
  colnames(comm) <- paste0("Axis", 1:ncol(comm))
  comm <- cbind(data.table(PLT_CN = as.integer64(rownames(comm))), comm)
  
  # Scale from -0.5 to 0.5
  comm[, 2:4] <- scale_uniform(comm[, 2:4], center = 0, max_range = 1)
  fwrite(comm, paste0(out_path, "/MDS_scores_comm_", type, "_total.csv"))
  
  # ----------------------------------------------------------------------------
  # Get species scores
  
  # Read species
  sp <- readRDS(species)
  vector_match <- match(comm_names, as.integer64(rownames(sp)))
  sp <- sp[vector_match,]
  sp1 <- decostand(as.matrix(sp), method = "hellinger")
  
  # Get scores
  sp_scores <- wascores(x = comm[, 2:4], sp1, expand = TRUE)
  sp_scores <- data.table(variable = rownames(sp_scores), sp_scores)
  fwrite(sp_scores, paste0(out_path, "/MDS_species-wascores_", type, "_total.csv"))
  rm(list = c("sp_scores", "sp1"))
  
  # ----------------------------------------------------------------------------
  # Get lineages scores

  # Aggregate species
  aggregation <- aggregate_batch(fread(sp_meta)[,c("Taxon2", "Family", "Genus", "mycorrhizal_association")], 
                                 decostand(sp, "total"))
  
  # Get scores
  sp_scores <- wascores(x = comm[, 2:4], aggregation, expand = TRUE)
  sp_scores <- data.table(variable = rownames(sp_scores), sp_scores)
  fwrite(sp_scores, paste0(out_path, "/MDS_lineages-wascores_", type, "_total.csv"))
  rm(list = c("sp_scores", "aggregation"))
  
  # ----------------------------------------------------------------------------
  # Get traits scores
  
  # Read traits
  traits <- fread(fia_traits)
  vector_match <- match(colnames(sp), traits$Taxon2)
  traits <- traits[vector_match, ]
  meta <- traits[, c("Taxon2", 
                     "HTq", "SLq", "WD", 
                     "shade", "drought", "waterlogging", 
                     "LMA", "N")]
  communities <- decostand(sp, "total")
  meta <- aggregate_batch(meta, communities)
  traits_scores <- wascores(x = comm[, 2:4], meta[,-1], expand = TRUE)
  traits_scores <- data.table(variable = rownames(traits_scores), traits_scores)
  fwrite(traits_scores, paste0(out_path, "/MDS_traits-wascores_", type, "_total.csv"))
  rm(list = c("traits_scores", "meta", "sp", "communities"))
  
  # ----------------------------------------------------------------------------
  # Get environmental scores
  
  # Read traits
  env <- na.exclude(env)
  vector_match <- match(env$PLT_CN, comm$PLT_CN)
  comm2 <- comm[vector_match,]
  env_scores <- wascores(x = comm2[, 2:4], env[,-1], expand = TRUE)
  env_scores <- data.table(variable = rownames(env_scores), env_scores)
  fwrite(env_scores, paste0(out_path, "/MDS_env-wascores_", type, "_total.csv"))
  rm(list = c("env_scores", "comm2", "env", "comm"))
  gc()
  
  # ----------------------------------------------------------------------------
  # Get goodness
  
  # Get distances from lmds
  lmds_distance <- euclidean_distance(as.matrix(ord_new), n_threads = threads)
  colnames(lmds_distance) <- rownames(ord_new)
  rownames(lmds_distance) <- rownames(ord_new)
  gc()
  
  # Get beta distance
  beta_total <- readRDS(beta)$Btotal
  match_values <- match.integer64(as.integer64(rownames(ord_new)), as.integer64(colnames(beta_total)))
  beta_total <- beta_total[match_values, match_values]
  gc()
  
  # Make a frame
  frame <- data.table(PLT_CN = as.integer64(rownames(ord_new)),
                      beta = c(beta_total), 
                      mds = c(lmds_distance),
                      z = 1)
  rm(list = c("beta_total", "lmds_distance"))
  gc()
  
  # Correlation and variance explained
  r <- cor(frame$beta, frame$mds)
  rsquared <- r * r
  
  # Degrees of freedom
  freedom = nrow(frame) - 2
  Fvalue <- rsquared /((1 - rsquared)/freedom)
  Pvalue <- pf(Fvalue, 1, freedom, lower.tail=FALSE)
  
  # Return
  goodness_of_fit <- data.table(r = r,
                                rsquared = rsquared,
                                Fvalue = Fvalue,
                                Pvalue = Pvalue)
  
  fwrite(goodness_of_fit, paste0(out_path, "/goodness_of_fit_", type, "_total.csv"))
  rm(list = c("goodness_of_fit", "env", "vector_match"))
  gc()
  
  #-----------------------------------------------------------------------------
  # Slope estimation
  sma_application <- function(X, plot_id, frame) {
    
    # Subset
    plot <- frame[PLT_CN == plot_id[X], ]
    
    # Test
    test <- sma(mds ~ beta, data = plot)
    elevation <- test$coef[[1]][1,1]
    slope <- test$coef[[1]][2,1]
    
    # Results
    return(data.table(PLT_CN = X,
                      elevation = elevation,
                      slope = slope,
                      rsquar = test$r2[[1]]))
    
  }
  
  sma <- mclapply(X = 1:length(unique.integer64(frame$PLT_CN)), 
                  FUN = sma_application,
                  plot_id = unique.integer64(frame$PLT_CN),
                  frame = frame,
                  mc.cores = 5)
  
  sma <- do.call(rbind.data.frame, sma)
  
  #-----------------------------------------------------------------------------
  # Data for stressplot
  
  vox <- voxels(frame[,-1], edge_length = c(0.001, 0.001, 0.001), threads = threads, obj.voxels = FALSE)
  rm(list = c("frame"))
  gc()
  
  # List to export
  to_export <- list(max = max(vox$N),
                    vox = vox,
                    linear = sma,
                    goodness_of_fit = goodness_of_fit)
  
  saveRDS(to_export, paste0(out_path, "/data_goodness_of_fit_", type, ".rds"))
  
}

#'------------------------------------------------------------------------------
#' @example

# root_path <- "/media/antonio/Work/MMDB"
root_path <- "/n/netscratch/jcavender_lab/Everyone/guzman/MMDB"
threads <- 99

# Taxonomic
# ordination <- paste0(root_path, "/Diversity/NMDS/NMDS_taxo_total.rds")
ordination <- paste0(root_path, "/Diversity/02_LMDS/lmds_taxo_VOL_hellinger.rds")
beta <- paste0(root_path, "/Diversity/01_beta/beta_taxo_VOL_decosand.rds")
fia_environment <- paste0(root_path, "/FIA/03_FIA-env/fia-env.csv")
sp_meta <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
species <- paste0(root_path, "/Diversity/01_beta/master_matched_comm_phylo_VOL.rds")
fia_traits <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
out_path <- paste0(root_path, "/Diversity/02_LMDS")
type <- "taxo"
postprocessing(ordination, beta, fia_environment, sp_meta, species, fia_traits, threads, out_path, type)

# Phylogenetic
ordination <- paste0(root_path, "/Diversity/02_LMDS/lmds_phylo_VOL_hellinger_norm.rds")
beta <- paste0(root_path, "/Diversity/01_beta/beta_phylo_VOL_hellinger_norm.rds")
fia_environment <- paste0(root_path, "/FIA/03_FIA-env/fia-env.csv")
sp_meta <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
species <- paste0(root_path, "/Diversity/01_beta/master_matched_comm_phylo_VOL.rds")
fia_traits <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
out_path <- paste0(root_path, "/Diversity/02_LMDS")
type <- "phylo"
postprocessing(ordination, beta, fia_environment, sp_meta, species, fia_traits, threads, out_path, type)

# Functional
ordination <- paste0(root_path, "/Diversity/02_LMDS/lmds_func_VOL_hellinger_norm.rds")
beta <- paste0(root_path, "/Diversity/01_beta/beta_func_VOL_hellinger_norm.rds")
fia_environment <- paste0(root_path, "/FIA/03_FIA-env/fia-env.csv")
sp_meta <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
species <- paste0(root_path, "/Diversity/01_beta/master_matched_comm_func_VOL.rds")
fia_traits <- paste0(root_path, "/FIA/02_FIA-traits/06_fia-traits.csv")
out_path <- paste0(root_path, "/Diversity/02_LMDS")
type <- "func"
postprocessing(ordination, beta, fia_environment, sp_meta, species, fia_traits, threads, out_path, type)

