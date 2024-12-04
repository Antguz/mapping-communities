#'------------------------------------------------------------------------------
#' @title Codes to visualize the MDS_predicted_mean values
#'------------------------------------------------------------------------------

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(outliers)

#'------------------------------------------------------------------------------
#' Working directory

root_path <- "/media/antonio/Work/MMDB"
sensor <- "DESIS"
folder <- paste0(root_path, "/Models/01_model_beta/test/", sensor)

taxo <- fread(paste0(folder, "/", sensor, "_", "taxo", "_predicted_training.csv"))
taxo$type <- "Taxonomic"
phylo <- fread(paste0(folder, "/", sensor, "_", "phylo", "_predicted_training.csv"))
phylo$type <- "Phylogenetic"
func <- fread(paste0(folder, "/", sensor, "_", "func", "_predicted_training.csv"))
func$type <- "Functional"

data <- rbind(taxo, phylo, func)

predicted_sd <- data[, c("Axes", "sample", "type", "predicted_sd")]
residuals <- data[, c("Axes", "sample", "type", "residuals")]
residuals_sd <- data[, c("Axes", "sample", "type", "residuals_sd")]

predicted_sd <- dcast(predicted_sd, Axes + sample ~ type)
residuals <- dcast(residuals, Axes + sample ~ type)
residuals_sd <- dcast(residuals_sd, Axes + sample ~ type)

top <- 0.99
low <- 0.01

predicted_sd <- predicted_sd[, .(sample = sample,
                                 Taxonomic = scores(Taxonomic, type = "z", prob = TRUE),
                                 Phylogenetic = scores(Phylogenetic, type = "z", prob = TRUE),
                                 Functional = scores(Functional, type = "z", prob = TRUE)), by = "Axes"]

residuals <- residuals[, .(sample = sample,
                           Taxonomic = scores(Taxonomic, type = "z", prob = TRUE),
                           Phylogenetic = scores(Phylogenetic, type = "z", prob = TRUE),
                           Functional = scores(Functional, type = "z", prob = TRUE)), by = "Axes"]

residuals_sd <- residuals_sd[, .(sample = sample,
                                 Taxonomic = scores(Taxonomic, type = "z", prob = TRUE),
                                 Phylogenetic = scores(Phylogenetic, type = "z", prob = TRUE),
                                 Functional = scores(Functional, type = "z", prob = TRUE)), by = "Axes"]

predicted_sd[Taxonomic >= top, Taxonomic := NA]
predicted_sd[Phylogenetic >= top, Phylogenetic := NA]
predicted_sd[Functional >= top, Functional := NA]

residuals[Taxonomic >= top, Taxonomic := NA]
residuals[Taxonomic <= low, Taxonomic := NA]
residuals[Phylogenetic >= top, Phylogenetic := NA]
residuals[Phylogenetic <= low, Phylogenetic := NA]
residuals[Functional >= top, Functional := NA]
residuals[Functional <= low, Functional := NA]

residuals_sd[Taxonomic >= top, Taxonomic := NA]
residuals_sd[Phylogenetic >= top, Phylogenetic := NA]
residuals_sd[Functional >= top, Functional := NA]

predicted_sd[is.na(Taxonomic) | is.na(Phylogenetic) | is.na(Functional), predicted_sd := TRUE]
predicted_sd[is.na(predicted_sd), predicted_sd := FALSE]
residuals[is.na(Taxonomic) | is.na(Phylogenetic) | is.na(Functional), residuals := TRUE]
residuals[is.na(residuals), residuals := FALSE]
residuals_sd[is.na(Taxonomic) | is.na(Phylogenetic) | is.na(Functional), residuals_sd := TRUE]
residuals_sd[is.na(residuals_sd), residuals_sd := FALSE]

predicted_sd <- dcast(predicted_sd[, c("Axes", "sample", "predicted_sd")], sample ~ Axes)
residuals <- dcast(residuals[, c("Axes", "sample", "residuals")], sample ~ Axes)
residuals_sd <- dcast(residuals_sd[, c("Axes", "sample", "residuals_sd")], sample ~ Axes)

predicted_sd[Axis1 | Axis2 | Axis3, predicted_sd := TRUE]
predicted_sd[is.na(predicted_sd), predicted_sd := FALSE]
residuals[Axis1 | Axis2 | Axis3, residuals := TRUE]
residuals[is.na(residuals), residuals := FALSE]
residuals_sd[Axis1 | Axis2 | Axis3, residuals_sd := TRUE]
residuals_sd[is.na(residuals_sd), residuals_sd := FALSE]

plot_to_keep <- cbind(predicted_sd[, c(1, 5)], residuals[, 5], residuals_sd[, 5])
plot_to_keep[predicted_sd == TRUE | residuals == TRUE | residuals_sd == TRUE, remove := TRUE]
plot_to_keep[is.na(remove), remove := FALSE]
nrow(plot_to_keep[remove == TRUE,])/nrow(plot_to_keep)
nrow(plot_to_keep) - nrow(plot_to_keep[remove == TRUE,])
plot_to_keep <- plot_to_keep[, c("sample", "remove")]

# Second

taxo <- fread(paste0(folder, "/", sensor, "_", "taxo", "_residuals_3D_training.csv"))
taxo$type <- "Taxonomic"
phylo <- fread(paste0(folder, "/", sensor, "_", "phylo", "_residuals_3D_training.csv"))
phylo$type <- "Phylogenetic"
func <- fread(paste0(folder, "/", sensor, "_", "func", "_residuals_3D_training.csv"))
func$type <- "Functional"

all <- cbind(taxo= taxo$residuals_3D, phylo= phylo$residuals_3D, func= func$residuals_3D)
mahal_dist <- mahalanobis(all, colMeans(all), cov(all))
threshold <- qchisq(0.99, df = ncol(all))
plot_to_keep[mahal_dist >= threshold, remove := TRUE]
nrow(plot_to_keep[remove == TRUE])

fwrite(plot_to_keep, paste0(root_path, "/Models/01_model_beta/test/", sensor, "_outliers.csv"))
