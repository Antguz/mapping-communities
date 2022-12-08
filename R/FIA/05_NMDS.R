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
library(ggrepel)
#library(BiodiversityR)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param stand
#' @param out_path

root_path <- "F:/DESIS"
root_path <- "/media/antonio/Extreme SSD/DESIS"
stand <- fread(paste0(root_path, "/04_match_spectra-stand/species.csv"))
out_path <- paste0(root_path, "/05_NMDS")

#'------------------------------------------------------------------------------
#' @analysis  

#Metadata
meta <- stand[,1:5]
data <- stand[, .SD, .SDcols=c(6:ncol(stand))]

#Remove plots without species and less than 5 species
nsp <- apply(data, 1, function(x) sum(x > 0))
meta <- meta[nsp != 0]
data <- data[nsp != 0]

#Remove rare species (presence <10 in plots)
presence <- colSums( data != 0)
presence <- names(presence[presence < 5])
remove <- NA

for(i in 1:nrow(data)) {
  
  value <- sum(data[i, ..presence])
  
  if(value > 0) {
    remove[i] <- TRUE
  } else {
    remove[i] <- FALSE
  }
  
}

data <- data[isFALSE(remove) == FALSE]
meta <- meta[isFALSE(remove) == FALSE]
data <- data[, !..presence]

#Transform
data_transfrom <- sqrt(data)
data_dist <- vegdist(data_transfrom, "bray")
PCOA <- pcoa(data_dist, correction = "lingoes", rn = meta$PLT_CN)
saveRDS(PCOA, "PCoA.rds")
PCOA <- readRDS("PCoA.rds")

biplot(PCOA, Y = data_transfrom, plot.axes = c(1,2))

plots_fia <- as.data.table(PCOA$vectors)
species <- wascores(x = PCOA$vectors, w = data_transfrom)
SPCD <- rownames(species)
species <- as.data.table(species)

plots_fia <- cbind(meta, plots_fia)
species <- cbind(data.table(species = SPCD), species)

#Export
fwrite(species, paste0(out_path, "/PCoA_species.csv"))
fwrite(plots_fia, paste0(out_path, "/PCoA_plots.csv"))

species <- fread(paste0(out_path, "/PCoA_species.csv"), header = TRUE)
plots_fia <- fread(paste0(out_path, "/PCoA_plots.csv"), header = TRUE)

plot(plots_fia$Axis.1, plots_fia$Axis.2)
points(species$Axis.20, species$Axis.10, col = "red")

ggplot() +
  geom_point(data = plots_fia, aes(x = Axis.1, y = Axis.2), shape = 18, colour = "grey", alpha = 0.5) +
  stat_density_2d(data = plots_fia, aes(x = Axis.1, y = Axis.2, fill = ..level..), alpha = 0.5, geom = "polygon", contour = TRUE) +
  #stat_density_2d(data = plots_fia, aes(x = Axis.1, y = Axis.2, fill = after_stat(level)), alpha = 0.2, geom = "polygon", contour = TRUE) +
  scale_fill_continuous(low="lavenderblush", high="red") +
  #geom_jitter(alpha=0.5, size = 1.1) +
  geom_point(data = species, aes(Axis.1, Axis.2), colour = "white", fill = "black", shape = 21) +
  theme_bw() +
  scale_x_continuous(limits = c(-.5, .5)) +
  scale_y_continuous(limits = c(-.5, .5)) +
  geom_text_repel(data = species, aes(x = Axis.1, y = Axis.2),
    label= species$species, 
    direction = "both", 
    box.padding = 0.5) + xlab("PCoA1") + ylab("PCoA2") +
  theme(legend.position="none")

ggplot() +
  geom_point(data = plots_fia, aes(x = Axis.1, y = Axis.3), shape = 18, colour = "grey", alpha = 0.5) +
  stat_density_2d(data = plots_fia, aes(x = Axis.1, y = Axis.3, fill = ..level..), alpha = 0.5, geom = "polygon", contour = TRUE) +
  #stat_density_2d(data = plots_fia, aes(x = Axis.1, y = Axis.2, fill = after_stat(level)), alpha = 0.2, geom = "polygon", contour = TRUE) +
  scale_fill_continuous(low="lavenderblush", high="red") +
  #geom_jitter(alpha=0.5, size = 1.1) +
  geom_point(data = species, aes(Axis.1, Axis.3), colour = "white", fill = "black", shape = 21) +
  theme_bw() +
  scale_x_continuous(limits = c(-.5, .5)) +
  scale_y_continuous(limits = c(-.5, .5)) +
  geom_text_repel(data = species, aes(x = Axis.1, y = Axis.3),
                  label= species$species, 
                  direction = "both", 
                  box.padding = 0.5) + xlab("PCoA1") + ylab("PCoA3") +
  theme(legend.position="none")

ggplot() +
  geom_point(data = plots_fia, aes(x = Axis.2, y = Axis.3), shape = 18, colour = "grey", alpha = 0.5) +
  stat_density_2d(data = plots_fia, aes(x = Axis.2, y = Axis.3, fill = ..level..), alpha = 0.5, geom = "polygon", contour = TRUE) +
  #stat_density_2d(data = plots_fia, aes(x = Axis.1, y = Axis.2, fill = after_stat(level)), alpha = 0.2, geom = "polygon", contour = TRUE) +
  scale_fill_continuous(low="lavenderblush", high="red") +
  #geom_jitter(alpha=0.5, size = 1.1) +
  geom_point(data = species, aes(Axis.2, Axis.3), colour = "white", fill = "black", shape = 21) +
  theme_bw() +
  scale_x_continuous(limits = c(-.5, .5)) +
  scale_y_continuous(limits = c(-.5, .5)) +
  geom_text_repel(data = species, aes(x = Axis.2, y = Axis.3),
                  label= species$species, 
                  direction = "both", 
                  box.padding = 0.5) + xlab("PCoA2") + ylab("PCoA3") +
  theme(legend.position="none")







grasslands.spe <- read.delim ('http://www.davidzeleny.net/anadat-r/data-download/grasslands-spe.txt', row.names = 1)
grasslands.env <- read.delim ('http://www.davidzeleny.net/anadat-r/data-download/grasslands-env.txt')
grasslands.spe.log <- log1p (grasslands.spe)

PCA <- rda (grasslands.spe.log)

vegtype <- as.numeric (as.factor (grasslands.env$classification))
ordiplot (PCA, display = 'sites', type = 'n')
points (PCA, pch = vegtype, col = vegtype)






#'------------------------------------------------------------------------------
# Test the significance of stress value for observed data

# Repetitions and dimensions
reps <- 1:10
k_dimensions <- 1:10
l <- seq(100, 10000, by = 250)

frame <- CJ(reps, k_dimensions, l)
#frame <- frame[1:26]

#Progress bar
###Progress bar -----------------------------------
pb <- txtProgressBar(min = 1, max = nrow(frame), style = 3) 
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

###Parallel---------------------------------------
cl <- makeCluster(25)
doSNOW::registerDoSNOW(cl)

data <- as.matrix(data_log)

repretition <- foreach(i = 1:nrow(frame), .combine = rbind, .packages = c("bigmds", "ecodist"), .options.snow = opts) %dopar% {
  
  result <- frame[i, ]
  
  mds <- try(fast_mds(x = data, 
                  l = 10000, 
                  s_points = 2*3, 
                  r = 3, 
                  n_cores = 12, 
                  dist_fn = ecodist::bcdist), silent = TRUE)
  
  if(class(mds) == "try-error") {
    
    result$GOF1 <- NA
    result$GOF2 <- NA
    
  } else {
    
    result$GOF1 <- mds$GOF[1]
    result$GOF2 <- mds$GOF[2]
    
  }
  
  return(result)
  
}

close(pb)
stopCluster(cl)


#Get evaluate different k 
for(i in 1:length(k_dimensions)) {
  
  setTxtProgressBar(pb, i)
  
  #stress <- oecosimu(comm = data,  method = "quasiswap_count", 
  #                   nestfun = metaMDS, autotransform = FALSE, k = k_dimensions[i], 
  #                   distance = "bray", nsimul = reps, parallel = 26, statistic = "stress", 
  #                   alternative = "less", trace = TRUE, maxit = 1000, 
  #                   trymax = 200, sratmax = 0.9999999)
  
  #clus <- makeCluster(26, type = "FORK") #, 
  #clusterEvalQ(clus, library(vegan))
  

  
  #stopCluster(clus)
  
  results[[i]] <- NMDS
  names(results)[i] <- paste0("k_", k_dimensions[i])
  
  Sys.sleep(5)
  
}






library(ape)
library(vegan)

#data(mite)
mite.log <- log(data+1)
mite.D <- vegdist(mite.log, "bray")

mite.log.st = apply(mite.log, 2, scale, center=TRUE, scale=TRUE)

par(mfrow=c(1,2))
biplot(res, mite.log)
biplot(res, mite.log.st)


data(varespec)
dist <- vegdist(varespec,  method = "bray")

# PCoA is not included in vegan. 
# We will use the ape package instead
library(ape)
PCOA <- pcoa(dist)
biplot.pcoa(PCOA)

# plot the eigenvalues and interpret
barplot(PCOA$values$Relative_eig[1:10])
# Can you also calculate the cumulative explained variance of the first 3 axes?

# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOA <- pcoa(dist, correction = "cailliez")

# Plot your results
biplot.pcoa(PCOA)

# You see what`s missing? 
# Indeed, there are no species plotted on this biplot. 
# That's because we used a dissimilarity matrix (sites x sites) 
# as input for the PCOA function. 
# Hence, no species scores could be calculated. 
#However, we could work around this problem like this:
biplot.pcoa(PCOA, varespec)
PCOAaxes <- PCOA$vectors[,c(1,2)]
par(mfrow = c(1, 2)) 
biplot.pcoa(PCOA)
plot(PCA)

PCOA <- readRDS("data/PCoA.rds")
plot(PCOA$values$Relative_eig[1:45])

wascores(PCOA$vectors, )



#Perform ordination
NMDS <- metaMDS(data, 
                distance = "bray", 
                k = 3, 
                parallel = 26, 
                autotransform = FALSE,
                try = 20, 
                trymax = 100)

saveRDS(NMDS, file = "data/NMDS2.rds")


#########################
#Data visualisation (THIS IS AN UPDATED VERSION OF THE SCRIPT, NOW USING GGPLOT)
library("ggplot2")

NMDS <- readRDS("data/NMDS.rds")

#Extract the axes scores
sites <- as.data.frame(scores(NMDS)$sites)  #extract the site scores
species <- as.data.frame(scores(NMDS)$species)

#plot
ggplot() +
  geom_point(data = sites, aes(x = NMDS1, y = NMDS2)) +
  coord_fixed(xlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1))+                                              
  theme_bw()+ 
  theme(legend.position="right",legend.text=element_text(size=10),legend.direction='vertical')




#####################
#Bootstrapping and testing for differences between the groups
fit <- adonis(data_1 ~ Habitat, data=data_2, permutations=999, method="bray")
fit


#####################
#Check assumption of homogeneity of multivariate dispersion
distances_data <- vegdist(data_1)
anova(betadisper(distances_data, data_2$Habitat))





ggplot() +
  geom_point(data = plots_fia, aes(x = Axis.1, y = Axis.2), shape = 18, colour = "grey", alpha = 0.5) 








