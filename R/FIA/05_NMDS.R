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
#library(BiodiversityR)
options(scipen = 100)

#'------------------------------------------------------------------------------
#' Arguments
#' @param stand
#' @param out_path

root_path <- "E:/DESIS/04_match_spectra-stand"
root_path <- "/media/antonio/Extreme SSD/DESIS/04_match_spectra-stand"
stand <- fread(paste0(root_path, "/species.csv"))
out_path <- paste0(root_path, "/04_match_spectra-stand")

#'------------------------------------------------------------------------------
#' @analysis  

#Metadata
meta <- stand[,1:5]

#Data for ordination
data <- stand[, .SD, .SDcols=c(6:ncol(stand))]
a <- rowSums(data)
data <- setnafill(data, fill = 0) #In case    
data <- data[a != 0]

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

data <- data[remove == FALSE]
data <- data[, !..presence]

#Transform
data_log <- log(data+1)

NMDS <- metaMDS(data_log, 
                distance = "bray", 
                k = 3, 
                autotransform = FALSE,
                try = 5, 
                trymax = 5,
                parallel = 8)

saveRDS(NMDS, "data/NMDS.rds")


NMDS_new <- metaMDS(test, 
                    distance = "bray", 
                    k = 3, 
                    autotransform = FALSE,
                    try = 5, 
                    trymax = 5,
                    previous.best = NMDS,
                    parallel = 16)











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

#Extract the axes scores
datascores <- as.data.frame(scores(NMDS)$sites)  #extract the site scores

#Add/calculate spider diagram
scores <- cbind(as.data.frame(datascores), Habitat = data_2$Habitat)
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Habitat, data = scores, FUN = mean)
seg <- merge(scores, setNames(centroids, c('Habitat','oNMDS1','oNMDS2')),
             by = 'Habitat', sort = FALSE)

#plot
ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = Habitat)) +
  geom_segment(data = seg,
               mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # add spiders
  geom_point(data = centroids, size = 4) +                    # add centroids
  geom_point() +                                              
  coord_fixed()+                                              
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













