#'------------------------------------------------------------------------------
#' @title Non-metric Multidimensional Scaling
#'------------------------------------------------------------------------------

#' @description It creates a new ordination of species abundance of plots
#' @return .txt file with the stack of results

#'------------------------------------------------------------------------------
#' Library

library(data.table)
library(vegan)
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
data <- stand[, c(6:72)]
data <- round(data, 2)

#Perform ordination
NMDS <- metaMDS(data, distance = "bray", k = 3, parallel = 6)



#########################
#Data visualisation (THIS IS AN UPDATED VERSION OF THE SCRIPT, NOW USING GGPLOT)
install.packages("ggplot2")
library("ggplot2")

#Extract the axes scores
datascores <- as.data.frame(scores(NMDS))  #extract the site scores

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

