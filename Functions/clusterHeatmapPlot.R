# Split at either " " or "-", then capitalise the first letter of each word
CapStr <- function(y) {
  c <-strsplit(y, "(?=[ -])", perl = TRUE)[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep ="", collapse="")
}

# %ni% function
'%ni%' <- function(x,y)!('%in%'(x,y))

# Load required libraries
library(data.table)
#library(magrittr)
#library(ggplot2)
#library(ggdendro)
#library(egg)
library(analogue)
library(clustsig)
library(compiler)
library(gplots)
library(Heatplus)
library(dendextend)

# Chi Square distance metric replacing fully 0 rows with distances of 0
Chi_dist <- function(x){d <- distance(x, method = "chi.square"
                                      , dist = TRUE)
d[is.na(d)] <- 0
return(d)
}

# customHeatmapPlot function
# Uses Chi^2 distance metric and Ward's D for the clustering

customHeatmapPlot <- function(metabolite
                              , colour_scale_labels =  c(0,10,100,1000)
                              , heatmap_colour = "#152736"
){

### Arguments for testing
#metabolite = "3,6-anhydro-D-galactose"
#metabolite = "vanillic acid"
#x_axis_size = 4
#y_axis_size = 4
#strip_text_size = 6
#dendrogram_left_margin = -0.7
#x_axis_vjust = 0.5
#colour_labels =  c(0,10,100,400)
#heatmap_colour = "#152736"
#colour_scale_labels = c(0, 10, 100, 1000)
#side_dend_width = 0.2
#top_dend_height = 0.2

##### Get the number of genes in each of the metabolite-correlated CAZy families ####
XyloseCAZyRichness <- lapply(1:nrow(cl[cl$Metabolite == metabolite,])
                             , FUN = function(i){
                               table(df$Phylum[df$CAZyme == cl[cl$Metabolite == metabolite,]$CAZymes[i]])
                             })
names(XyloseCAZyRichness) <- cl[cl$Metabolite == metabolite,]$CAZymes

# Create a matrix of CAZymes vs Orders for xylose
xCAZm <- data.frame(matrix(NA
                           , nrow = length(AllPhyla)
                           , ncol = length(cl[cl$Metabolite == metabolite,]$CAZymes)
))
names(xCAZm) <- cl[cl$Metabolite == metabolite,]$CAZymes
rownames(xCAZm) <- AllPhyla

# Put the data about cazyme-containing contig abundance into a matrix 
OCzR <- lapply(seq_along(XyloseCAZyRichness), FUN = function(x){
  
  if (length(XyloseCAZyRichness[[x]]) == 0) {
    
    Rch <- data.frame(Phylum = rownames(xCAZm)
                      , rep(0, nrow(xCAZm)))
    names(Rch)[2] <- names(XyloseCAZyRichness)[x]
    return(Rch)
    
  } else {
    
    # Reference data frame to get all phyla in the dataset
    PhylumRef <- data.frame(Phylum = rownames(xCAZm))
    
    # Richness data for the CAZyme in question
    RichnessData <- data.frame(XyloseCAZyRichness[[x]])
    
    # Merge the two datasets, keeping all the phyla
    Rch <- merge(PhylumRef, RichnessData
                 , by.x = "Phylum", by.y = "Var1"
                 , all.x = TRUE)
    
    # Replace any NA values with 0
    Rch[2][is.na(Rch[2])] <- 0
    
    # Rename the colmns nicely
    names(Rch) <- c("Phylum", names(XyloseCAZyRichness)[x])
    return(Rch)
  }
})

# Stick them side by side, ensuring that values of Phylum match
XyloseRichness <- purrr::reduce(OCzR, dplyr::left_join, by = 'Phylum')

#### Set up data for plotting a heatmap with dendogram ####

# Create a CAZyme abundance matrix
mat <- XyloseRichness[2:length(XyloseRichness)]

# Remove CAZymes with no hits
SUMS <- colSums(mat)
ZeroSumColumns <- names(SUMS[SUMS == 0])
mat <- mat[, !names(mat) %in% ZeroSumColumns]

# Add 1 so log transformation works
mat <- mat + 1 

# Format the matrix
rownames(mat) <- XyloseRichness$Phylum
CAZyme_names <- colnames(mat)

#### Significance tests for clusters of phyla ####

# SIMPROF analysis on phyla--which phyla have similar collections of genes?
# Chi square distance
simprofResults <- simprof(data = mat-1
                          , method.cluster = 'ward.D2'
                          , method.distance = Chi_dist
                          , num.expected = 10000
                          , num.simulated = 10000
                          )

# Obtain the phylum dendrogram
dend <- as.dendrogram(simprofResults$hclust)
#simprof.plot(simprofResults)
#plot(dend)

# Create a data frame identifying which cluster each phylum is in
Clusters <- rbindlist(lapply(seq_along(simprofResults$significantclusters)
                             , FUN = function(x){
                               data.frame(Cluster = x
                                          , Phylum = simprofResults$significantclusters[[x]]
                                          )
                               }))

# If any phyla aren't in a cluster, add them to the dataframe 
# so they are in cluster 0
if(length(rownames(mat)[rownames(mat) %ni% Clusters$Phylum]) != 0){
  Clusters <- rbind(Clusters,
                    data.frame(Cluster = 0
                               , Phylum =  rownames(mat)[rownames(mat) %ni% Clusters$Phylum]
                    ))
}

##### Sort a dendrogram for each CAZy type: ####

# Work with the transposed data
mat_horiz <- t(mat)

# SIMPROF analysis on CAZy families--which families are present 
# in the same phyla?
simprofResults_horiz <- simprof(data = mat_horiz-1
                                , method.cluster = 'ward.D2'
                                , method.distance = Chi_dist
                                , num.expected = 10000
                                , num.simulated = 10000)

# Obtain the dendrogram
dend_horiz <- as.dendrogram(simprofResults_horiz$hclust)
#simprof.plot(simprofResults_horiz)
#plot(dend_horiz)

# Create a data frame identifying which cluster each phylum is in
Clusters_horiz <- rbindlist(lapply(seq_along(simprofResults_horiz$significantclusters), FUN = function(x){
  data.frame(Cluster = x,
             CAZyme = simprofResults_horiz$significantclusters[[x]])
}))

# If any phyla aren't in a cluster, add them to the dataframe 
# so they are in cluster 0
if(length(rownames(mat)[rownames(mat) %ni% Clusters$Phylum]) != 0){
  Clusters_horiz <- rbind(Clusters_horiz,
                        data.frame(Cluster = 0
                                   , CAZyme =  rownames(mat_horiz)[rownames(mat_horiz) %ni% Clusters_horiz$CAZyme]
                        ))}

#### Edit the Phylum dendrogam as necessary ####

# Create a vector of clusters and CAZymes in the order 
# they need to be in for the plot
Phylum_Clusters <- Clusters[match(labels(dend)
                                  , as.character(Clusters$Phylum)),]
# Colour the branches
Phylum_dend <- color_branches(dend, clusters = Phylum_Clusters$Cluster)

# Make the branches all the same width
Phylum_dend <- set(Phylum_dend, "branches_lwd", 3)

#### Edit the  CAZyme dendrogam as necessary ####

# Create a vector of clusters and CAZymes in the order 
# they need to be in for the plot
CAZyme_Clusters <- Clusters_horiz[match(labels(dend_horiz)
                                  , as.character(Clusters_horiz$CAZyme)
),]

# Colour the branches
CAZyme_dend <- color_branches(dend_horiz, clusters = CAZyme_Clusters$Cluster)

# Make the branches all the same width
CAZyme_dend <- set(CAZyme_dend, "branches_lwd", 3)

#### Node labels ####

# Find the highest node containing the cluster
# ClusterNodes <- unlist(lapply(seq_along(unique(CAZyme_Clusters$Cluster))
#                        , FUN = function(i){
#                          which_node(CAZyme_dend
#                                     , CAZyme_Clusters$CAZyme[CAZyme_Clusters$Cluster == i])
#                          }))
#        
# # Get the xy coordinates of every node
# xy <- as.data.frame(CAZyme_dend %>% get_nodes_xy())
# # Note what each cluster is
# xy[,3] <- NA
# xy[ClusterNodes,3] <- 1:length(ClusterNodes)
# names(xy) <- c("x","y","cluster")
# 
# # Set the node attributes
# CAZyme_dend <- set(CAZyme_dend, "nodes_pch", as.character(xy$cluster))
# CAZyme_dend <- set(CAZyme_dend, "nodes_cex", 1.5)

#### The final Plot ####

Palette <- colorRampPalette(c("white", heatmap_colour), space = "rgb")(500)

# Desired actual values for colour labels


TickFun <- function() {
  #breaks <- breaksparent.frame()$breaks
  breaks <- log10(colour_scale_labels+1)
  
  return(list(
    at = parent.frame()$scale01(breaks)
    , labels = as.character(colour_scale_labels)#as.character((10^breaks)-1)
    
  ))
}

heatmap.2(as.matrix(log10(mat))
          , Rowv = Phylum_dend
          , Colv = CAZyme_dend
          
          #, RowSideColors = as.character(Clusters$Cluster)
          #, ColSideColors = as.character(Clusters_horiz$Cluster)
          
          , col = Palette
          , trace = "none"  
          , main = CapStr(metabolite)
          , density.info = "none"
          , key.title = NA
          , key.xlab = "Number of Genes"
          , key.xtickfun = TickFun
          )
}


#plot(heatMapObj);par(mfrow=c(1,1))
