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

prepareClusterHeatmap <- function(metabolite
                             #, colour_scale_labels =  c(0,10,100,1000)
                             #, heatmap_colour = "#152736"
){

### Arguments for testing
#metabolite = "3,6-anhydro-D-galactose"
#metabolite = "xylose"
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
  
  ##### Sort a dendrogram for the CAZymes ####
  
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
  
  # If any CAZymes aren't in a cluster, add them to the dataframe 
  # so they are in cluster 0
  if(length(rownames(mat_horiz)[rownames(mat_horiz) %ni% Clusters_horiz$CAZyme]) != 0){
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
  
  #### Set up the phylum dendrogam for ggplot2 ####
  
  # Find the nodes which contain each cluster
  ClusterNodes_Phylum <- unlist(lapply(seq_along(unique(Phylum_Clusters$Cluster[Phylum_Clusters$Cluster != 0]))
                                       , FUN = function(i){
                                         which_node(Phylum_dend
                                                    , Phylum_Clusters$Phylum[Phylum_Clusters$Cluster == i])
                                       }))
  
  dend_data <- dendro_data(Phylum_dend)
  
  # Setup the data, so that the layout is inverted (this is more 
  # "clear" than simply using coord_flip())
  segment_data <- with(
    segment(dend_data), 
    data.frame(x = y, y = x, xend = yend, yend = xend))
  
  # Use the dendrogram label data to position the Phylum labels
  Phylum_pos_table <- with(
    dend_data$labels, 
    data.frame(y_center = x, Phylum = as.character(label), height = 1))
  
  # Add data about which cluster the different phyla were in
  Phylum_pos_table <- merge(Phylum_pos_table, Phylum_Clusters, by = "Phylum")
  
  # Get further information about the dendrogram
  Phylum_dend2 <- as.ggdend(Phylum_dend)
  Phylum_dend2$segments$lwd <- 2
  Phylum_dend2$nodes$text <- ""
  
  Phylum_dend2$nodes$text[ClusterNodes_Phylum] <- 1:length(ClusterNodes_Phylum)
  
  segment_data$col <- Phylum_dend2$segments$col
  
  # ggplot(Phylum_dend2) + 
  #   geom_text(data = Phylum_dend2$nodes[ClusterNodes_Phylum,]
  #             , aes(x = x-0.5, y = y+2
  #                   , label = text)
  #   )
  
  
  #### Set up the CAZyme data for plotting ####
  ClusterNodes_CAZyme <- unlist(lapply(seq_along(unique(CAZyme_Clusters$Cluster[CAZyme_Clusters$Cluster != 0]))
                                       , FUN = function(i){
                                         which_node(CAZyme_dend
                                                    , CAZyme_Clusters$CAZyme[CAZyme_Clusters$Cluster == i])
                                       }))
  
  dend_data_horiz <- dendro_data(CAZyme_dend)
  
  # Setup the data, so that the layout is inverted (this is more 
  # "clear" than simply using coord_flip())
  segment_data_horiz <- with(
    segment(dend_data_horiz), 
    data.frame(x = x, y = y, xend = xend, yend = yend))
  
  # Get further information about the dendrogram
  CAZyme_dend2 <- as.ggdend(CAZyme_dend)
  CAZyme_dend2$segments$lwd <- 2
  CAZyme_dend2$nodes$text <- ""
  CAZyme_dend2$nodes$text[ClusterNodes_CAZyme] <- 1:length(ClusterNodes_CAZyme)
  
  segment_data_horiz$col <- CAZyme_dend2$segments$col
  
  #plot(CAZyme_dend)
  
  # ggplot(CAZyme_dend2) + 
  #   geom_text(data = CAZyme_dend2$nodes[ClusterNodes_CAZyme,]
  #              , aes(x = x-1, y = y+2
  #                    , label = text)
  #             )
  
  #####
  
  # Use the dendrogram label data to position the Phylum labels
  CAZyme_pos_table <- with(
    dend_data_horiz$labels, 
    data.frame(x_center = x, CAZyme = as.character(label), width = 1))
  
  # Add data about which cluster the different phyla were in
  CAZyme_pos_table <- merge(CAZyme_pos_table, CAZyme_Clusters, by = "CAZyme"
  )
  
  # Make the CAZyme matrix ggplot2able
  # add a Phylum column
  mat_2 <- XyloseRichness[2:length(XyloseRichness)] + 1
  mat_2$Phylum <- row.names(mat)
  
  heatmap_data <- mat_2 %>% 
    reshape2::melt(value.name = "Abundance"
                   , varnames = c("Phylum", "CAZyme")) %>%
    dplyr::left_join(Phylum_pos_table) 
  names(heatmap_data)[2] <- 'CAZyme'
  heatmap_data <- heatmap_data %>%
    dplyr::left_join(CAZyme_pos_table, by = "CAZyme")
  
  # Limits for the vertical axes
  Phylum_axis_limits <- with(
    Phylum_pos_table, 
    c(min(y_center - 0.5 * height)
      , max(y_center + 0.5 * height))
  ) + 
    0.1 * c(-1, 1) # extra spacing: 0.1
  
  # Limits for the horizontal axis
  CAZyme_axis_limits <- with(
    CAZyme_pos_table, 
    c(min(x_center - 0.5 * width)
      , max(x_center + 0.5 * width))
  ) + 
    0.1 * c(-1, 1) # extra spacing: 0.1
  
  # Stop errors because of missing values where they aren't needed
  heatmap_data$x_center[is.na(heatmap_data$x_center)] <- 1
  heatmap_data$width[is.na(heatmap_data$width)] <- 1
  
  # The function returns all of the information needed for plotting 
  # the heatmap, but means that the plot can be quickly modified :)
  return(
    list(
      metabolite = metabolite
      , heatmap = list(
        heatmap_data = heatmap_data
        , CAZyme_pos_table = CAZyme_pos_table
        , Phylum_pos_table = Phylum_pos_table
        , Phylum_axis_limits = Phylum_axis_limits)
      , phylum_dend = list(
        segment_data = segment_data
        , Phylum_pos_table = Phylum_pos_table
        , Phylum_axis_limits = Phylum_axis_limits
        , Phylum_dend_node_data = Phylum_dend2$nodes
        , ClusterNodes_Phylum = ClusterNodes_Phylum
      )
      , cazyme_dend = list(
        segment_data_horiz = segment_data_horiz
        , CAZyme_pos_table = CAZyme_pos_table
        , CAZyme_axis_limits = CAZyme_axis_limits
        , CAZyme_dend_node_data = CAZyme_dend2$nodes
        , ClusterNodes_CAZyme = ClusterNodes_CAZyme
      )
    ))
}

prepareClusterHeatmap <- cmpfun(prepareClusterHeatmap)
