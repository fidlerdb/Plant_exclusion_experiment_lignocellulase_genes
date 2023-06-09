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
library(analogue)
library(clustsig)
library(compiler)
library(gplots)
library(Heatplus)
library(dendextend)
library(ggdendro)

# Chi Square distance metric replacing fully 0 rows with distances of 0
Chi_dist <- function(x){d <- distance(x, method = "chi.square"
                                      , dist = TRUE)
d[is.na(d)] <- 0
return(d)
}

#### prepareClusterHeatmap function ####

prepareClusterHeatmap <- function(inputMatrix, simprof_expected = 10000
                                  , simprof_simulated = 10000
                                  , distance_metric = 'euclidean'
                                  #, colour_scale_labels =  c(0,10,100,1000)
                                  #, heatmap_colour = "#152736"
){
  
  ### Arguments for testing
  #metabolite = "3,6-anhydro-D-galactose"
  # metabolite = "4-hydroxybenzoic acid"
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
  
  # Create a CAZyme abundance matrix
  mat <- inputMatrix
  
  # Remove CAZymes with no hits
  SUMS <- colSums(mat)
  ZeroSumColumns <- names(SUMS[SUMS == 0])
  mat <- mat[, !names(mat) %in% ZeroSumColumns]
  
  # Add 1 so log transformation works
  #mat <- mat + 1 
  
  # Format the matrix
  #rownames(mat) <- XyloseRichness$Phylum
  #CAZyme_names <- colnames(mat)
  
  # Make the CAZyme matrix ggplot2able
  # add a Phylum column -- this will become the heatmap_data object
  mat_2 <- mat
  #mat_2$Phylum <- row.names(mat)
  
  #### Significance tests for clusters of phyla ####
  
  # SIMPROF analysis on phyla--which phyla have similar collections of genes?
  # Chi square distance
  simprofResults <- simprof(data = mat
                            , method.cluster = 'ward.D2'
                            , method.distance = distance_metric#Chi_dist
                            , num.expected = simprof_expected
                            , num.simulated = simprof_simulated
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
  
  ##### Sort a dendrogram for the CAZymes, or not if custom column order is provided ####
  
  # Work with the transposed data
  
  mat_horiz <- t(mat)
  
  # SIMPROF analysis on CAZy families--which families are present 
  # in the same phyla?
  simprofResults_horiz <- simprof(data = mat_horiz#-1
                                  , method.cluster = 'ward.D2'
                                  , method.distance = distance_metric#Chi_dist
                                  , num.expected = simprof_expected
                                  , num.simulated = simprof_simulated)
  
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
  
  
  # Create the final heatmap_data object
  mat_2$Phylum <- rownames(mat_2)
  mat_2 <- data.table(mat_2)
  
  heatmap_data <- mat_2 %>% 
    data.table::melt(id = "Phylum"
                     , value.name = "Abundance"
                     , varnames = c("Phylum", "CAZyme", "Abundance")) %>%
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
      #metabolite = metabolite
      heatmap = list(
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

############ plotClusterHeatmap function #############

plotClusterHeatmap <- function(df, x_axis_size = 8, y_axis_size = 8
                               , dendrogram_left_margin = -0.6, x_axis_vjust = 0.5, x_axis_hjust = 0.5
                               , colour_labels =  c(0,25,50, 75, 100), heatmap_colour = "#152736"
                               , colour_scale_labels = c(0, 25, 50, 75, 100), colour_scale_title = "Number of Genes"
                               , side_dend_width = 0.5
                               , top_dend_height = 0.5, side_dendrogram_margin = c(0, 0.2, 0.2, dendrogram_left_margin)
                               , heatmap_margin = c(0.2, -0.7, 0.2, 0.2)
                               , heatmap_ylab = "Phylum"
                               , column_order = NULL){

library(ggrepel)

### Arguments for testing
# x_axis_size = 8
# y_axis_size = 8
# strip_text_size = 6
# dendrogram_left_margin = -0.7
# x_axis_vjust = 0.5
# colour_labels =  c(0,25,50, 75, 100)
# heatmap_colour = "#152736"
# colour_scale_labels = c(0, 25, 50, 75, 100)
# side_dend_width = 0.2
# top_dend_height = 0.2

#df <- clust_plot_object

#### Heatmap plot ####

if(!is.null(column_order)){
  column_order <- data.frame(CAZyme = column_order, column_order = 1:length(column_order))
  # A quick check
  #dplyr::left_join(df$heatmap$heatmap_data[,c("CAZyme", "x_center")], column_order, by = "CAZyme")$CAZyme == df$heatmap$heatmap_data[,c("CAZyme")]
    
  # Join the two datasets
    
  df$heatmap$heatmap_data <- dplyr::left_join(df$heatmap$heatmap_data, column_order, by = "CAZyme")
  df$heatmap$heatmap_data$x_center <- df$heatmap$heatmap_data$column_order
  
  # Now work on CAZyme_pos_table as well
  
  df$heatmap$CAZyme_pos_table <- dplyr::left_join(df$heatmap$CAZyme_pos_table, column_order, by = "CAZyme")
  df$heatmap$CAZyme_pos_table$x_center <- df$heatmap$CAZyme_pos_table$column_order
    
}
  
plt_hmap <- with(df$heatmap, (ggplot(heatmap_data, 
                                     aes(x = factor(x_center), y = y_center
                                         , fill = Abundance, 
                                         height = height, width = width)) + 
                                geom_tile() +
                                scale_fill_gradient(colour_scale_title
                                                    , high = heatmap_colour, low = "white"
                                                    #, trans = 'log10'
                                                    , breaks = colour_labels
                                                    , labels = colour_labels) +
                                scale_x_discrete(breaks = factor(CAZyme_pos_table$x_center), 
                                                 labels = CAZyme_pos_table$CAZyme, 
                                                 expand = c(0, 0)) + 
                                # For the y axis, alternatively set the labels as: Phylum_position_table$Phylum
                                scale_y_continuous(breaks = Phylum_pos_table[, "y_center"], 
                                                   labels = Phylum_pos_table$Phylum,
                                                   limits = Phylum_axis_limits, 
                                                   expand = c(0, 0)) + 
                                labs(x = "CAZyme Cluster"
                                     , y = if(heatmap_ylab == "NULL"){element_blank()} else {heatmap_ylab}
                                     ) +
                                theme_bw() +
                                theme(axis.text.x = element_text(size = x_axis_size
                                                                 #, hjust = x_axis_hjust
                                                                 , angle = 90
                                                                 , vjust = x_axis_vjust
                                                                 , hjust = x_axis_hjust),
                                      axis.text.y = element_text(size = y_axis_size
                                                                 , colour = Phylum_pos_table$Colour),
                                      # margin: top, right, bottom, and left
                                      plot.margin = unit(heatmap_margin, "cm"), 
                                      panel.grid = element_blank(),
                                      legend.position = 'bottom',
                                      panel.spacing = unit(0, "lines")
                                ))
)

#### Dendrogram plot ####
plt_dendr <- with(df$phylum_dend, ggplot(segment_data) + 
                    geom_segment(aes(x = x, y = y
                                     , xend = xend, yend = yend
                                     , col = col)) + 
                    #scale_x_reverse(expand = c(0, 0.5)) + 
                    scale_x_continuous(expand = c(0, 0.5)) +
                    scale_y_continuous(breaks = Phylum_pos_table$y_center, 
                                       limits = Phylum_axis_limits, 
                                       expand = c(0, 0)) + 
                    labs(x = "", y = "", colour = "", size = "") +
                    theme_bw() + 
                    theme(panel.grid = element_blank()
                          , axis.text = element_blank()
                          , rect = element_blank()
                          , axis.ticks = element_blank()
                          # trbl
                          , plot.margin = unit(side_dendrogram_margin, "cm")
                          , legend.position = "none"
                    ) #+
                  # geom_text(data = Phylum_dend_node_data[ClusterNodes_Phylum,]
                  #           , aes(x = y + 1, y = x + 1
                  #                 , label = text
                  #                 , size = 6
                  #                 , fontface = "bold"
                  #                 #, col = na.omit(unique(segment_data$col))
                  #           ), colour = 'grey50'
                  # )
                  #geom_text_repel(data = Phylum_dend_node_data[ClusterNodes_Phylum,]
                  #                , aes(x = y , y = x
                  #                      , label = text
                  #                      , size = 6
                  #                      , fontface = "bold"
                  #                ), colour = 'grey50'
                  #                , nudge_x = 1
                  #, segment.size  = 0.2
                  #, segment.color = "grey50"
                  #, direction     = "x"
                  #)
)

#### Top dendrogram  ####
plt_dendr_horiz <- with(df$cazyme_dend, 
                        ggplot(segment_data_horiz) + 
                          geom_segment(aes(x = x 
                                           , y = y
                                           , xend = xend 
                                           , yend = yend
                                           , col = col
                                           )) + 
                          scale_x_continuous(breaks = CAZyme_pos_table$x_center, 
                                             limits = CAZyme_axis_limits, #c(0.5, 50),
                                             expand = c(0, 0)) + 
                          labs(x = "", y = "", colour = "", size = "") +
                          theme_bw() + 
                          theme(panel.grid = element_blank()
                                , axis.text = element_blank()
                                , rect = element_blank()
                                , axis.ticks = element_blank()
                                # trbl
                                , plot.margin = unit(c(0.2, 0, -0.7, 0.0), "cm")
                                , legend.position = "none"
                          ) +
                          #ggtitle(CapStr(df$metabolite)) +
                          geom_text_repel(data = CAZyme_dend_node_data[ClusterNodes_CAZyme,]
                                          , aes(x = x , y = y
                                                , label = text
                                                , size = 6
                                                , fontface = "bold"
                                                
                                                
                                                #, col = na.omit(unique(segment_data_horiz$col))
                                          ), colour = 'grey50'
                                          #, nudge_y = 1
                                          #, segment.size  = 0.2
                                          #, segment.color = "grey50"
                                          #, direction     = "x"
                          )
)

#### Blank Plot ####
BLANK <- ggplot() + cowplot::theme_nothing() +
  #trbl
  theme(plot.margin = unit(c(0.2, 0.2, -0.7, -0.7)
                           , "cm"))

##### Put them all together ####

# metabolitePlot <- plot_grid(plt_dendr+scale_x_reverse()
#                             , plt_hmap
#                             , align = "h"
#                             , axis = "tb"
#                             , rel_widths = c(side_dend_width, 1))


metabolitePlot <- plot_grid(plt_hmap
                            , plt_dendr
                            , align = "h"
                            , axis = "tb"
                            , rel_widths = c(1, side_dend_width))
 return(list(metabolitePlot = metabolitePlot
             , heatmap = plt_hmap
             , dendrogram = plt_dendr))
}
