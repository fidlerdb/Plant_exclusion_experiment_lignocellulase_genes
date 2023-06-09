
# Split at either " " or "-", then capitalise the first letter of each word
CapStr <- function(y) {
  c <-strsplit(y, "(?=[ -])", perl = TRUE)[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep ="", collapse="")
}

# %ni% function
'%ni%' <- function(x,y)!('%in%'(x,y))

library(data.table)
library(magrittr)
library(ggplot2)
library(ggdendro)
library(egg)

# customHeatmapPlot function
# Uses Chi^2 distance metric and Ward's D for the clustering

customHeatmapPlot <- function(metabolite
                              , x_axis_size = 4
                              , y_axis_size = 4
                              , strip_text_size = 6
                              , dendrogram_left_margin = -0.7
                              , x_axis_vjust = 0.5
                              , colour_labels =  c(0,10,100,400)
                              , heatmap_colour = "#152736"
                              , side_dend_width = 0.2
                              , top_dend_height = 0.2
){
  
  
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
  
  #mat <- log10(mat + 1)
  SUMS <- colSums(mat)
  ZeroSumColumns <- names(SUMS[SUMS == 0])
  mat <- mat[, !names(mat) %in% ZeroSumColumns]
  mat <- mat + 1
  
  rownames(mat) <- XyloseRichness$Phylum
  CAZyme_names <- colnames(mat)
  
  # Obtain the dendrogram
  library(analogue)
  d <- distance(mat-1, method = "chi.square", dist = TRUE)
  d[is.na(d)] <- 0
  
  dend <- as.dendrogram(hclust(d, method = "ward.D2"))
  dend_data <- dendro_data(dend)
  
  ##### Sort a dendrogram for each CAZy type: ####
  mat_horiz <- t(mat)
  
  # Seperate the CAZyme types
  CT <- gsub("endo-1,4-β-xylanase", "", rownames(mat_horiz))
  CT <- gsub("[0-9]*", "", CT)
  CT <- gsub("β-Glucosidase", "", CT)
  CT <- gsub("_", "", CT)
  CT[CT == "CBM.CBM.GH"] <- "CBM.GH"
  CT[CT == "CBM.GH.GH"] <- "CBM.GH"
  
  # Get the data for each cazy type
  dend_data_horiz <- by(mat_horiz, CT, function(x){
    if(nrow(x) == 1){
      return(NA)
    } else {
      d <- distance(x-1, method = "chi.square", dist = TRUE)
      d[is.na(d)] <- 0
      
      dend <- as.dendrogram(hclust(d, method = "ward.D2"))
      dend_data <- dendro_data(dend)
      
      #dend <- as.dendrogram(hclust(dist(x)))
      #dend_data <- dendro_data(dend)
      return(dend_data)}
  })
  
  #### Set up the phylum dendrogam ####
  
  # Setup the data, so that the layout is inverted (this is more 
  # "clear" than simply using coord_flip())
  segment_data <- with(
    segment(dend_data), 
    data.frame(x = y, y = x, xend = yend, yend = xend))
  
  # Use the dendrogram label data to position the Phylum labels
  Phylum_pos_table <- with(
    dend_data$labels, 
    data.frame(y_center = x, Phylum = as.character(label), height = 1))
  
  #### Set up the CAZyme data for plotting ####
  
  nFamilies <- unlist(lapply(dend_data_horiz, FUN = function(x){
    if(is.na(x)){return(1)}
    max(segment(x)[,1])
  }))
  
  toAdd_x <- cumsum(nFamilies)
  toAdd_x <- data.frame(toAdd_x = c(0, toAdd_x[1:length(toAdd_x) -1])
                        , CAZy_Type =  names(toAdd_x))
  
  dend_data_horiz_2 <- lapply(dend_data_horiz, FUN = function(x){
    if(is.na(x)){return()}
    merge(segment(x), label(x)[c(1,3)]
          , by.x = 'xend', by.y = 'x'
          , all.x = TRUE)}
  )
  
  segment_data_horiz <- rbindlist(dend_data_horiz_2
                                  , idcol = "CAZy_Type")
  
  segment_data_horiz <- merge(segment_data_horiz, toAdd_x
                              , all.x = TRUE)
  
  # Use the dendrogram label data to position the Phylum labels
  
  CAZyme_pos_table <- rbindlist(lapply(dend_data_horiz, FUN = function(x){
    if(is.na(x)){return(
      data.frame(x_center = 1
                 , CAZyme = NA#as.character(names(x))
                 , width = 1)
    )
    }
    with(x$labels,
         data.frame(x_center = x
                    , CAZyme = as.character(label)
                    , width = 1))
  }))
  
  NoDendrogram <- names(mat)[names(mat) %ni% c("phylum", as.character(CAZyme_pos_table$CAZyme))] # need to add these
  CAZyme_pos_table[is.na(CAZyme_pos_table$CAZyme),]$CAZyme <- NoDendrogram
  
  # Add a CAZy_Type column to classify the domains
  CAZyme_pos_table$CAZy_Type <- gsub("endo-1,4-β-xylanase", "", CAZyme_pos_table$CAZyme)
  CAZyme_pos_table$CAZy_Type <- gsub("[0-9]*", "", CAZyme_pos_table$CAZy_Type)
  CAZyme_pos_table$CAZy_Type <- gsub("β-Glucosidase", "", CAZyme_pos_table$CAZy_Type)
  CAZyme_pos_table$CAZy_Type <- gsub("_", "", CAZyme_pos_table$CAZy_Type)
  
  # Make the ones which have repeat domain types in lump with the others 
  CAZyme_pos_table$CAZy_Type[CAZyme_pos_table$CAZy_Type == "CBM.CBM.GH"] <- "CBM.GH"
  CAZyme_pos_table$CAZy_Type[CAZyme_pos_table$CAZy_Type == "CBM.GH.GH"] <- "CBM.GH"
  
  # Add data aboout which facet they will be in
  CAZyme_pos_table <- merge(CAZyme_pos_table, toAdd_x)
  
  
  # Make the CAZyme matrix ggplot2able
  # add a Phylum column
  mat$Phylum <- row.names(mat)
  
  heatmap_data <- mat %>% 
    reshape2::melt(value.name = "Abundance"
                   , varnames = c("Phylum", "CAZyme")) %>%
    dplyr::left_join(Phylum_pos_table) 
  names(heatmap_data)[2] <- 'CAZyme'
  heatmap_data <- heatmap_data %>%
    dplyr::left_join(CAZyme_pos_table)
  
  # Add a CAZy_Type column to classify the domains
  heatmap_data$CAZy_Type <- gsub("endo-1,4-β-xylanase", "", heatmap_data$CAZyme)
  heatmap_data$CAZy_Type <- gsub("[0-9]*", "", heatmap_data$CAZy_Type)
  heatmap_data$CAZy_Type <- gsub("β-Glucosidase", "", heatmap_data$CAZy_Type)
  heatmap_data$CAZy_Type <- gsub("_", "", heatmap_data$CAZy_Type)
  
  # Make the ones which have repeat domain types in lump with the others 
  heatmap_data$CAZy_Type[heatmap_data$CAZy_Type == "CBM.CBM.GH"] <- "CBM.GH"
  heatmap_data$CAZy_Type[heatmap_data$CAZy_Type == "CBM.GH.GH"] <- "CBM.GH"
  
  
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
      
      , max(x_center + toAdd_x + 0.5 * width))
  ) + 
    0.1 * c(-1, 1) # extra spacing: 0.1
  
  #CAZyme_axis_limits[2] <- CAZyme_axis_limits[2] +2
  
  # Stop errors because of missing values where they aren't needed
  heatmap_data$x_center[is.na(heatmap_data$x_center)] <- 1
  heatmap_data$width[is.na(heatmap_data$width)] <- 1
  
  # Ensure labels plot correctly
  heatmap_data$x_center_true <- with(heatmap_data, x_center + toAdd_x)
  
  CAZyme_pos_table$x_center_true <- CAZyme_pos_table$x_center + 
    CAZyme_pos_table$toAdd_x
  
  #### Heatmap plot ####
  
  plt_hmap <- ggplot(heatmap_data, 
                     aes(x = factor(x_center_true), y = y_center, fill = Abundance, 
                         height = height, width = width)) + 
    geom_tile() +
    scale_fill_gradient(expression("Number of Genes")
                        , high = heatmap_colour, low = "white"
                        , trans = 'log10'
                        , breaks = colour_labels + 1
                        , labels = colour_labels) +
    scale_x_discrete(breaks = factor(CAZyme_pos_table$x_center_true), 
                     labels = CAZyme_pos_table$CAZyme, 
                     expand = c(0, 0)) + 
    # For the y axis, alternatively set the labels as: Phylum_position_table$Phylum
    scale_y_continuous(breaks = Phylum_pos_table[, "y_center"], 
                       labels = Phylum_pos_table$Phylum,
                       limits = Phylum_axis_limits, 
                       expand = c(0, 0)) + 
    labs(x = "CAZy Family", y = "Phylum") +
    theme_bw() +
    theme(axis.text.x = element_text(size = x_axis_size
                                     #, hjust = x_axis_hjust
                                     , angle = 90
                                     , vjust = x_axis_vjust),
          axis.text.y = element_text(size = y_axis_size),
          strip.text.x = element_text(size = strip_text_size, angle = 90), 
          # margin: top, right, bottom, and left
          plot.margin = unit(c(0.2, -0.7, 0.2, 0.2), "cm"), 
          panel.grid = element_blank(),
          legend.position = 'bottom',
          panel.spacing = unit(0, "lines")
    ) +
    facet_grid( ~ CAZy_Type
                , scales = "free"
                , space = 'free') 
  
  #### Dendrogram plot ####
  plt_dendr <- ggplot(segment_data) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
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
          , plot.margin = unit(c(0, 0.2, 0.2, dendrogram_left_margin), "cm")) 
  
  #### Top dendrogram  ####
  plt_dendr_horiz <- ggplot(segment_data_horiz) + 
    geom_segment(aes(x = x + toAdd_x, y = y
                     , xend = xend + toAdd_x, yend = yend)) + 
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
    ) +
    ggtitle(CapStr(metabolite))
  
  #### Blank Plot ####
  BLANK <- ggplot() + cowplot::theme_nothing() +
    #trbl
    theme(plot.margin = unit(c(0.2, 0.2, -0.7, -0.7)
                             , "cm"))
  
  ##### Put them all together ####
  
  metabolitePlot <- ggarrange(plt_dendr_horiz, BLANK,
                              plt_hmap, plt_dendr
                              , ncol = 2
                              , widths = c(1, side_dend_width)
                              , heights = c(top_dend_height, 1)
  )
  return(metabolitePlot)
}
