
plotClusterHeatmap <- function(df, x_axis_size = 4, y_axis_size = 4
         , dendrogram_left_margin = -0.6, x_axis_vjust = 0.5
         , colour_labels =  c(0,10,100,400), heatmap_colour = "#152736"
         , colour_scale_labels = c(0, 10, 100, 1000), side_dend_width = 0.5
         , top_dend_height = 0.5){

library(ggrepel)

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

#### Heatmap plot ####

plt_hmap <- with(df$heatmap, (ggplot(heatmap_data, 
                   aes(x = factor(x_center), y = y_center
                       , fill = Abundance, 
                       height = height, width = width)) + 
  geom_tile() +
  scale_fill_gradient(expression("Number of Genes")
                      , high = heatmap_colour, low = "white"
                      , trans = 'log10'
                      , breaks = colour_labels + 1
                      , labels = colour_labels) +
  scale_x_discrete(breaks = factor(CAZyme_pos_table$x_center), 
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
        axis.text.y = element_text(size = y_axis_size
                                   , colour = Phylum_pos_table$Colour),
        # margin: top, right, bottom, and left
        plot.margin = unit(c(0.2, -0.7, 0.2, 0.2), "cm"), 
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
        , plot.margin = unit(c(0, 0.2, 0.2, dendrogram_left_margin), "cm")
        , legend.position = "none"
  ) +
  # geom_text(data = Phylum_dend_node_data[ClusterNodes_Phylum,]
  #           , aes(x = y + 1, y = x + 1
  #                 , label = text
  #                 , size = 6
  #                 , fontface = "bold"
  #                 #, col = na.omit(unique(segment_data$col))
  #           ), colour = 'grey50'
  # )
  geom_text_repel(data = Phylum_dend_node_data[ClusterNodes_Phylum,]
                  , aes(x = y , y = x
                        , label = text
                        , size = 6
                        , fontface = "bold"
                        ), colour = 'grey50'
                  , nudge_x = 1
                  #, segment.size  = 0.2
                  #, segment.color = "grey50"
                  #, direction     = "x"
  )
)

#### Top dendrogram  ####
plt_dendr_horiz <- with(df$cazyme_dend, 
                        ggplot(segment_data_horiz) + 
  geom_segment(aes(x = x 
                   , y = y
                   , xend = xend 
                   , yend = yend
                   , col = col)) + 
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
  ggtitle(CapStr(df$metabolite)) +
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

metabolitePlot <- ggarrange(plt_dendr_horiz, BLANK,
                            plt_hmap, plt_dendr
                            , ncol = 2
                            , widths = c(1, side_dend_width)
                            , heights = c(top_dend_height, 1)
)
return(metabolitePlot)
}
