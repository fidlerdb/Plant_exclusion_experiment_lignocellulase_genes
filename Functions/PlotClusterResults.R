PlotClusterResults <- function(GeneCount, nrow = max(gcl$Cluster), ncol = ncol){
  library(ggplot2)
  gcl <- melt(data.frame(Cluster = 1:length(unique(rownames(GeneCount)))
                       , Phylum = rownames(GeneCount)
                       , GeneCount)
            , id.vars = c("Phylum", "Cluster"))
  gcl$Phylum <- sub(" .*", "", as.character(gcl$Phylum))
  gcl$CAZy_Cluster <- sub("\\..*", "", as.character(gcl$variable))

  ggplot(gcl, aes(x = CAZy_Cluster
                  , y = value
                  )) + 
  geom_bar(stat = "identity") +
  facet_wrap(.~ Phylum, nrow = nrow, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90
                                   , hjust = 1
                                   , vjust = 0.5))
}
  
