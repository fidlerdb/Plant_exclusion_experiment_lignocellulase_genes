library(compiler)

# Make simprof a distance function it can use
Chi_dist <- function(x){d <- distance(x, method = "chi.square", dist = TRUE)
d[is.na(d)] <- 0
return(d)
}

Chi_dist <-  cmpfun(Chi_dist)

PhylumClustering <- function(metabolite){
  
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
  
  # Set up data for plotting a heatmap with dendogram
  
  # Create a CAZyme abundance matrix
  mat <- XyloseRichness[2:length(XyloseRichness)]
  
  #mat <- log10(mat + 1)
  SUMS <- colSums(mat)
  ZeroSumColumns <- names(SUMS[SUMS == 0])
  mat <- mat[, !names(mat) %in% ZeroSumColumns]
  
  rownames(mat) <- XyloseRichness$Phylum
  CAZyme_names <- colnames(mat)
  
  #### (3) Create the distance matrix for phyla  ####
  
  # calculate the distance matrix
  d <- distance(mat, method = "chi.square", dist = TRUE)
  d[is.na(d)] <- 0 # Doesn't like fully zero rows
  #, so replace these distances with 0
  # what the metric is trying to do:
  # sqrt((0-0)^2/0+0))
  
  #### (3) Perform the heirarchical clustering and make the dendrogram ####
  
  metaboliteClust <- hclust(d, method = "ward.D2")
  dend <- as.dendrogram(metaboliteClust)
  
  #plot(dend)
  
  #### (4) Significance testing of clusters ####
  
  # Perform the test
  simprofResults <- simprof(data = mat, method.cluster = 'ward.D2'
                            , method.distance = Chi_dist
  )
  
  # Easier viewing
  Clusters <- rbindlist(lapply(seq_along(simprofResults$significantclusters), FUN = function(x){
    data.frame(Cluster = x,
               Phylum = simprofResults$significantclusters[[x]])
  }))
  
  
  # Easy table-format for XyloseCAZyRichness
  Gene_Richness <- rbindlist(lapply(XyloseCAZyRichness
                           , function(y) {lapply(split(y,names(y))
                                          , function(x) {Reduce("+", x)})
                             })
                    , idcol = "CAZyme"
                    , fill = TRUE)
  Total_Genes <- colSums(Gene_Richness[,!"CAZyme", with = FALSE]
                         , na.rm = TRUE)
  
  # What is the mean number of genes in each cluster?
  Cluster_Mean_Genes <- lapply(simprofResults$significantclusters
                               , FUN = function(x){
                                 return(
                                   mean(Total_Genes[names(Total_Genes) %in% x])
                                   )
                                 })
  
  
  
  # Analysis results
  Results <- list(Clusters = Clusters
                  
                  # How many of which gene did each phylum have?
                  , Gene_Richness = Gene_Richness
                  
                  # Total number of correlated CAZyme genes
                  , Total_Genes = Total_Genes
                  
                  # Mean number of genes in each cluster
                  , Cluster_Mean_Genes = Cluster_Mean_Genes
                  
                  # Plottable dendrogram and Cluster composition
                  , simprofResults = simprofResults
                  )
  
  return(Results)
}

PhylumClustering <-  cmpfun(PhylumClustering)

