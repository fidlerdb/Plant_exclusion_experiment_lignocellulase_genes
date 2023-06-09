library(compiler)

# Make simprof a distance function it can use
Chi_dist <- function(x){d <- distance(x, method = "chi.square", dist = TRUE)
d[is.na(d)] <- 0
return(d)
}

Chi_dist <-  cmpfun(Chi_dist)

CAZymeClustering <- function(metabolite){
  
  #metabolite = "3,6-anhydro-D-galactose"
  
  #### (1) Count how many genes each phylum has for each CAZyme ####
  XyloseCAZyRichness <- lapply(1:nrow(cl[cl$Metabolite == metabolite,])
                               , FUN = function(i){
                                 table(df$Phylum[df$CAZyme == cl[cl$Metabolite == metabolite,]$CAZymes[i]])
                               })
  names(XyloseCAZyRichness) <- cl[cl$Metabolite == metabolite,]$CAZymes
  
  #### (2) Format the data nicely ####
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
  
  dimnames(mat)[[1]] <- XyloseRichness$Phylum
  
  SUMS <- colSums(mat)
  ZeroSumColumns <- names(SUMS[SUMS == 0])
  mat <- mat[, !dimnames(mat)[[2]] %in% ZeroSumColumns]
  
  ROWSUMS <- rowSums(mat)
  mat <- mat[ROWSUMS != 0,]
  
  mat <- t(mat) # Now the data is cleaned we can turn the matrix the 
  # right way round for clustering
  
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
               CAZyme = simprofResults$significantclusters[[x]])
  }))
  
  CAZymeRichness <- data.frame(Phylum = dimnames(mat)[[2]], t(mat))
  
  #### (5) Results ####
  
  # 5.1 How many genes were there in each cluster for each phylum?
  
  ClusterCAZymeRichness <- lapply(simprofResults$significantclusters
                                  , FUN = function(x){
                                    rowSums(data.frame(CAZymeRichness[,names(CAZymeRichness) %in% x]))}
  )
  
  names(ClusterCAZymeRichness) <- lapply(simprofResults$significantclusters
                                         , FUN = function(x){
                                           paste(x, collapse = " ")
                                         })
  
  ClusterCAZymeRichness <- do.call(cbind, ClusterCAZymeRichness)
  
  # 5.2 What was the total number of genes in each CAZyme cluster?
  TotalClusterRichness <- colSums(ClusterCAZymeRichness)
  
  #####(6) Final analysis results ####
  Results <- list(Clusters = Clusters
                  
                  # How many of which gene did each phylum have?
                  , Gene_Richness = ClusterCAZymeRichness
                  
                  # Total number of correlated CAZyme genes
                  , Total_Genes = TotalClusterRichness
                  
                  # Plottable dendrogram and Cluster composition
                  , simprofResults = simprofResults
  )
  
  return(Results)
}

CAZymeClustering <-  cmpfun(CAZymeClustering)