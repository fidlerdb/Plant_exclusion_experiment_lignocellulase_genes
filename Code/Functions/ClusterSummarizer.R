ClusterSummarizer <- function(metaboliteColumn
                              , metaboliteCAZyResults){
#metaboliteColumn = PhylumClusters$glucose
#metaboliteCAZyResults = GlucoseResults

  #### (1) Percentage of genes that a CAZy cluster takes up ####
  GeneCount <- lapply(unique(metaboliteColumn)
                     , FUN = function(x){
                       Phyla <- PhylumClusters$Phylum[metaboliteColumn == x]
                       NumberofGenes <- metaboliteCAZyResults$Gene_Richness[
                         row.names(metaboliteCAZyResults$Gene_Richness) %in% Phyla,]
                       
                       if(class(NumberofGenes) == "matrix"){
                         NumberofGenes <- data.frame(NumberofGenes)
                         Results <- colSums(NumberofGenes)
                       } else {Results <- NumberofGenes}
                       
                       Results <- matrix(Results
                                         , ncol = ncol(metaboliteCAZyResults$Gene_Richness)
                       )
                       Results <- data.frame(Results)
                       names(Results) <- dimnames(metaboliteCAZyResults$Gene_Richness)[[2]]
                       rownames(Results) <- paste(Phyla, collapse = " ")
                       return(Results)
                     })
  
  PhylumLabels <- lapply(GeneCount, rownames)
  GeneCount <- data.table::rbindlist(GeneCount)
  
  # Find percentage of genes from a phylum cluster which belong
  # to each CAZyme cluster
  PercentageGenes <- (GeneCount/rowSums(GeneCount))*100
  rownames(PercentageGenes) <- PhylumLabels
  
  # Find the mean number of genes per phylum in each cluster--10 is a cutoff
  TotalGenes <- rowSums(GeneCount)
  #rownames(TotalGenes) <- PhylumLabels
  NumberofTaxa <- unlist(lapply(strsplit(unlist(PhylumLabels), " "), length))
  MeanGenes <- TotalGenes/NumberofTaxa
  names(MeanGenes) <- unlist(PhylumLabels)
  names(TotalGenes) <- unlist(PhylumLabels)
  

  return(list(PercentageGenes = PercentageGenes
              , MeanGenes = MeanGenes
              , TotalGenes = TotalGenes))
}

