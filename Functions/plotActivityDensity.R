#### Create a function to plot activity frequency to match up with the heatmaps produced earlier ####

# Clear the workspace
#rm(list = ls())

# Load packages
library(ggplot2)
library(data.table)



# Import data
#ct <- fread("CAZy/05_CAZyme_Functional_information/OutputData/CAZy_Cluster_Activity_Summary.csv")



# Make the data ggplottable
#cl <- melt(ct, id = c("Cluster", "nFamilies"))

# Select a metabolite

plotActivityDensity <- function(metabolite){
  
  #clg <- ctl[grep(metabolite, ctl$Cluster),]
  #clg$Order <- as.numeric(sub(".*CazC", "", clg$Cluster)) # Make data table ordering easy
  #clg <- clg[order(Order)] # Order the data sensibly (for indexing in a bit)
  
  # stop benzoic acid matching 4-hydroxybenzoic acid
  if(metabolite == "benzoic acid"){metabolite <- paste0("^", metabolite)}
   
  # Keep only relevant data
  clg <- ctl[grep(metabolite, ctl$Cluster, perl = TRUE),]
  
  
  clg$Order <- sub(".*CazC", "", clg$Cluster) # Make data table ordering easy
  
  if(any(grepl("_", clg$Order))){
    # Find any that are not in a cluster
    NoCluster <- unique(sub("_", "", clg[!grep("^[0-9]", Order)]$Order))
    
    # Replace cluster order with the right value 
    clg[grep("^[0-9]", clg$Order),]$Order <-  setdiff(1:length(NoCluster), NoCluster)
    clg$Order <- as.numeric(sub("_", "", clg$Order))
    clg <- clg[order(Order)] # Order the data sensibly (for indexing in a bit)
    
  } else {
    clg <- ctl[grep(metabolite, ctl$Cluster),]
    clg$Order <- as.numeric(sub(".*CazC", "", clg$Cluster)) # Make data table ordering easy
    clg <- clg[order(Order)] # Order the data sensibly (for indexing in a bit)
  }
  # stop benzoic acid not matching 4-hydroxybenzoic acid being a problem
  if(metabolite == "^benzoic acid"){metabolite <- "benzoic acid"}
  
  # Repeat the values nFamilies times so i can be added to the heatmap graphic. 
  # Add an index for plotting
  clg <- clg[rep(1:.N,nFamilies)][,xval:=1:.N, by = variable]
  
  # Ensure things plot how I want them to
  clg$variable <- factor(clg$variable
                         , levels=c("Oligosaccharides", "Hemicellulose", "Cellulose", "Lignin"))
  
  # Create the plot
  p <- ggplot(data = clg
         , aes(x = xval, y = variable, fill = value
               , height = 1, width = 1)
  ) + geom_tile() +
    scale_fill_gradient(low = "white", high = 'black') +
    scale_x_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(axis.text.x = element_blank()
          , axis.title = element_blank()
          , panel.grid = element_blank()
          , plot.margin = unit(c(0.0, -0.7, 0.0, 0.2), "cm")
          , panel.spacing = unit(0, "lines")
          #, legend.position = "none"
          , panel.border = element_rect(fill = NA)
          ) +
    guides(fill = guide_colourbar(direction = 'horizontal'
                                , barwidth = 4
                                #, barheight = 3
                                , title = "% CAZy families with\nrelated activities"
                                , title.position = "top"
                                ))
  return(p)
}


# "#0D066C" # Cellulose
# "#6C1706" # Lignin
# "#066C0D" # Hemicellulose
# "cyan3" # Oligosaccharides