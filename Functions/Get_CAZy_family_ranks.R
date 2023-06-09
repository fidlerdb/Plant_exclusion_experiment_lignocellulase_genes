
# Create a function to order CAZy genes nicely

if(!require(data.table)){require(data.table)}
if(!require(magrittr)){require(magrittr)}


#cellulases <- c("GH9", "GH8", "GH44", "GH5_35", "GH5_7", 
#  "GH5_8", "GH5_13", "GH5_5", "GH6", "GH5_46", "GH5_26", "GH5", 
#  "GH5_24", "GH5_40", "GH5_19", "GH5_25", "GH5_4", "GH5_28", "GH5_27", 
#  "GH12", "GH5_36", "GH45", "CBM2|GH5_1", "GH5_22", "CBM4|GH9", 
#  "CBM2|GH44", "CBM8|GH44", "CBM6|GH5_46", "GH5_2|GH5_1"
#  , "GH6_1|GH5_2", "GH6|GH5_1|GH5_2", "GH6|GH5_1|GH5_1")

get_CAZy_order <- function(gene_families
                           , class_of_interest){
  if(!require(magrittr)){library(magrittr)} # Load package if necessary
  # Split the gene into domains
  domains <- strsplit(gene_families, "\\|") %>%
    lapply(., FUN = function(x){
      x[grep(class_of_interest, x)]
    })
  
  # Domain numbers
  domain_numbers <- sub(class_of_interest, "", domains[[1]]) %>% # Get rid of the enzyme class
    sub("_.*", "", .) %>% # Get rid of the subfamily
    as.numeric(.)
  
  subfamilies <- sub(class_of_interest, "", domains[[1]]) %>% # Get rid of the enzyme class
    sub("^[A-Z0-9]*", "", .) %>%
    sub("_", "", .) %>%
    as.numeric(.)
  subfamilies[is.na(subfamilies)] <- 0 # Any family without a subfamily is assigned as the top
  
  data.table(Domains = gene_families
             , Family_number = domain_numbers
             , Subfamily_number = subfamilies
             , Family_rank = domain_numbers + (subfamilies / 1e2)
             , Domain_number = 1:length(domain_numbers))
  }

get_gene_order <- function(gene_families, class_of_interest){
  # Get the rank of each domain, and what domain number it was
  family_ranks <- lapply(gene_families, get_CAZy_order, class_of_interest = "GH") %>% 
  rbindlist(.)
  
  # Make the nth domain much less important than the 1st domain
  family_ranks[, Domain_rank := Family_rank ^ (-Domain_number)]

  # Combine the information about each CAZy gene, rankign nthem based off of
  # Family number. Then spit them out in the order they need to be in
  out <- family_ranks[, sum(Domain_rank), by = "Domains"]
  setkey(out, "V1")
  

  # Spit out the corect order of genes
  return(rev(out[, Domains]))
}

