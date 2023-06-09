### Create a plot showing CAZy richness per phylum for cellulases etc., as   ###
### well as number of species per phylum                                     ###

#### Data import and cleaning ####
# .rs.restartR()
rm(list = ls()) # Clean the environment

# Import useful packages
library(data.table)
library(ggplot2)
library(magrittr)
#library(scales)

tf <- fread("Data/Henfaes_Species_Multi-Method_TPM_contigs_1kbpmin.csv", nThread = 3) # Taxonomy of all contigs
cfc <- fread("Data/Henfaes_CAZy_Family_SP_TPM_of_all_DNA_contigs.csv") # Total abundance of each CAZy family in a sample

cfc
samplecols <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7"          
, "c1", "c2", "c3", "c4", "c5", "c6", "c7")

### Keep only the relevant data for the CAZymes
colstokeep_cfc <- c("GenestoSearch", "CAZyme", "Chr"
                    #, "n1", "b1", "b2", "b3", "b4", "b5", "b6", "b7"          
                    #, "c1", "c2", "c3", "c4", "c5", "c6", "c7"
                    )
cfc <- cfc[, ..colstokeep_cfc]
setnames(cfc, old = c("GenestoSearch", "Chr"), new = c("Gene_ID", "Contig"))

### Clean the taxonomy data

tf2 <- tf # Create a column tyo view the number of species per phylum

tf2[is.na(Phylum),]
tf2[Phylum == "UNASSIGNED",]
tf2[Phylum == "UNASSIGNED",]

# remove superfluous columns
# Keep the reads for these contigs. If a contig with a CAZy gene on was mapped to
# it makes sense that that gene was present in the sample
names(tf)
colstokeep_tf <- c("Contig", "SuperKingdom", "Phylum", "Class"
                   , "Order", "Family", "Genus", "Species"
                   , "n1", "b1", "b2", "b3", "b4", "b5", "b6", "b7"          
                   , "c1", "c2", "c3", "c4", "c5", "c6", "c7"
                   )
tf[, names(tf)[c(!names(tf) %in% colstokeep_tf)] := NULL]
setnames(tf, old = "SuperKingdom", new = "Domain")

taxcols <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Make all unassigned taxonomic ranks obvious
tf[, (taxcols) := lapply(.SD, FUN = function(x){ifelse(x == ""
                                                       , yes = "UNASSIGNED"
                                                       , no = x)})
   , .SDcols = taxcols]

### Merge the taxonomy data to the CAZy data
tc <- merge(tf, cfc, all = TRUE)
tc[is.na(Domain),] # NAs are from contigs which weren't given any classification
tc[is.na(Phylum),] # NAs are from contigs which weren't given any classification
tc[!is.na(Domain),]

tc[,Genus]

100* (nrow(tc[!is.na(Domain),])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Phylum != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Class != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Order != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Family != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Genus != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Species != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes


tc[, Taxonomy := apply(.SD, 1, FUN = function(x){paste(x, collapse = ";")}), .SDcols = taxcols]
tc

sum(tc$n1, na.rm = TRUE) # I have done the cleaning job properly :)

# Get rid of the contigs that don't have a CAZyme on them
tc <- tc[!is.na(CAZyme), ]
tc <- tc[!is.na(Domain), ]
tc <- tc[Domain != "UNASSIGNED",]
tc <- tc[Phylum != "UNASSIGNED",]

### Find the number of genes mapped to from each phylum across all treatments

# Sum abundance of common cellulase CAZy families
common_cellulases <- c("GH5", "GH6", "GH7", "GH8", "GH9", "GH12"
                       , "GH44", "GH45", "GH48")

# GH8 (11%, 5%), GH10 (96%, 0.5%), GH11 (100%, 0.3%), GH30 (37%, 11%)
common_xylanases <- c("GH10", "GH11", "GH8", "GH30")

# https://www.pnas.org/doi/10.1073/pnas.2008888118#sec-1
LPMOs <- c("AA9", "AA10", "AA11", "AA13", "AA14", "AA15")

Auxiliary_cols <- unique(cfc$CAZyme)[grep("AA", unique(cfc$CAZyme))]

CBM_cols <- unique(cfc$CAZyme)[grep("CBM", unique(cfc$CAZyme))]

strict_lignocellulases <- c(common_cellulases, common_xylanases
                            , LPMOs, Auxiliary_cols
                            #, CBM_cols
                            )

get_regex <- function(CAZy_fams){sapply(CAZy_fams, FUN = function(x){
  paste0(x,"$", "|" # family end of string
         , x, "\\|", "|" # Family another family
         , x, "_" # family subfamily
  )
})}

regex_list <- lapply(list(regex_cel = common_cellulases
            , regex_xyl = common_xylanases
            , regex_LP = LPMOs
            , regex_AA = Auxiliary_cols
            , regex_CBM = CBM_cols
            , regex_ligno = strict_lignocellulases), FUN = function(s){
              get_regex(s)
            })


active_families_list <- lapply(regex_list,
                               FUN = function(x){
                                 out <- unlist(sapply(x, FUN = function(r){
                                   unique(cfc$CAZyme)[grep(r, unique(cfc$CAZyme))]
                                 }))
                                 names(out) <- NULL 
                                 return(out)
                               })
names(active_families_list) <- c("Cellulases", "Xylanases"
                                 , "LPMOs", "Auxiliary Activities"
                                 , "CBMs"
                                 , "Lignocellulases")

unlist(active_families_list)[anyDuplicated(unlist(active_families_list))]
# One shared cellulolytic and xylanolytic family (GH8)

# The object lignocellulases is the vector of all CAZy families containing one of the lignocellulase families in some way

## Count the number of genes in each of the specified CAZy families, at the phylum level

# Make the funtion to do this
get_Richness <- function(data, CAZy_to_use){
  # Make the dataset boolean at the gene level
  g <- data[CAZyme %in% CAZy_to_use,]
  
  g[, (samplecols) := lapply(.SD, function(x){ifelse(x == 0, yes = 0, no = 1)})
    , .SDcols = samplecols]

  # Boolean "was this CAZyme in this phylum in this sample?" column
  g[, Overall := ifelse(rowSums(.SD) == 0, yes = 0, no = 1)
     , .SDcols = samplecols]

  # Sum the number of CAZy genes within each phylum
  g <- g[, .(Richness = sum(Overall))
           , by = c("Phylum", "Domain")]
  
  return(g)

  }

# Count the number of gene
cr <- list(Cellulases = get_Richness(data = tc, CAZy_to_use = active_families_list$Cellulases)
     , Xylanases = get_Richness(data = tc, CAZy_to_use = active_families_list$Xylanases)
     , LPMOs = get_Richness(data = tc, CAZy_to_use = active_families_list$LPMOs)
     , AAs = get_Richness(data = tc, CAZy_to_use = active_families_list$`Auxiliary Activities`)
     , CBMs = get_Richness(data = tc, CAZy_to_use = active_families_list$CBMs)
     , Lignocellulases = get_Richness(data = tc, CAZy_to_use = active_families_list$Lignocellulases)
     , Species = tc[, .(Richness = length(unique(Taxonomy))), by = c("Domain", "Phylum")])

# Fix the names for a multi-merge
lapply(seq_along(cr), FUN = function(i){
  setnames(cr[[i]], old = "Richness", new = names(cr)[i])
})

# Stick these variables side-by-side
cr <- Reduce(f = function(dt1, dt2){merge(dt1, dt2, by = c("Phylum", "Domain")
                                          , all = TRUE)}
       , x = cr)

# Put it in log-format for visualization, and tidy it up for plotting
CAZy_cols <- c("Cellulases", "Xylanases", "LPMOs", "AAs", "CBMs", "Lignocellulases", "Species")
cr[, (CAZy_cols) := lapply(.SD, FUN = function(x){ifelse(is.na(x), yes = 0, no = x)})
   , .SDcols = CAZy_cols] # Remove NAS
cr_w <- cr
cr[, colorder := frank(Species, ties.method = "random")] # Order the columns by species abundance
cr <- melt(cr, id = c("Domain", "Phylum", "colorder")
           , variable.name = "Variable", value.name = "Richness") # Long format

CAZy_cols <- c("Cellulase genes", "Xylanase genes"
               , "LPMO genes", "AA genes"
               , "CBM genes", "Lignocellulase genes", "Species")

# Tidy up names
cr[Variable == "Cellulases", Variable := "Cellulase genes"]
cr[Variable == "Xylanases", Variable := "Xylanase genes"]
cr[Variable == "AAs", Variable := "AA genes"]
cr[Variable == "CBMs", Variable := "CBM genes"]
cr[Variable == "Lignocellulases", Variable := "Lignocellulase genes"]
cr[, Variable := factor(Variable, levels = c("Cellulase genes", "Xylanase genes"
                                             , "LPMO genes", "AA genes"
                                             , "CBM genes", "Lignocellulase genes"
                                             , "Species"))]

cr[Variable == "Species", Richness := Richness +1]

#design <- matrix(c(1,2,5,3,4,5), 3, 2)
design <- matrix(c(1,3,5,2,4,6), 3, 2)

scales <- list(
  scale_x_continuous(expand = c(0, 0), limits = c(0,55)),
  scale_x_continuous(expand = c(0, 0), limits = c(0,55)),
  scale_x_continuous(expand = c(0, 0), limits = c(0,55)),
  scale_x_continuous(expand = c(0, 0), limits = c(0,500)),
  scale_x_continuous(expand = c(0, 0), limits = c(0,125)),
  scale_x_log10()#expand = c(0, 0))#, limits = c(0,500), trans = "log10")
)


library(ggh4x)
library(cowplot)
# Palette
DomainPalette <- c('#e41a1c'# Archaea
                   , '#377eb8'# Bacteria
                   , '#4daf4a'# Eukaryotes
)

p <- ggplot(cr[Variable != "LPMOs",], aes(x = Richness, y = reorder(Phylum, colorder), fill = Domain)) +
  geom_bar(stat = "identity") +
  facet_manual(~Variable, design = design, scales = "free_x"
               , axes = "y", remove_labels = "y") +
  facetted_pos_scales(x = scales) +
  scale_fill_manual(values = DomainPalette) +
  ylab("Phylum") +
  theme_cowplot() + 
  theme(legend.position = "bottom"
        , axis.text.y = element_text(size = 9)
        , axis.text.x = element_text(size = 10)
        , legend.text = element_text(size = 10)
        , legend.title = element_text(size = 12)
        , axis.title = element_text(size = 12)
        , plot.margin = margin(2,5,2,2, "mm")
        , panel.grid.major.x = element_line(colour = "grey70")
        )
  
p

# Save the figure
tiff(filename = 'Figure_4_Lignocellulase_gene_origins.TIFF'
     , width = 15.5, height = 22
     , units = "cm"
     , res = 300)
p
dev.off()

cr_l <- melt(cr_w[, c('Phylum', 'Domain', 'Cellulases', 'Xylanases', 'AAs' , 'CBMs', 'Species')]
     , id = c("Domain", "Phylum", "Species")
     , variable.name = "Gene_type", value.name = "Richness")

ggplot(cr_l, aes(x = Species + 1, y = Richness + 1, label = Phylum)) +
  geom_smooth(method = "lm") +
  geom_text() +
  scale_y_log10() +
  scale_x_log10() + 
  facet_wrap(~Gene_type) 


# species-CAZy relationships
ggplot(cr_w, aes(x = Species+1, y = Cellulases+1, label = Phylum)) + 
  geom_smooth(method = 'lm') +
  geom_text() +
  scale_y_log10() +
  scale_x_log10()
# Acidobacteria have a high cellulase:species

ggplot(cr_w, aes(x = Species+1, y = Xylanases+1, label = Phylum)) + 
  geom_smooth(method = 'lm') +
  geom_text()+
  scale_y_log10() +
  scale_x_log10()
# Acidobacteria have a high xylanase:species

ggplot(cr_w, aes(x = Species+1, y = AAs+1, label = Phylum)) + 
  geom_smooth(method = 'lm') +
  geom_text() +
  scale_y_log10() +
  scale_x_log10()
# Most points follow the line very well

ggplot(cr_w, aes(x = Species+1, y = CBMs+1, label = Phylum)) + 
  geom_smooth(method = 'lm') +
  geom_text()+
  scale_y_log10() +
  scale_x_log10()
# Actinobacteria have a high CBM:species, as do acidobacteria

#### 2022-09-30 Taxonomy of specific CAZy genes ####

# Get the gene families we are interested in
tc
gh5_subfamiles <- c(1,5,25,46,36,13,19)
get_regex <- function(family, subfamily){
  paste0("GH", family, "_", subfamily,"$|GH",family, "_", subfamily, "\\|")
}
gh5_regex <- lapply(gh5_subfamiles, FUN = get_regex, family = "5")

tc_gh5_1 <- tc[grep(gh5_regex[[1]], CAZyme), ] # GH5_1
tc_gh5_5 <- tc[grep(gh5_regex[[2]], CAZyme), ] # GH5_1
tc_gh5_25 <- tc[grep(gh5_regex[[3]], CAZyme), ] # GH5_1
tc_gh5_46 <- tc[grep(gh5_regex[[4]], CAZyme), ] # GH5_1
tc_gh5_36 <- tc[grep(gh5_regex[[5]], CAZyme), ]
tc_gh5_13 <- tc[grep(gh5_regex[[6]], CAZyme), ]
tc_gh5_19 <- tc[grep(gh5_regex[[7]], CAZyme), ]

# Subset to CAZy genes in these familes
tcr <- list(tc_gh5_1[, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc_gh5_5[, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc_gh5_25[, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc_gh5_46[, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc_gh5_46[, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc_gh5_36[, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc_gh5_13[, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc_gh5_19[, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("GH30$|GH30_1|GH30_3|GH30\\|", CAZyme),][, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("AA10$|AA10_|AA10\\|", CAZyme),][, .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("GH8$|GH8_|GH8\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("GH11$|G11_|G11\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("GH10$|G10_|G10\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("GH30_2$|GH30_2\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("AA3$|AA3\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("AA3_2$|AA3_2\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("AA5$|AA5\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("AA6$|AA6\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("AA7$|AA7\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("AA10$|AA10\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
     , tc[grep("AA12$|AA12\\|", CAZyme), .N, keyby = c("CAZyme", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
) %>% rbindlist(.)
tcr

#fwrite(tcr, "Species_Lignocellulases_Number_of_contigs.csv")

# How many contigs were there?
tcr[, sum(N)]

# Now find the number of genes in any given taxonomic rank, for each CAZy subfamily
get_unique_classifications_genes <- function(Rank){
  tcr[, .(Unique_classifications = length(unique(paste(Phylum,  Class, Order, Family, Genus, Species)))) 
      , by = c("CAZyme", Rank)] #%>%
    #dcast(., CAZyme ~ get(Rank)
          #, fill = 0
          #, value.var = 'Unique_classifications')
}

B10y_czm <- get_unique_classifications_genes("Phylum") %>%
  tidyr::complete(., Phylum, CAZyme) %>% setDT(.)

source("Functions/Get_CAZy_family_ranks.R")
B10y_czm[, CAZyme := factor(CAZyme, levels = get_gene_order(CAZyme))]

p <- ggplot(B10y_czm, aes(x = Phylum, y = CAZyme, fill = Unique_classifications)) +
  #geom_tile(fill = "white", colour = "black") +
  geom_tile(colour = "black"
            , na.rm = TRUE) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_continuous(name = "Number of taxa\nwith gene", na.value = 'salmon') +
  ylab("CAZy classification")

p

tiff(filename = 'Figure_5_10-year_bare_associated_CAZy.TIFF'
     , width = 15.5, height = 15.5
     , units = "cm"
     , res = 600)
p
dev.off()

#
tcr[CAZyme == "GH5_36",]
tcr[CAZyme == "GH5_13",]

tcr[CAZyme == "GH5_13", sum(N), by = c("Phylum", "Class", "Order", "Family")]
tcr[CAZyme == "GH5_13", sum(N), by = c("Phylum")]

tcr[CAZyme == "GH5_19",]

#### Where did the 10-year bare treatment's cellulases come from? ####
tcr[CAZyme %in% c("GH5_1", 'GH5_5', "GH5_25", "GH5_46", "GH5_36", "GH5_13", "GH5_19")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order")]

# More specifically within acidobacteria
tcr[CAZyme %in% c("GH5_1", 'GH5_5', "GH5_25", "GH5_46", "GH5_36", "GH5_13", "GH5_19")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Acidobacteria",]

tcr[CAZyme %in% c("GH5_1", 'GH5_5', "GH5_25", "GH5_46", "GH5_36", "GH5_13", "GH5_19")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Actinobacteria",]

tcr[CAZyme %in% c("GH5_1", 'GH5_5', "GH5_25", "GH5_46", "GH5_36", "GH5_13", "GH5_19")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Bacteroidetes",]

tcr[CAZyme %in% c("GH5_1", 'GH5_5', "GH5_25", "GH5_46", "GH5_36", "GH5_13", "GH5_19")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Planctomycetes",]

tcr[CAZyme %in% c("GH5_1", 'GH5_5', "GH5_25", "GH5_46", "GH5_36", "GH5_13", "GH5_19")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Proteobacteria",]

#### Where did the 10-year bare treatment's xylanases come from? ####
tcr[CAZyme %in% c("GH30", 'GH30_1', "GH30_3")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order")]

# More specifically within acidobacteria
tcr[CAZyme %in% c("GH30", 'GH30_1', "GH30_3")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Acidobacteria",]

tcr[CAZyme %in% c("GH30", 'GH30_1', "GH30_3")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Actinobacteria",]

tcr[CAZyme %in% c("GH30", 'GH30_1', "GH30_3")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Bacteroidetes",]

tcr[CAZyme %in% c("GH30", 'GH30_1', "GH30_3")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Firmicutes",]

tcr[CAZyme %in% c("GH30", 'GH30_1', "GH30_3")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order", "Family", "Genus", "Species")][Phylum == "Proteobacteria",]

#### Where did the 10-year bare treatment's auxiliary activities come from? ####
tcr[CAZyme %in% c("AA10|CBM73")
    , sum(N)
    , keyby = c("Phylum", "Class", "Order")]

#### Make an interpretable table version of this ####

tt <- tcr[CAZyme %in% c("GH30", 'GH30_1', "GH30_3", "AA10|CBM73"
                        , "GH5_1", 'GH5_5', "GH5_25", "GH5_46", "GH5_36", "GH5_13", "GH5_19")]

tt[, `Activity grouping` := "Cellulase"]
tt[CAZyme %in% c("GH30", 'GH30_1', "GH30_3"), `Activity grouping` := "Xylanase"]
tt[CAZyme %in% c("AA10|CBM73"), `Activity grouping` := "Auxiliary Activity"]


tt_o <- tt[, .(`Number of genes` =  sum(N))
    , keyby = c("Activity grouping", "Phylum", "Class", "Order", "Family", "Genus", "Species")]
tt_o <- tt_o[, lapply(.SD, FUN = function(x){ifelse(x == "UNASSIGNED", yes = "Unassigned", no = x)})]

tt_o



