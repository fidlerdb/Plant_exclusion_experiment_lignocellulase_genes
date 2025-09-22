### Create a plot showing CAZy richness per phylum for cellulases etc., as   ###
### well as number of species per phylum                                     ###

#### Data import and cleaning ####
# .rs.restartR()
rm(list = ls()) # Clean the environment

# Import useful packages
library(data.table)
library(ggplot2)
library(magrittr)
#library(glmmTMB)
#library(emmeans)
library(vegan)
library(patchwork)
#library(scales)

tf <- fread("Taxonomy/Multi-Method/Cleaned_Data/Henfaes_Species_Multi-Method_TPM_contigs_1kbpmin.csv", nThread = 3) # Taxonomy of all contigs

cfc <- fread("CAZy/Cleaned_Data/Henfaes_CAZy_Family_SP_TPM_of_all_DNA_contigs.csv") # Total abundance of each CAZy family in a sample
cfc_full <- copy(cfc)
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

# Ensure everything has some sort of classification
# tf2[SuperKingdom == "", SuperKingdom := "UNASSIGNED"]
# tf2[Phylum == "", Phylum := "UNASSIGNED"]
# tf2[Class == "", Class := "UNASSIGNED"]
# tf2[Order == "", Order := "UNASSIGNED"]
# tf2[Family == "", Family := "UNASSIGNED"]
# tf2[Genus == "", Genus := "UNASSIGNED"]
# tf2[Species == "", Species := "UNASSIGNED"]

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

# Add values from the contig file below where we have NA values(and reads not mapped to negative)
tc[Contig == 'k141_2024747', (samplecols) := .(b1 = 0.270937, b2 =0, b3=0
                                               ,b4=0,b5= 0.2189873,b6= 0.3863135,b7= 0.2322574
                                               ,c1 =0,c2 =0,c3=0
                                               ,c4=0,c5=0,c6=0,c7=0)]

# tc[CAZyme == 'GH30_4', Taxonomy := "UNASSIGNED;UNASSIGNED;UNASSIGNED;UNASSIGNED;UNASSIGNED;UNASSIGNED;UNASSIGNED"]

# Replace NA with unassigned
tc[, (taxcols) := lapply(.SD, FUN = function(x){ifelse(is.na(x)
                                                       , yes = "UNASSIGNED"
                                                       , no = x)})
   , .SDcols = taxcols]

head(tc[,Genus])

100* (nrow(tc[!is.na(Domain),])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Phylum != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Class != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Order != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Family != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Genus != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes
100* (nrow(tc[Species != "UNASSIGNED",])/nrow(tc)) # Percentage assigned CAZymes


tc[, Taxonomy := apply(.SD, 1, FUN = function(x){paste(x, collapse = ";")}), .SDcols = taxcols]
tc
tc[CAZyme == 'GH30_4',]# Taxonomy]
sum(tc$n1, na.rm = TRUE) # I have done the cleaning job properly :)

# Find our Halogeometricum
tc[grep("Halogeometricum borinquense", Taxonomy),]

#### How many taxa had lignocellulases? ####

Ntaxa <- tc[, unique(Taxonomy)] %>% length(.)
Ntaxa_CAZyme <- tc[!is.na(CAZyme), unique(Taxonomy)] %>% length(.)

genes_lignocellulases <- c("GH5_35", "GH5_7", "GH5_8", "GH5_13", "GH5_5", "GH5_46", "GH5_26", 
                           "GH5", "GH5_24", "GH5_40", "GH5_19", "GH5_25", "GH5_4", "GH5_28", 
                           "GH5_27", "GH5_36", "CBM2|GH5_1", "GH5_22", "CBM6|GH5_46", "GH6", 
                           "GH8", "GH9", "CBM4|GH9", "GH12", "GH44", "CBM2|GH44", "CBM8|GH44", 
                           "GH45", "GH10", "GH11", "CBM60|GH11", "GH8", "GH30", "GH30_2", 
                           "GH30_4", "GH30_1", "CBM32|GH30_3", "AA10", "AA10|CBM73", "AA12", 
                           "AA7", "AA10", "AA10|CBM73", "AA3", "AA3_2", "AA3_2", "AA6", 
                           "AA5", "AA10", "AA10|CBM73")

Ntaxa_Lignocellulase <- tc[!is.na(CAZyme) & 
                                 CAZyme %in% genes_lignocellulases
                               , unique(Taxonomy)] %>% length(.)
100*(Ntaxa_Lignocellulase/Ntaxa) %>% round(., 3)

#### Back to data cleaning ####

# Get rid of the contigs that don't have a CAZyme on them

tc_all <- copy(tc)
tc <- tc[!is.na(CAZyme), ]
#tc <- tc[!is.na(Domain), ]
#tc <- tc[Domain != "UNASSIGNED",]
#tc <- tc[Phylum != "UNASSIGNED",]

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

#### Abundance increased in bare 2024-07-24 ####

# Get the lignocellulase genes for analysis (abundance) and format this for analysis
tc_c <- tc_all[CAZyme %in% active_families_list$Lignocellulases,]
colstokeep <- c(samplecols, 'CAZyme')
tc_c <- melt(tc_c[, ..colstokeep], id = "CAZyme"
             , value.name = "CPM"
             , variable.name = "Sample")
tc_c[, Treatment := "10-Year Bare"]
tc_c[grep("b[4-7]", Sample), Treatment := "1-Year Bare"]
tc_c[grep("c[1-3]", Sample), Treatment := "10-Year Grassland"]
tc_c[grep("c[4-7]", Sample), Treatment := "1-Year Grassland"]
tc_c[, Treatment := factor(Treatment, levels = c("1-Year Bare"
                                                 , "10-Year Bare"
                                                 , "1-Year Grassland"
                                                 , "10-Year Grassland"))]
tc_c$Treatment <- relevel(tc_c$Treatment, ref = "10-Year Bare")
# Add in zeros where they are appropriate
tc_c[is.na(CPM), CPM := 0]

# Sumamrise CAZyme reads at the sample level
tc_c <- tc_c[, .(CPM = sum(CPM)), keyby = .(CAZyme, Sample, Treatment)]

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

# Remove CAZy families whicha re fully zero
tc_c[CAZyme %in% c("AA10", "AA10|CBM73", "AA6", "GH5_8"),]

tc_c <- tc_c[!CAZyme %in% c("AA10", "AA10|CBM73", "AA6", "GH5_8"),]


#### End of old code ####

# Import indval results

inds <- fread("CAZy_2/outputData/Blackout_CAZyme_indvals.csv")
inds <- inds[, .(Row.names, group_name, indval, pvalue)]

setnames(inds
         , c("CAZyme", "trt", "indval", "pvalue"))
inds[, CAZyme := gsub("\\.", "|", CAZyme)]
inds[trt == 'bare_old', Treatment := '10-Year\nBare']
inds[trt == 'bare_new', Treatment := '1-Year\nBare']
inds[trt == 'control_new', Treatment := '1-Year\nGrassland']
inds[trt == 'control_old', Treatment := '10-Year\nGrassland']


#### Model indval CAZyme abundance by treatment ####


d_pb <- tc_c[CAZyme %in% inds[#trt %in% c("bare_old", "bare_new")
                              , CAZyme],]

d_pb[CAZyme %in% active_families_list$Cellulases, Substrate := "Cellulose"]
d_pb[CAZyme %in% active_families_list$Xylanases, Substrate := "Xylan"]
d_pb[CAZyme %in% active_families_list$`Auxiliary Activities`, Substrate := "Lignin"]
d_pb[, Substrate := factor(Substrate, levels = c("Cellulose", "Xylan", "Lignin"))]

#### Strip colouring / Taxonomic origins of CAZymes ####

# Importing code from below so I can be lazy
# Summarise the samples 
tct <- tc[CAZyme %in% inds[!is.na(pvalue),]$CAZyme # inds[grepl("bare", trt)]$CAZyme, #names(int_genes)[interesting_genes], 
          , lapply(.SD, sum)
          , .SDcols = samplecols
          , keyby = .(CAZyme, Taxonomy, Domain, Phylum, Class, Order, Family, Genus, Species)]


#tct <- tct[complete.cases(tct),]

tctl <- melt(tct, id = c("CAZyme", "Taxonomy", taxcols)
             , variable.name = "Sample", value.name = "CPM")
tctl[, Treatment := "10-Year Bare"]
tctl[grep("b[4-7]", Sample), Treatment := "1-Year Bare"]
tctl[grep("c[1-3]", Sample), Treatment := "10-Year Grassland"]
tctl[grep("c[4-7]", Sample), Treatment := "1-Year Grassland"]
tctl[, Treatment := factor(Treatment, levels = c("10-Year Bare"
                                                 , "1-Year Bare"
                                                 , "1-Year Grassland"
                                                 , "10-Year Grassland"))]

#tct <- tct[complete.cases(tct),]

# Plot the taxonomic origins of the CAZymes, as well as their abundance

se <- function(x){x <- na.omit(x);se <- sd(x)/sqrt(length(x));return(se)}

tctls <- tctl[, .(CPM = mean(CPM), CPM_se = se(CPM), CPM_sd = sd(CPM))
              , keyby = .(CAZyme, Treatment, Taxonomy, Phylum)]
tctls[, CPM_s2n := CPM / CPM_sd]
tctls[is.nan(CPM_s2n), CPM_s2n := 0.01]

hist(tctls[, CPM_s2n], breaks = 50);abline(v = 1, col = 'red')
min(tctls[!is.nan(CPM_s2n), CPM_s2n])


tctls[CAZyme %in% active_families_list$Cellulases, Substrate := "Cellulose"]
tctls[CAZyme %in% active_families_list$Xylanases, Substrate := "Xylan"]
tctls[CAZyme %in% active_families_list$`Auxiliary Activities`, Substrate := "Lignin"]
tctls[, Substrate := factor(Substrate, levels = c("Cellulose", "Xylan", "Lignin"))]

tx <- tctls[, lapply(strsplit(Taxonomy, ";"), FUN = function(x){
  data.table(SuperKingdom = x[1], Phylum = x[2], Class = x[3], Order = x[4]
             , Family = x[5], Genus = x[6], Species = x[7])
}) %>% rbindlist(.)]

source("Functions/get_highest_resolution_taxonomy.R")

tctls[, Best_classification := get_highest_resolution_taxonomy(tx)]

# TODO: Include as many CAZymes as possible

strip_CAZyme <- c("CBM2|GH5_1", "CBM4|GH9", "GH45", "GH5", "GH5_19", "GH5_25"
                  ,"GH5_35", "GH5_46", "GH6", "CBM60|GH11", "GH30", "GH30_4"
                  , "GH30_2", "AA10|CBM73", "AA12")

# Check we are using the right CAZymes
strip_CAZyme[strip_CAZyme %in% inds[pvalue < 0.05, CAZyme]]
strip_CAZyme[!strip_CAZyme %in% inds[pvalue < 0.05, CAZyme]]
inds[pvalue < 0.05, CAZyme][inds[pvalue < 0.05, CAZyme] %in% strip_CAZyme]
inds[pvalue < 0.05, CAZyme][!inds[pvalue < 0.05, CAZyme] %in% strip_CAZyme]

strip_trt <- inds[match(strip_CAZyme, CAZyme), .(CAZyme, Treatment)][, Treatment]

inds[, Treatment := factor(Treatment
                           , levels = c("10-Year\nBare"
                                        , "1-Year\nBare"
                                        , "1-Year\nGrassland"
                                        , "10-Year\nGrassland"))]
#inds$Treatment
inds_ordered <- inds[CAZyme %in% strip_CAZyme, .SD, keyby = .(Treatment, CAZyme)]

# Order nicely
CAZ_ordered <- data.table(trt = inds_ordered$Treatment
                          , CAZyme = c(c(bare_old = "AA12"
                                         , bare_old = "CBM2|GH5_1"
                                         , bare_old = "GH5_19"
                                         , bare_old = "GH5_25"
                                         , bare_old = "GH5_46"
                                         , bare_old = "GH30"
                                         )
                                       , c(bare_new = 'GH5'
                                           , bare_new = 'GH6'
                                           , bare_new = 'CBM4|GH9' 
                                           , bare_new = 'GH30_2'
                                           , bare_old = "GH30_4")
                                       , c(control_new = 'CBM60|GH11'
                                           , control_new = 'GH45') 
                                       , c(control_old = 'GH5_35')
                          )
                          , Palette = BlackoutPalette[inds_ordered$Treatment]
)

tctls <- merge.data.table(tctls, inds_ordered[, .(CAZyme, Indicator = Treatment)]
                          , by = 'CAZyme')

CAZ_ordered_2 <- unique(tctls[, .(CAZyme, Indicator, Substrate)])
setkeyv(CAZ_ordered_2, c("Substrate", "Indicator", "CAZyme"))
CAZ_ordered_2[, Palette := BlackoutPalette[CAZ_ordered_2$Indicator]]
CAZ_ordered_2
tctls[, CAZyme := factor(CAZyme, levels = CAZ_ordered_2$CAZyme)]


# Find missing CAZy families
unique(d_pb$CAZyme) 
inds_ordered$CAZyme
d_pb <- d_pb[CAZyme %in% inds_ordered$CAZyme, ]

# Finally we can order the CAZymes nicely
d_pb[, CAZyme := factor(CAZyme, levels = CAZ_ordered_2$CAZyme)]

ridiculous_strips_bar <- ggh4x::strip_nested(
  # Horizontal strips 
  background_x = ggh4x::elem_list_rect(fill = 
                                         c(# Top layer
                                           dplyr::case_match(unique(as.character(tctls$Substrate))
                                                             , 'Lignin' ~ 'grey80'
                                                             , 'Cellulose' ~ 'grey80'
                                                             , 'Xylan' ~ 'grey80')
                                           , "grey80"
                                           # Bottom layer
                                           , dplyr::case_match(CAZ_ordered_2$CAZyme,
                                                               CAZ_ordered_2[Palette == '#7b3294', CAZyme] ~ "#7b3294"
                                                               , CAZ_ordered_2[Palette == '#c2a5cf', CAZyme] ~ "#c2a5cf"
                                                               , CAZ_ordered_2[Palette == '#a6dba0', CAZyme] ~ "#a6dba0"
                                                               , CAZ_ordered_2[Palette == '#008837', CAZyme] ~ "#008837")

                                         )
  )
  , text_x = ggh4x::elem_list_text(colour = c(rep("black", 4), rep('white', length(CAZ_ordered$CAZyme))))
  , by_layer_y = FALSE#TRUE
)

#### Plotting bare -CAZymes ####

d_pb

pb <- ggplot(d_pb
             , aes(x = Treatment, y = CPM, fill = Treatment)) + 
  stat_summary(fun = 'mean', geom = 'col') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0.3) +
  geom_point(size = 0.2, position = position_jitter(width = 0.2)) +
  scale_colour_manual(breaks = levels(tc_c$Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(tc_c$Treatment), values = BlackoutPalette) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  ggh4x::facet_nested_wrap(~ Substrate + CAZyme
                           , nrow = 2
                           , scales = 'free_y'
                           , strip = ridiculous_strips_bar) +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank()
        , legend.position = 'none') 


pb
# Okay, so 10 gene families are indicators of bare plots

# This points to these gene families being the ones most involved in carbon degradation in plantless soils

#FIXME: Use cfc_full for this plot
# cfc_full
c_sum <- cfc_full[CAZyme %in% inds[pvalue < 0.05,]$CAZyme
                  , lapply(.SD, sum), keyby = CAZyme, .SDcols = samplecols]
c_sum <- melt(c_sum, id = c("CAZyme"), variable.name = "Sample", value.name = "CPM")

c_sum[, Treatment := "10-Year Bare"]
c_sum[grep("b[4-7]", Sample), Treatment := "1-Year Bare"]
c_sum[grep("c[1-3]", Sample), Treatment := "10-Year Grassland"]
c_sum[grep("c[4-7]", Sample), Treatment := "1-Year Grassland"]
c_sum[, Treatment := factor(Treatment, levels = c("10-Year Bare"
                                                 , "1-Year Bare"
                                                 , "1-Year Grassland"
                                                 , "10-Year Grassland"))]

c_sum[CAZyme %in% active_families_list$Cellulases, Substrate := "Cellulose"]
c_sum[CAZyme %in% active_families_list$Xylanases, Substrate := "Xylan"]
c_sum[CAZyme %in% active_families_list$`Auxiliary Activities`, Substrate := "Lignin"]
c_sum[, Substrate := factor(Substrate, levels = c("Cellulose", "Xylan", "Lignin"))]
c_sum <- merge.data.table(c_sum, inds_ordered[, .(CAZyme, Indicator = Treatment)]
                          , by = 'CAZyme')

CAZ_ordered_3 <- unique(c_sum[, .(CAZyme, Indicator, Substrate)])
setkeyv(CAZ_ordered_3, c("Substrate", "Indicator", "CAZyme"))
CAZ_ordered_3[, Palette := BlackoutPalette[CAZ_ordered_3$Indicator]]
CAZ_ordered_3
c_sum[, CAZyme := factor(CAZyme, levels = CAZ_ordered_3$CAZyme)]


ridiculous_strips_bar <- ggh4x::strip_nested(
  # Horizontal strips 
  background_x = ggh4x::elem_list_rect(fill = 
                                         c(# Top layer
                                           dplyr::case_match(unique(as.character(c_sum$Substrate))
                                                             , 'Lignin' ~ 'grey80'
                                                             , 'Cellulose' ~ 'grey80'
                                                             , 'Xylan' ~ 'grey80')
                                           , rep("grey80", 2)
                                           # Bottom layer
                                           , dplyr::case_match(CAZ_ordered_3$CAZyme,
                                                               CAZ_ordered_3[Palette == '#7b3294', CAZyme] ~ "#7b3294"
                                                               , CAZ_ordered_3[Palette == '#c2a5cf', CAZyme] ~ "#c2a5cf"
                                                               , CAZ_ordered_3[Palette == '#a6dba0', CAZyme] ~ "#a6dba0"
                                                               , CAZ_ordered_3[Palette == '#008837', CAZyme] ~ "#008837")
                                           
                                         )
  )
  , text_x = ggh4x::elem_list_text(colour = c(rep("black", 5)
                                              #, rep('white', length(CAZ_ordered_3$CAZyme))
                                              , rep("white", 4)
                                              , rep("black", 4)
                                              , rep("white", 2)
                                              , rep("black", 3)
                                              , rep("white", 1)
                                              )
                                   , size = c(rep(8, 5), rep(6, length(CAZ_ordered_3$CAZyme)))
                                   )
  , by_layer_y = FALSE#TRUE
)

unique(c_sum[, .(CAZyme, Indicator)])

pb2 <- ggplot(c_sum
             , aes(x = Treatment, y = CPM, fill = Treatment)) + 
  stat_summary(fun = 'mean', geom = 'col') +
  stat_summary(fun.data = 'mean_cl_boot', geom = 'errorbar', width = 0.3) +
  geom_point(size = 0.2, position = position_jitter(width = 0.2)) +
  scale_colour_manual(breaks = levels(tc_c$Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(tc_c$Treatment), values = BlackoutPalette) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  ggh4x::facet_nested_wrap(~ Substrate + CAZyme
                           , nrow = 3
                           , scales = 'free_y'
                           , strip = ridiculous_strips_bar) +
  theme(axis.text.x = element_blank()
        , axis.ticks.x = element_blank()
        , axis.title.y = element_text(vjust = 1)
        , legend.position = 'bottom' #'none'
        ) 


pb2

tiff(filename = "Figures_2024/Review_updates/Figure_S7_Gene_abundances.tiff"
     , res = 600
     , units = "in"
     , width = 8.3, height = 8.3
);pb2;dev.off() 



#### Strip editing ####


ridiculous_strips <- ggh4x::strip_nested(
  # Horizontal strips 
  background_x = ggh4x::elem_list_rect(fill = 
                                         c(# Top layer
                                           dplyr::case_match(unique(as.character(tctls$Substrate))
                                                             , 'Lignin' ~ 'grey80'
                                                             , 'Cellulose' ~ 'grey80'
                                                             , 'Xylan' ~ 'grey80')
                                           # Bottom layer
                                           , dplyr::case_match(CAZ_ordered_2$CAZyme,
                                                               CAZ_ordered_2[Palette == '#7b3294', CAZyme] ~ "#7b3294"
                                                               , CAZ_ordered_2[Palette == '#c2a5cf', CAZyme] ~ "#c2a5cf"
                                                               , CAZ_ordered_2[Palette == '#a6dba0', CAZyme] ~ "#a6dba0"
                                                               , CAZ_ordered_2[Palette == '#008837', CAZyme] ~ "#008837")
                                         )
  )
  , text_x = ggh4x::elem_list_text(colour = c(rep("black", 3)
                                              #, rep('white', length(CAZ_ordered$CAZyme))
                                              , rep('white', 4)
                                              , rep('black', 4)
                                              , rep('white', 2)
                                              , rep('black', 3)
                                              , 'white'
                                              )
                                   , size = c(rep(6, 3), rep(6, length(CAZ_ordered$CAZyme)))
                                   )
  , by_layer_y = FALSE#TRUE
)

#### Finding missing CPM ####
# README: This was a contamination issue

# tctls[, .(TC = sum(CPM)), by = .(CAZyme, Treatment)]
# t_summ <- tctl[, .(TC = sum(CPM)), by = .(CAZyme, Sample)]
# c_summ <- c_sum[CAZyme %in% tctl$CAZyme, .(CS = sum(CPM)), by = .(CAZyme, Sample)]
# 
# c_summ <- c_sum[CAZyme %in% tctl$CAZyme, .(CS = sum(CPM)), by = .(CAZyme, Sample)]
# c_summ[CAZyme == "GH30_4",]
# c_sum[CAZyme == "GH30_4",]
# cfc_full[CAZyme == "GH30_4",]
# tc_comp <- merge(t_summ, c_summ, by = c("CAZyme", "Sample"))
# 
# # Read in all contigs
# cont <- fread("C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/TPM_contigs/Blackout_contigs_TPM.csv"
#                , nThread = 3)
# # Keep only CAZy contigs
#cont[Contig %in% 'k141_2024747',] # :) Cool. This one (GH30_4) can stay.
# tc[Contig == 'k141_2024747', (samplecols) := .(b1 = 0.270937, b2 =0, b3=0
#                                                ,b4=0,b5= 0.2189873,b6= 0.3863135,b7= 0.2322574
#                                                ,c1 =0,c2 =0,c3=0
#                                                ,c4=0,c5=0,c6=0,c7=0)]
# tc[Contig == 'k141_2024747',]
# cont <- cont[Contig %in% sub("_[0-9]+$", "", cfc_full$GenestoSearch),]
# 
# cont[!Contig %in% tc$Contig, ]
# tc[!Contig %in% cont$Contig,]
# cont[Contig == "k141_6342750",]
# tc[Contig == "k141_6342750",]
# # 
# # tc[Contig == "k141_6342750"
# #    , (samplecols) := cont[Contig == "k141_6342750", .(b1, b2, b3, b4, b5, b6)]]
# 
# tc[Contig == "k141_6342750", b1 := cont[Contig == "k141_6342750",]$b1]
# tc[Contig == "k141_6342750", b2 := cont[Contig == "k141_6342750",]$b2]
# tc[Contig == "k141_6342750", b3 := cont[Contig == "k141_6342750",]$b3]
# tc[Contig == "k141_6342750", b4 := cont[Contig == "k141_6342750",]$b4]
# tc[Contig == "k141_6342750", b5 := cont[Contig == "k141_6342750",]$b5]
# tc[Contig == "k141_6342750", b6 := cont[Contig == "k141_6342750",]$b6]
# tc[Contig == "k141_6342750", b7 := cont[Contig == "k141_6342750",]$b7]
# tc[Contig == "k141_6342750", c1 := cont[Contig == "k141_6342750",]$c1]
# tc[Contig == "k141_6342750", c2 := cont[Contig == "k141_6342750",]$c2]
# tc[Contig == "k141_6342750", c3 := cont[Contig == "k141_6342750",]$c3]
# tc[Contig == "k141_6342750", c4 := cont[Contig == "k141_6342750",]$c4]
# tc[Contig == "k141_6342750", c5 := cont[Contig == "k141_6342750",]$c5]
# tc[Contig == "k141_6342750", c6 := cont[Contig == "k141_6342750",]$c6]
# tc[Contig == "k141_6342750", c7 := cont[Contig == "k141_6342750",]$c7]
# 
# tc[Contig == "k141_6342750",]
# 
# tctls

#### Plotting ####

tctls[Phylum == 'UNASSIGNED', Phylum := "Unassigned"]

fwrite(tctls, "CAZy_taxonomy_2022/Outfiles/Taxonomy_CAZy_links.csv")

# EDIT 2025/00/22: wrap the CAZy names
tctls[, CAZyme := gsub("\\|", "\n|", CAZyme)]

pb_ts <- ggplot(tctls, aes(x = Treatment, y = Taxonomy
                           , size = CPM
                           , colour = Treatment
)) +
  geom_point(aes(alpha = CPM_s2n)) +
  scale_colour_manual(breaks = levels(tctls$Treatment), values = BlackoutPalette) +
  scale_size_area(max_size = 3) +
  scale_alpha_continuous(name = "CPM Signal to noise"
                         , trans = 'log10') +
  theme_bw() +
  theme(legend.position = 'bottom'
        , strip.text.y = element_text(angle = 0, size = 5)
        , strip.text.x = element_text(size = 5)
        , axis.text.y = element_text(size = 3, face = "italic")
        , axis.text.x = element_blank()
        , axis.ticks.x = element_blank()
  ) + 
  ggh4x::facet_nested(Phylum ~ Substrate + CAZyme
                      , scales = 'free', space = 'free'
                      , strip = ridiculous_strips
  ) +
  scale_y_discrete(name = "", breaks = tctls[, unique(.(Taxonomy, Best_classification))]$V1
                   , labels = tctls[, unique(.(Taxonomy = Taxonomy
                                               , Best_classification = Best_classification))]$V2 #unique(tctls$Taxonomy) %>% gsub('.*;', "", .)
  )

pb_ts




#### Create final versions of the plots to combine ####

pb_f <- pb2 + theme(axis.title.y = element_text(vjust = -15
                                                )
                   , title = element_text(size = 10)
                   , axis.text = element_text(size = 6)
                   )# +

# Removed "species
pb_ts_f <- pb_ts + theme(panel.grid = element_line(linewidth = 0.01)
                         , strip.background = element_rect(colour = 'black'
                                                           , fill = 'grey80')
                         , legend.position = 'bottom'
                         , legend.box = "vertical"
                         , legend.margin = margin()
                         , strip.text.y = element_text(angle = 0, size = 6) # 6 is okay
                         , strip.text.x = element_text(size = 10) # 5,6,7 too small -- not working
                         , axis.text.y = element_text(size = 6, face = "italic") # 4, 5 too small
                         , axis.text.x = element_blank()
                         , axis.ticks.x = element_blank()
                         , plot.margin = unit(c(0,0,0,0), units = "cm")
                         )

#### Save the figure ####

# top_layer <- ((pb_f|pg_f) + 
#                 patchwork::plot_layout(widths = c(5,3)
#                                        #, heights = c(1,3.5,0.8)
#                 )
#               
# )
# 
# (pb_f/(pb_ts_f)/guide_area() + 
#     patchwork::plot_layout(guides = 'collect'
#                            #, widths = c(5,3,8,8)
#                            , heights = c(1,3.5,0.8))
# ) +
#   plot_annotation(tag_levels = 'A')

pb_ts_f

# tiff(filename = "Figures_2024/Review_updates/Figure_4_Gene_abundances_and_origins.tiff"
#      , res = 600
#      , units = "in"
#      , width = 8.3, height = 11.7
#      );pb_ts_f;dev.off() 
tiff(filename = "Figures_environ_microb/Figure_4_Gene_abundances_and_origins.tiff"
     , res = 600
     , units = "in"
     , width = 8.3, height = 11.7
);pb_ts_f;dev.off()

pb_f
#|pg_f


#### Summarise the results ####

tctls

tctls[Indicator == '10-Year\nBare' & CPM != 0, unique(Taxonomy)] %>% length(.)
tctls[Indicator == '10-Year\nGrassland' & CPM != 0, unique(Taxonomy)] %>% length(.)
tctls[Indicator == '1-Year\nGrassland' & CPM != 0, unique(Taxonomy)] %>% length(.)
tctls[Indicator == '1-Year\nBare' & CPM != 0, unique(Taxonomy)] %>% length(.)

tctls[Indicator == '10-Year\nBare' & CPM != 0, unique(CAZyme)] %>% length(.)




