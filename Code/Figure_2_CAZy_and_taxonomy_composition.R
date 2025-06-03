#### Metabolome ####

# Metabolite NMDS

# Set things up
rm(list = ls())

library(data.table)
library(magrittr)
library(vegan)
library(ggplot2)
#library(ggpubr)
source("Functions/stat_chull.R")

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

sample_cols <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7'
                 , 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7')

Treatment <- factor(c(rep("b_old", 3), rep("b_new", 4)
                      , rep("c_old", 3), rep("c_new", 4))
                    , levels = c("b_old", "b_new", "c_new", "c_old"))

##### (1) Import the data #####
mf <- read.csv("Metabolomics/Data/Metabolome_Blackout_Raw.csv"); setDT(mf)
KEGG <- fread("Metabolomics/Cleaned Data/Henfaes_Metabolome_Aggregated_KEGG_V2.csv")

# Get total metabolites for scaling later
sample_totals <- mf[, lapply(.SD, sum), .SDcols = sample_cols] %>% 
  transpose(., keep.names = "Sample")
setnames(sample_totals, old = "V1", new = "Total_metabolites")

# Name unknown compounds something short
mf[grep("^[0-9]+$", BinBase.name), BinBase.name := paste0(".", .I)]

# Get the compound matrix
cols_to_keep <- c("BinBase.name", sample_cols)
mw <- transpose(mf[, ..cols_to_keep], make.names = "BinBase.name"
                #, keep.names = "Sample"
)
# Scale by sample totals
mwr <- mw/sample_totals$Total_metabolites

named_compounds <- mf$BinBase.name[!grepl("^\\.", mf$BinBase.name)]

# Do an NMDS
nmds1 <- metaMDS(log(mwr+1)#decostand(mwr, method = 'log')
)

# plot the NMDS
plot(nmds1, type = 'n')
text(vegan::scores(nmds1, display = 'species')
     , labels = row.names(vegan::scores(nmds1, display = 'species'))
     , col = "grey50")
points(vegan::scores(nmds1, display = 'sites')
       , col = BlackoutPalette[Treatment]
       , pch = 16
)

# TODO: Add large differences in chemical composition using envfit


nmds_obj_metabolome <- data.table(sample_totals, Treatment)
nmds_obj_metabolome[, NMDS1 := vegan::scores(nmds1, display = "sites", choices = 1)]
nmds_obj_metabolome[, NMDS2 := vegan::scores(nmds1, display = "sites", choices = 2)]

Metabolome_comp <- ggplot(nmds_obj_metabolome
                         , aes(x = NMDS1
                               , y = NMDS2
                               , colour = Treatment
                               , fill = Treatment
                               , group = Treatment)) +
  stat_chull(alpha = 0.7, colour = NA) +
  geom_point(size = 3) +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  theme_classic() +
  theme(text = element_text(size = 16)
        , legend.position = "bottom")# +

Metabolome_comp

#### Envfit of aggregated chemicals ####

KEGG 
names(KEGG)


# C00212__Adenosine -> nucleosides

KEGG2 <- copy(KEGG)

setnames(KEGG2
         , old = c("SL01_Acylaminosugars"#, "FA03_Eicosanoids"
                   , "FA05_Fatty_alcohols"
                   , "FA08_Fatty_amides", "FA01_Fatty_Acids_and_Conjugates"
                   , "PK12_Flavonoids")
         , new = c("Acylaminosugars"#, "Eicosanoids"
                   , "Fatty_alcohols"
                   , "Fatty_amides", "Fatty_Acids_and_Conjugates"
                   , "Flavonoids")
)

#KEGG2[, Nucleosides := Nucleosides + C00212__Adenosine]
KEGG_cols <- c("Carboxylic_acids", "Amino_acids", "Acylaminosugars", 
               "Vitamins", "Oligosaccharides", "Bases", "Monosaccharides", 
               "Nucleosides", "Rhodopsin_family", "Fatty_acids", 
               "Fatty_Acids_and_Conjugates", "Alkaloids_derived_from_nicotinic_acid", 
               "Alkaloids_derived_from_lysine", "Amines",
               "Eicosanoids"#, "Eicosanoids"
               , "Monolignols", 
               "Fatty_alcohols", "Fats", "24_Carbon_atoms", 
               "Isoflavonoids", "Flavonoids", "Fatty_amides"
)

KEGG2[, ..KEGG_cols]

# Make data proportional and then scale
KEGG2[, (KEGG_cols) := lapply(.SD, FUN = function(x){x/Total_metabolites})
      , .SDcols = KEGG_cols]
KEGG2[, (KEGG_cols) := lapply(.SD, scale)
      , .SDcols = KEGG_cols]

#TODO scale

fit <- envfit(nmds1 ~ 
                Carboxylic_acids                  
              + Acylaminosugars #? # Saccharolipids group 1
              + Bases                                   
              #+ C00212__Adenosine      #?       Nucleic acids - nucleosides - should be a double-count of nucleosides         
              + Fatty_acids                          
              #+ Others                                
              # + Oxidoreductases__EC1_   #?     enzyme?           
              # + DG02883__SLC22A7_inhibitor    #?     Human gene     
              + Fats                                      
              #+ Isoflavonoids     # Same as flavanoids                       
              + Oligosaccharides                   
              #+ V1                                
              + Rhodopsin_family                
              + Alkaloids_derived_from_nicotinic_acid # Possibly aggregate
              + Amines                              
              + Eicosanoids     #?     Fatty acids        
              + Fatty_alcohols    #?      Fatty alcohols      
              # + ST04_Bile_acids_and_derivatives     
              + Fatty_amides #? Fastty amides
              + Amino_acids                    
              + Vitamins                       
              + Monosaccharides                
              + Nucleosides                  
              + Fatty_Acids_and_Conjugates
              + Alkaloids_derived_from_lysine  # Possibly_aggregate
              #+ Eicosanoids                     
              + Monolignols                    
              + `24_Carbon_atoms`  #?          
              + Flavonoids #?
              , data = KEGG2
              
)

fit

# # Just check that the breakdown products aren't correlated
# PermMANOVA3 <- adonis2(species_matrix_rel ~ 
#                          Percentage_Carbon 
#                        #+ Percentage_Nitrogen # (2) 0.49
#                        # + PC # (1) 0.46
#                        + NC # 
#                        # + Phosphorous # (3) 0.55
#                        + Total_cations #
#                        # + Hemicellulose_breakdown_products (6) 0.38
#                        # + Lignin_breakdown_products (5) 0.21
#                        #+ glucose (4) 0.136
#                        , data = mf
#                        , by = "margin" # Assess marginal significance of terms
#                        , parallel = 3 # Make it faster
#                        , permutations = 60000)

#PermMANOVA3 # Nope. The original was best. Include these in the envfit. 

# MAke these able to be used with ggplot
en_coord_cont_Metabolome = as.data.frame(vegan::scores(fit, "vectors")) * ordiArrowMul(fit) #* 0.8
en_coord_cont_Metabolome$Treatment <- "b_old"#"10-Year\nBare"
setDT(en_coord_cont_Metabolome)
en_coord_cont_Metabolome[, Variable := row.names(vegan::scores(fit, "vectors"))
              #c("Carbon", "N:C", "Total cations", "Hemicellulose breakdown\nproducts"
              # , "Lignin breakdown\nproducts", "Glucose")
]
en_coord_cont_Metabolome[, p_value := fit$vectors$pvals]
en_coord_cont_Metabolome[, p_value_adj := p.adjust(p_value, method = "fdr")]
en_coord_cont_Metabolome[, Significance := "N.S.",]
en_coord_cont_Metabolome[p_value_adj < 0.1, Significance := ifelse(p_value_adj < 0.05
                                                        , yes = "p < 0.05"
                                                        , no = "p < 0.1"),]
en_coord_cont_Metabolome
en_coord_cont_Metabolome_s <- en_coord_cont_Metabolome[p_value_adj < 0.1,]

# Plot the output
Metabolome_comp_e <- Metabolome_comp +
  geom_segment(data = en_coord_cont_Metabolome_s
               , aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2
                     , linetype = Significance)
               , size = 1, alpha = 0.5, colour = "grey30"
               , arrow = arrow(ends = "last", length = unit(0.1, "inches"))
  ) +
  geom_text(data = en_coord_cont_Metabolome_s, aes(x = NMDS1*1.05, y = NMDS2*1.05), 
            label = en_coord_cont_Metabolome_s$Variable
            , colour = "navy", fontface = "bold", size = 3) +
  annotate("text", 0.4, y = -0.2
           , label = paste0("Stress = ", round(nmds1$stress, 2))
           , size = 4) +
  xlim(-0.6, 0.6)
#geom_text(aes(x = -0.2, y = 0.1
#              , label = paste0("Stress = ", round(nmds1$stress, 2)))
#          , colour = "black")

Metabolome_comp_e

# Write about not in plot
# Monolignols
# Monosaccharides
# Amino acids

#### CAZy ####

spc<-read.csv("Taxonomy/Multi-Method/Cleaned_Data/Henfaes_Species_Multi-Method_TPM_1kbpmin.csv")
tax<-spc[,1:8]
spc<-spc[,9:ncol(spc)]
spc<-t(spc)
trt<-c("bare_old","bare_old","bare_old","bare_new",
       "bare_new","bare_new","bare_new",
       "control_old","control_old","control_old",
       "control_new","control_new","control_new","control_new")

meta<-read.csv("Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv",row.names=1)

library(vegan)
mod1<-metaMDS(spc)
plot(mod1,"sites")
ordispider(mod1,trt,label=T)

mod2<-metaMDS(decostand(spc,"total"))
plot(mod2,"sites")
ordispider(mod2,trt,label=T)

mod3<-metaMDS(decostand(spc,"frequency"))
plot(mod3,"sites")
ordispider(mod3,trt,label=T)

library(devtools)
source_gist("https://gist.github.com/robiwangriff/79633a738b128e64226bf58381232da7")


aa <- read.csv("CAZy_2/outputData/Blackout_aa_matrix.csv")
cl <- read.csv("CAZy_2/outputData/Blackout_cellulase_matrix.csv")
xy <- read.csv("CAZy_2/outputData/Blackout_xylanase_matrix.csv")

# AA10.CBM73 was a contaminant contig. It should be removed from the analysis.
# Evidence in these files:
# "OneDrive - Bangor University/Blackout_David_Nov_24/TPM_contigs/Blackout_contigs_TPM.csv"
# "CAZy/Cleaned_Data/Henfaes_CAZy_Family_SP_TPM_of_all_DNA_contigs.csv"
aa <- aa[, names(aa) != "AA10.CBM73"]

enz <- do.call("cbind", list(aa, cl, xy))
enz <- enz[, !duplicated(names(enz))] # Remove duplicated columns

mod5<-metaMDS(enz)
plot(mod5,"sites")
ordispider(mod5,trt,label=T)

inds <- read.csv("CAZy_2/outputData/Blackout_CAZyme_indvals.csv")

#library(devtools)
#source_gist("https://gist.github.com/robiwangriff/79633a738b128e64226bf58381232da7")

#library(labdsv)
#inds<-indies2df(indvals= indval(enz,trt),taxa=data.frame(row.names=names(enz)),return.alltax=TRUE)

#write.csv(inds, "CAZy_2/outputData/Blackout_CAZyme_indvals.csv")

# inds2<-indies2df(indvals = indval(enz, sub("_.*", "", trt))
#                  ,taxa=data.frame(row.names=names(enz))
#                  ,return.alltax=TRUE)
# 
# inds$Row.names %in% inds2$Row.names
# inds2$Row.names %in% inds$Row.names
# 
# inds2$Row.names == inds$Row.names
# inds2$group_name == sub("_.*", "", inds$group_name)
# itest <- data.table(CAZyme = inds$Row.names, full = inds$group_name, simple = inds2$group_name)
# itest[!is.na(full) & !is.na(simple),]
# itest[!is.na(simple),]
# 
# inds3<-indies2df(indvals = indval(enz[trt!='bare_new'], sub("_.*", "", trt[trt!='bare_new']))
#                  ,taxa=data.frame(row.names=names(enz[trt!='bare_new']))
#                  ,return.alltax=TRUE)
# 
# inds3

# Create a color vector that maps each factor level to a color
colors.spc<-c("bare_new"="plum4","bare_old"="purple4","control_new"="springgreen","control_old"="springgreen4")[as.character(inds$group_name)]
colors.spc[is.na(inds$group_name)] <- "gray"

colors.site<-c("bare_new"="plum4","bare_old"="purple4","control_new"="springgreen","control_old"="springgreen4")[trt]


plot(mod5,"species",type="n")
points(mod5,"sites",col=colors.site,pch=18)
ordispider(mod5,trt,label=T)
#points(mod5, "species",col=colors.spc,pch=8,cex=0.7)
text(mod5,"species",col=colors.spc,cex=0.7)

library(ggplot2)
Sample <- row.names(spc)
blc <- data.frame(Sample, Treatment = trt)

blc <- cbind(blc, vegan::scores(mod5, display = 'sites', choices = 1:2))

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837'
                     , 'grey40')

spc <- vegan::scores(mod5, display = 'species', choices = 1:2) 
spc <- as.data.frame(spc)
spc$CAZyme <- row.names(spc)
spc$CAZyme <- sub("\\.", "|", spc$CAZyme)

trtlevs <- c("bare_old",    "bare_new",  "control_new",  "control_old" )

CAZy_comp <- ggplot(blc, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_point(size = 3) +
  cowplot::theme_cowplot() + 
  geom_text(data = spc
            , aes(label = CAZyme
                  , colour = names(colors.spc))
            , size = 2#2.5#3
            ) +
  ggpubr::stat_chull(aes(fill = Treatment
                         , group = Treatment)
                     , geom = 'polygon'
                     , alpha = 0.7, colour = NA) +
  scale_colour_manual(breaks = c(trtlevs, 'NA'), values = BlackoutPalette
                      , labels = c("10-Year\nBare", "1-Year\nBare", "1-Year\nGrassland", "10-Year\nGrassland", "")) +
  scale_fill_manual(breaks = c(trtlevs, 'NA'), values = BlackoutPalette
                    , labels = c("10-Year\nBare", "1-Year\nBare", "1-Year\nGrassland", "10-Year\nGrassland", "")) +
  #theme_classic() +
  theme(text = element_text(size = 16)
        , legend.position = "bottom") +
  scale_y_continuous("MDS2") 

CAZy_comp

# Do the envfit

mf <- fread("Soil_data/Soil_properties.csv", drop = "Notes")
names(mf)
mcols <- names(mf)[!names(mf) %in% "Sample"]

mf[, NC := Percentage_Nitrogen / Percentage_Carbon]
mf[, PC := Phosphorous / Percentage_Carbon]

####
#### 2022-05-18
####

me <- fread("Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv"
            , select = c("Sample", "xylose", "vanillic.acid"
                         , "glucose"
                         , "fucose"
                         , "benzoic.acid" 
                         , "X4.hydroxybenzoic.acid" 
                         , "X3.6.anhydro.D.galactose"
            ))

setnames(me
         , old = c("X3.6.anhydro.D.galactose", "vanillic.acid", "X4.hydroxybenzoic.acid", "benzoic.acid" )
         , new = c('3,6-anhydro-D-galactose', 'vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid')
)

me

mf <- merge(mf, me, by = "Sample")
# scale everything

mcols <- names(mf)[!names(mf) %in% c("Sample", "Treatment")]
mf[, (mcols) := lapply(.SD, scale), .SDcols = mcols]
#mf <- merge(dt, mf, by = "Sample")

fit <- envfit(mod5 ~ #`3,6-anhydro-D-galactose` + 
                #`4-hydroxybenzoic acid` 
                + xylose #+ `vanillic acid`# + glucose 
              + fucose #+ `benzoic acid`
              + Percentage_Carbon 
              + Percentage_Nitrogen
              + Total_cations 
              #+ Ammonium
              + NC
              
              , data = mf
)

fit
plot(mod5, type="p", display = "sites")
plot(fit, p.max = 0.05)

p.adjust(fit$vectors$pvals, method = "holm") # p-values for metabolites < 0.1
p.adjust(fit$vectors$pvals, method = "BY") # p-values for metabolites < 0.1
p.adjust(fit$vectors$pvals, method = "fdr") # p-values still significant

####
####
####

# Basic plots of community composition

# Get envfit arrow positions
en_coord_cont_CAZy = as.data.frame(vegan::scores(fit, "vectors")) * ordiArrowMul(fit) * 5
en_coord_cont_CAZy$Treatment <- "10-Year\nBare"

row.names(en_coord_cont_CAZy) <- sub("Percentage_", "", row.names(en_coord_cont_CAZy))
row.names(en_coord_cont_CAZy) <- sub("_", " ", row.names(en_coord_cont_CAZy))
library(stringr)
row.names(en_coord_cont_CAZy) <- stringr::str_to_sentence(row.names(en_coord_cont_CAZy))
row.names(en_coord_cont_CAZy)[row.names(en_coord_cont_CAZy) == "Nc"] <- "N:C"
en_coord_cont_CAZy$pval <- p.adjust(fit$vectors$pvals, method = "fdr") 
en_coord_cont_CAZy$Linetype <- ifelse(en_coord_cont_CAZy$pval < 0.05, yes = 'p < 0.05', no = 'p < 0.1')

CAZy_comp <- CAZy_comp +
  #geom_segment(data = en_coord_cont_CAZy
  #             , aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2
  #                   , linetype = Linetype)
  #             , size = 1, alpha = 0.5, colour = "grey30"
  #             , arrow = arrow(ends = "last", length = unit(0.1, "inches"))
  #) +
  #geom_text(data = en_coord_cont_CAZy, aes(x = NMDS1, y = NMDS2), 
  #          label = row.names(en_coord_cont_CAZy)
  #          , colour = "navy", fontface = "bold") +
  annotate(geom = 'text', x = -0.35, y = 0.35
           , label = paste0("Stress = ", round(mod5$stress, 2))
           , colour = "black")

CAZy_comp

#### Taxa ####

# Import useful packages
library(data.table)
library(magrittr)
library(FSA)
library(vegan)
library(multcompView)
library(ggplot2)
library(scales)
library(ggh4x)

# Import the taxonomy data
df <- fread("Taxonomy/Multi-Method/Cleaned_Data/Henfaes_Species_Multi-Method_TPM_1kbpmin.csv")

samplecols <- c("b1", "b2", "b3"
                , "b4", "b5", "b6", "b7"
                , "c1", "c2", "c3"
                , "c4", "c5", "c6", "c7")

Sample <- samplecols
Treatment <- factor(
  c(paste(substr(Sample, 1, 1)
          , ifelse(substr(Sample, 2, 2) < 4
                   , yes = "Old", no = "New")
          , sep = "-")))
levels(Treatment) <- c("1-Year\nBare"
                       , "10-Year\nBare"
                       , "1-Year\nGrassland"
                       , "10-Year\nGrassland")
Treatment <- relevel(Treatment, ref = "10-Year\nBare")


# Transpose the matrix for analysis
colstokeep <- c(samplecols, "Resolved_taxonomy")
df2 <- df[! SuperKingdom %in% c("Viruses", "UNASSIGNED"),]#[! SuperKingdom %in% c("Viruses"),]#

unc <- df[ SuperKingdom %in% c("Viruses", "UNASSIGNED"), lapply(.SD, sum)
           , keyby = SuperKingdom
           , .SDcols = samplecols
][, .SD, .SDcols = samplecols, keyby = SuperKingdom] %>% 
  melt(., id = "SuperKingdom")
unc[, Treatment := rep(Treatment, each = 2)]

#### Tidy species names ####
# Get the highest taxonomic classification for all 
# taxonomic units

source("Functions/get_highest_resolution_taxonomy.R")
df2[, Resolved_taxonomy := get_highest_resolution_taxonomy(.SD)]

df2

species_matrix <- transpose(df2[, ..colstokeep], make.names = "Resolved_taxonomy")

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

dt <- data.table(Sample = Sample, Treatment = Treatment)

abundant_species <- names(species_matrix)[log10(colMeans(species_matrix)) >= 0]

# Scale the species abundance data to help the NMDS
# Standardise does not work with bray curtis.
# Used frequency which divides each value by the total of that
# column (% abundance for that species), then multiplies by N non-zero
# cases. The result is that the mean of the species column has a value of 1

# Changed to total
species_matrix_rel <- decostand(species_matrix#[, ..abundant_species]
                                , method = "frequency" #"total"
)

nmds1 <- metaMDS(species_matrix_rel)
# Shepards test/goodness of fit
goodness(nmds1) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds1) # Produces a Shepards diagram
nmds_obj_taxa <- data.table(Sample
                            , Treatment)  #data.table(nmds1$points)
nmds_obj_taxa$MDS1 <- vegan::scores(nmds1, display = 'sites', choices = 1)
nmds_obj_taxa$MDS2 <- vegan::scores(nmds1, display = 'sites', choices = 2)

summary_table_rel <- as.data.table(species_matrix_rel)

all_tax_cols <- c(c('SuperKingdom', 'Phylum', 'Class', 'Order', 'Family'
                    , 'Genus', 'Species'), 'Resolved_taxonomy')
st_r <- transpose(summary_table_rel, keep.names = "Resolved_taxonomy") %>%
  merge(., df2[, ..all_tax_cols], by = 'Resolved_taxonomy')


summarise_at_rank <- function(data, rank = "Genus"){
  tax_rank <- c('SuperKingdom', 'Phylum', 'Class', 'Order', 'Family'
                , 'Genus', 'Species')
  tax_rank <- tax_rank[1:which(tax_rank == rank)]
  
  st_rs <- data[, lapply(.SD, sum), keyby = tax_rank
                , .SDcols = names(data)[grep("^V", names(data))]]
  st_rs[, apply(.SD, 1, paste, collapse = ';'), .SDcols = tax_rank]
  st_rs[, Resolved_taxonomy := get_highest_resolution_taxonomy(.SD
                                                               , tax_cols = tax_rank)
  ]
  cols_to_keep <- c("Resolved_taxonomy", names(st_rs)[grep("^V", names(st_rs))])
  output_matrix <- transpose(st_rs[, ..cols_to_keep], make.names = "Resolved_taxonomy")
  return(output_matrix)
}

st_rs <- summarise_at_rank(st_r, "Species")
st_rg <- summarise_at_rank(st_r, "Genus")
st_rf <- summarise_at_rank(st_r, "Family")
st_ro <- summarise_at_rank(st_r, "Order")
st_rc <- summarise_at_rank(st_r, "Class")
st_rp <- summarise_at_rank(st_r, "Phylum")
st_rd <- summarise_at_rank(st_r, "SuperKingdom")

# How many samples were they in?
species_cols <- !grepl("Unassigned", names(st_rs))
freq_s <- apply(st_rs[, ..species_cols]>0, 2, sum)

genus_cols <- !grepl("Unassigned", names(st_rg))
freq_g <- apply(st_rg[, ..genus_cols]>0, 2, sum)

family_cols <- !grepl("Unassigned", names(st_rf))
freq_f <- apply(st_rf[, ..family_cols]>0, 2, sum)

order_cols <- !grepl("Unassigned", names(st_ro))
freq_o <- apply(st_ro[, ..order_cols]>0, 2, sum)

class_cols <- !grepl("Unassigned", names(st_rc))
freq_c <- apply(st_rc[, ..class_cols]>0, 2, sum)

phylum_cols <- !grepl("Unassigned", names(st_rp))
freq_p <- apply(st_rp[, ..phylum_cols]>0, 2, sum)

domain_cols <- !grepl("Unassigned", names(st_rd))
freq_d <- apply(st_rd[, ..domain_cols]>0, 2, sum)


hist(freq_s);unique(freq_s)
hist(freq_g);unique(freq_g)
hist(freq_f);unique(freq_f)
hist(freq_o);unique(freq_o)

# species_wa <- data.table(vegan::wascores(nmds_obj_taxa, st_rs[, ..species_cols]), keep.rownames = "Taxonomy")
# genus_wa <- data.table(vegan::wascores(nmds_obj_taxa, st_rg[, ..genus_cols]), keep.rownames = "Taxonomy")
# family_wa <- data.table(vegan::wascores(nmds_obj_taxa, st_rf[, ..family_cols]), keep.rownames = "Taxonomy")
# order_wa <- data.table(vegan::wascores(nmds_obj_taxa, st_ro[, ..order_cols]), keep.rownames = "Taxonomy")
# class_wa <- data.table(vegan::wascores(nmds_obj_taxa, st_rc[, ..class_cols]), keep.rownames = "Taxonomy")
# phylum_wa <- data.table(vegan::wascores(nmds_obj_taxa, st_rp[, ..phylum_cols]), keep.rownames = "Taxonomy")
# domain_wa <- data.table(vegan::wascores(nmds_obj_taxa, st_rd[, ..domain_cols]), keep.rownames = "Taxonomy")


get_taxon_cv_by_treat <- function(data, taxon_cols, treatment
                                  , cutoff = 0.4){
  cv <- data[, lapply(.SD, FUN = function(x){sd(x)/mean(x)}), .SDcols = taxon_cols
             , keyby = treatment]
  taxon_cols_n <- names(data[, ..taxon_cols])
  
  # hist(unlist(cv[, ..phylum_cols_n]))
  cv[, lapply(.SD, min), .SDcols = taxon_cols_n] %>% 
    unlist(.) %>% 
    hist(., breaks = 50)
  
  taxon_cols_low_cv <- cv[, lapply(.SD, FUN = function(x){x <= cutoff}), .SDcols = taxon_cols_n
  ][, sapply(.SD, any)]
  
  taxon_cols_low_cv_n <- names(taxon_cols_low_cv[taxon_cols_low_cv])
  return(taxon_cols_low_cv_n)
}

phylum_to_use <- get_taxon_cv_by_treat(st_rp, taxon_cols = phylum_cols, treatment = Treatment
                                       , cutoff = 0.1)
class_to_use <- get_taxon_cv_by_treat(st_rc, taxon_cols = class_cols, treatment = Treatment
                                      , cutoff = 0.2)
order_to_use <- get_taxon_cv_by_treat(st_ro, taxon_cols = order_cols, treatment = Treatment
                                      , cutoff = 0.2)
families_to_use <- get_taxon_cv_by_treat(st_rf, taxon_cols = family_cols, treatment = Treatment
                                         , cutoff = 0.2)
genera_to_use <- get_taxon_cv_by_treat(st_rg, taxon_cols = genus_cols, treatment = Treatment
                                       , cutoff = 0.3)
species_to_use <- get_taxon_cv_by_treat(st_rs, taxon_cols = species_cols, treatment = Treatment
                                        , cutoff = 0.3)


Community_comp_taxa <- ggplot(nmds_obj_taxa
                              , aes(x = MDS1
                                    , y = MDS2
                              )) +
  # geom_text(data = #species_wa[Taxonomy %in% species_to_use,]
  #  #genus_wa[Taxonomy %in% genera_to_use,]
  # family_wa[Taxonomy %in% families_to_use,]
  # #class_wa[Taxonomy %in% class_to_use,]
  # #phylum_wa[Taxonomy %in% phylum_to_use,]
  #         , aes(label = Taxonomy)
  #         , size = 2) +
  # ggrepel::geom_text_repel(data = 
  #                           # species_wa[Taxonomy %in% species_to_use,]
  #                           # genus_wa[Taxonomy %in% genus_to_use,]
  #                           # family_wa[Taxonomy %in% family_to_use,]
#                           # class_wa[Taxonomy %in% class_to_use,]
#                           # phylum_wa[Taxonomy %in% phylum_to_use,]
#                          , aes(label = Taxonomy)
#                          , size = 2
#                          , max.overlaps = 500
#                          , point.size = NA
#                          , segment.size = 0.1
#                          , min.segment.length = 0.1
#                          ) +
stat_chull(aes(colour = Treatment
               , fill = Treatment
               , group = Treatment)
           , alpha = 0.7, colour = NA) +
  geom_point(aes(colour = Treatment
                 , fill = Treatment
                 , group = Treatment)
             , size = 3) +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  theme_classic() +
  theme(text = element_text(size = 16)
        , legend.position = "bottom") 
#geom_text(aes(x = -0.15, y = 0.22, label = stresslabel), colour = "black")

Community_comp_taxa

tdf <- data.frame(Treatment = Treatment)



mf <- fread("Soil_data/Soil_properties.csv", drop = "Notes")
names(mf)
mcols <- names(mf)[!names(mf) %in% "Sample"]
mf[, NC := Percentage_Nitrogen / Percentage_Carbon]
mf[, PC := Phosphorous / Percentage_Carbon]

me <- fread("Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv"
            , select = c("Sample", "xylose", "vanillic.acid"
                         , "glucose"
                         , "fucose"
                         , "benzoic.acid" 
                         , "X4.hydroxybenzoic.acid" 
                         , "X3.6.anhydro.D.galactose"
            ))

setnames(me
         , old = c("X3.6.anhydro.D.galactose", "vanillic.acid", "X4.hydroxybenzoic.acid", "benzoic.acid" )
         , new = c('3,6-anhydro-D-galactose', 'vanillic acid', '4-hydroxybenzoic acid', 'benzoic acid')
)

me

mf <- merge(mf, me, by = "Sample")
mcols <- names(mf)[!names(mf) %in% "Sample"]

# scale everything

mf[, (mcols) := lapply(.SD, scale), .SDcols = mcols]
mf <- merge(dt, mf, by = "Sample")

fit <- envfit(nmds1 ~ #`3,6-anhydro-D-galactose` + 
                #`4-hydroxybenzoic acid` 
                + xylose #+ `vanillic acid`# + glucose 
              + fucose #+ `benzoic acid`
              + Percentage_Carbon 
              + Percentage_Nitrogen
              + Total_cations 
              #+ Ammonium
              + NC
              
              , data = mf)

plot(nmds1, type="p", display = "sites")
plot(fit, p.max = 0.05)

fit_phyla <- envfit(nmds1 ~ .
                    , data = st_rp)

# Basic plots of community composition

# Get envfit arrow positions
en_coord_cont_taxa = as.data.frame(vegan::scores(fit, "vectors")) * ordiArrowMul(fit) * 2.5 #3 #4 #5
en_coord_cont_taxa$Treatment <- "10-Year\nBare"

row.names(en_coord_cont_taxa) <- sub("Percentage_", "", row.names(en_coord_cont_taxa))
row.names(en_coord_cont_taxa) <- sub("_", " ", row.names(en_coord_cont_taxa))
library(stringr)
row.names(en_coord_cont_taxa) <- stringr::str_to_sentence(row.names(en_coord_cont_taxa))
row.names(en_coord_cont_taxa)[row.names(en_coord_cont_taxa) == "Nc"] <- "N:C"

en_coord_cont_taxa$p_value_adj <- p.adjust(fit$vectors$pvals, method = 'fdr')
en_coord_cont_taxa$Linetype <- ifelse(en_coord_cont_taxa$p_value_adj < 0.05, yes = 'p < 0.05', no = 'p < 0.1')

en_coord_cont_taxa <- en_coord_cont_taxa[en_coord_cont_taxa$p_value_adj < 0.1,]

# Phylum
en_coord_cont_phylum = as.data.frame(vegan::scores(fit_phyla, "vectors")) * ordiArrowMul(fit_phyla) * 4
en_coord_cont_phylum$Treatment <- "10-Year\nBare"

row.names(en_coord_cont_phylum) <- sub("Percentage_", "", row.names(en_coord_cont_phylum))
row.names(en_coord_cont_phylum) <- sub("_", " ", row.names(en_coord_cont_phylum))
library(stringr)
row.names(en_coord_cont_phylum) <- stringr::str_to_sentence(row.names(en_coord_cont_phylum))
row.names(en_coord_cont_phylum)[row.names(en_coord_cont_phylum) == "Nc"] <- "N:C"

en_coord_cont_phylum$p_value_adj <- p.adjust(fit_phyla$vectors$pvals, method = 'fdr')
en_coord_cont_phylum$Linetype <- ifelse(en_coord_cont_phylum$p_value_adj < 0.05, yes = 'p < 0.05', no = 'p < 0.1')

en_coord_cont_phylum <- en_coord_cont_phylum[en_coord_cont_phylum$p_value_adj < 0.1,]


Taxa_comp <- Community_comp_taxa +
    geom_text(data = en_coord_cont_phylum[en_coord_cont_phylum$p_value_adj < 0.05,]
            , aes(x = NMDS1, y = NMDS2)
            , label = row.names(en_coord_cont_phylum[en_coord_cont_phylum$p_value_adj < 0.05,])
            , size = 2
            ) +
  geom_segment(data = en_coord_cont_taxa
               , aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2)
               , size = 1, alpha = 0.5, colour = "grey30"
               , arrow = arrow(ends = "last", length = unit(0.1, "inches"))
  ) +
  geom_text(data = en_coord_cont_taxa, aes(x = NMDS1, y = NMDS2), 
            label = row.names(en_coord_cont_taxa)
            , colour = "navy", fontface = "bold") +
  annotate(geom = 'text', x = 0.1 #0.2
           , y = 0.2
           , label = paste0("Stress = ", round(nmds1$stress, 2))
           , size = 4)
#, colour = "black")
  # geom_text(aes(x = -0.12, y = -0.055
  #               , label = paste0("Stress = ", round(nmds1$stress, 2)))
  #           , colour = "black")

Taxa_comp
#### Everything ####

# Taxa_comp
# Metabolome_comp_e
# CAZy_comp

library(patchwork)

nmds_obj_taxa
nmds_obj_metabolome
blc

fig_metabolome <- Metabolome_comp_e + 
  guides(colour = 'none', fill = 'none') +
  theme(legend.text = element_text(size = 6)
        , legend.title = element_text(size = 6.5))
  #theme(legend.position = 'none')


fig_taxa <- Taxa_comp + theme(legend.position = 'none') +
  scale_x_continuous("MDS1"
                     , limits = c(-0.3, 0.25)
  #, expansion(mult = c(0.1, 0.2))
  )
 # xlim(-0.4, 0.4) #+ ylim(-0.4,0.4)
# theme(legend.text = element_text(size = 6)
#                               , legend.title = element_text(size = 6.5))


fig_CAZy <- CAZy_comp + theme(legend.text = element_text(size = 6)
                              , legend.title = element_text(size = 6.5)
                              , legend.direction = 'vertical') +
  scale_x_continuous("MDS1", limits = c(-0.5, 0.55#-0.55, 0.55
  )
  #, expansion(mult = c(0.2, 0.2))
  )
  #xlim(-0.7, 0.7)
  #+ theme(legend.position = 'none')
  

allplots <- (fig_taxa / fig_CAZy) # / fig_metabolome)


tiff("Figures_2024/Figure_2_Taxa_and_Lignocellulase_composition.tiff"
     , units = "in"
     , width = 6, height = 6
     , res = 600)
allplots + #| guide_area() +
  plot_layout(guides = 'collect'
              #, heights = c(1,1, 0.3)
  ) +
  plot_annotation(#title = c("Metagenome taxonomy", "Lignocellulases", "Metabolome"),
    tag_levels = "A") & theme(plot.margin = unit(c(0,0,0,0)
                                                 , units = "in")
    )
dev.off()



tiff("C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/Multi_panel_NMDS/Multipanel_NMDS.tiff"
     , units = "in"
     , width = 6, height = 6
     , res = 600)
 allplots + #| guide_area() +
  plot_layout(guides = 'collect'
              #, heights = c(1,1, 0.3)
              ) +
  plot_annotation(#title = c("Metagenome taxonomy", "Lignocellulases", "Metabolome"),
                  tag_levels = "A") & theme(plot.margin = unit(c(0,0,0,0)
                                                               , units = "in")
                  )
dev.off()



#### Metabolome + Fermentation ####

fig_metabolome
