#### Henfaes Species community composition analysis ####

# Questions to ask are: 
# (1) Did species richness differ between treatments?
# (2) Did species diversity differ between treatments?
# (3) Did beta diversity (NMDS distances) differ between treatments?
# (4) Which taxa changed between treatments, if there was a community-level effect?

#### Data and function import and cleaning ####
rm(list = ls()) # Clean the environment

# Import useful packages
library(data.table)
library(FSA)
library(vegan)
library(multcompView)
library(ggplot2)
library(scales)
library(ggh4x)

source("Functions/univariate_kw_dunn_plot.R")
source("Functions/stat_chull.R")

# Palette
DomainPalette <- c('#e41a1c'# Archaea
  , '#377eb8'# Bacteria
  , '#4daf4a'# Eukaryotes
)

# Import the taxonomy data
df <- fread("Data/Henfaes_Species_Multi-Method_TPM_1kbpmin.csv")
df
sum(df$n1)
samplecols <- c("b1", "b2", "b3"
                , "b4", "b5", "b6", "b7"
                , "c1", "c2", "c3"
                , "c4", "c5", "c6", "c7")

# Transpose the matrix for analysis
colstokeep <- c(samplecols, "Taxonomy")
df2 <- df[! SuperKingdom %in% c("Viruses", "UNASSIGNED"),]
species_matrix <- transpose(df2[, ..colstokeep], make.names = "Taxonomy")

df[, 1:5] # Rows are samples in order

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

# Create a sensible 
mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}

BlackoutPalette <- c('#7b3294'
  ,'#c2a5cf'
  ,'#a6dba0'
  ,'#008837')

dt <- data.table(Sample = Sample, Treatment = Treatment)

#### (1) Did species richness differ between treatments? ####

dt[, Species_richness := specnumber(species_matrix)] 

Species_richness_results <- univariate_kw_dunn_plot("Species_richness"
                                                    , y_lab = "Species richness"
                                                    , Letter_height = 5220
                                                    , p_type = "unadjusted"
)
Species_richness_results$Kruskal_wallis
Species_richness_results$Dunn
Species_richness_results$Plot

#### (2) Did species diversity differ between treatments? ####

dt[, Simpson := diversity(species_matrix, index = "simpson")] 

Simpson_results <- univariate_kw_dunn_plot("Simpson"
                                           , y_lab = "Simpson's D"
                                           , Letter_height = 0.990#1
                                           , p_type = "unadjusted"
)
Simpson_results$Kruskal_wallis
Simpson_results$Dunn
Simpson_results$Plot #+ ylim(0,1)
Simpson_results$Effect

dt[, Shannon := diversity(species_matrix, index = "shannon")] 

Shannon_results <- univariate_kw_dunn_plot("Shannon"
                                           , y_lab = "Shannon's H′"
                                           , Letter_height = 6.35
                                           , p_type = "unadjusted"
)
Shannon_results$Kruskal_wallis
Shannon_results$Dunn
Shannon_results$Plot #+ ylim(0,1)


#### (3) Did beta diversity (NMDS distances) differ between treatments?

# Scale the species abundance data to help the NMDS
# Standardise does not work with bray curtis.
# Used frequency which divides each value by the total of that
# column (% abundance for that species), then multiplies by N non-zero
# cases. The result is that the mean of the species column has a value of 1
species_matrix_rel <- decostand(species_matrix, method = "frequency")

nmds1 <- metaMDS(species_matrix_rel)
# Shepards test/goodness of fit
goodness(nmds1) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds1) # Produces a Shepards diagram
nmds_obj <- data.table(nmds1$points)
names(df)

Community_comp <- ggplot(nmds_obj
                         , aes(x = MDS1
                               , y = MDS2
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
#geom_text(aes(x = -0.15, y = 0.22, label = stresslabel), colour = "black")

Community_comp

### Check for significant differences in  microbial community composition

# Check PermMNOVA assumptions 
Group_disp <- betadisper(vegdist(species_matrix_rel, method = "bray")
                         , group = Treatment
                         , type = "median")
plot(Group_disp, label = FALSE)
anova(Group_disp) # Nonsignificant differences in dispersion between groups 
permutest(Group_disp) 
boxplot(Group_disp)

Group_disp <- betadisper(vegdist(species_matrix_rel, method = "bray")
                         , group = Treatment
                         , type = "median"
                         , bias.adjust = TRUE)
plot(Group_disp, label = FALSE)
anova(Group_disp) # Significant differences between groups
permutest(Group_disp) # Significant differences between groups
boxplot(Group_disp)

# The design is nearly balanced and the dispersion between groups is okay, so 
# PermMANOVA should be okay :)

PermMANOVA <- adonis2(species_matrix_rel ~ Treatment
                      , by = "margin" # Assess marginal significance of terms
                      , parallel = 3 # Make it faster
                      , permutations = 60000)

PermMANOVA

# There was a significant difference in community composition due to treatment


#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#?pairwiseAdonis::pairwise.adonis2
tdf <- data.frame(Treatment = Treatment)

#Pairwise_PermMANOVA <- pairwiseAdonis::pairwise.adonis2(species_matrix_rel ~ Treatment
#                                 , by = "margin" # Assess marginal significance of terms
#                                 , data = tdf
#                                 , parallel = 3 # Make it faster
#                                 , permutations = 60000)

mf <- fread("Data/Soil_properties.csv", drop = "Notes")
names(mf)
mcols <- names(mf)[!names(mf) %in% "Sample"]
mf[, NC := Percentage_Nitrogen / Percentage_Carbon]
mf[, PC := Phosphorous / Percentage_Carbon]

####
#### 2022-05-18
####

me <- fread("Data/Henfaes_Taxonomy_Metabolome.csv"
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

fit <- envfit(nmds1 ~ Percentage_Carbon 
              + Percentage_Nitrogen 
              # + Nitrate (4) p = 0.145
              # + Soil_moisture_DW # (2) p = 0.683
              + Total_cations 
              # + Phosphorous  # (1) p = 0.845
              + Ammonium
              + NC
              # + PC # (3) p = 0.44
              , data = mf)

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
fit

fit
plot(nmds1, type="p", display = "sites")
plot(fit, p.max = 0.05)

p.adjust(fit$vectors$pvals, method = "holm") # p-values for metabolites < 0.1
p.adjust(fit$vectors$pvals, method = "BY") # p-values for metabolites < 0.1
p.adjust(fit$vectors$pvals, method = "fdr") # p-values still significant

####
####
####

# Basic plots of community composition

# Get envfit arrow positions
en_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit)# * 0.3
en_coord_cont$Treatment <- "10-Year\nBare"

row.names(en_coord_cont) <- sub("Percentage_", "", row.names(en_coord_cont))
row.names(en_coord_cont) <- sub("_", " ", row.names(en_coord_cont))
library(stringr)
row.names(en_coord_cont) <- stringr::str_to_sentence(row.names(en_coord_cont))
row.names(en_coord_cont)[row.names(en_coord_cont) == "Nc"] <- "N:C"

Community_comp <- Community_comp +
  geom_segment(data = en_coord_cont
               , aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2)
               , size = 1, alpha = 0.5, colour = "grey30"
                 , arrow = arrow(ends = "last", length = unit(0.1, "inches"))
  ) +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), 
            label = row.names(en_coord_cont)
            , colour = "navy", fontface = "bold") +
  geom_text(aes(x = -0.12, y = -0.055
                , label = paste0("Stress = ", round(nmds1$stress, 2)))
            , colour = "black")

Community_comp

# Save a figure for SBB
tiff("Figure_2_Microbial_community_composition.tif"
     , width = 16.5, height = 17
     , units = "cm"
     , res = 600)
Community_comp
dev.off()

#### (4) Which taxa changed between treatments? ####

### Initial attempts with species were too complicated to interpret. 
### Trying again at the genus/family level, with the idea of focussing 
### on particularly affected groups

samplecols
df_g <- df2
df_g

# Get the taxonomy of the species, to genus level
df_g[, Taxonomy_to_genus := sapply(strsplit(Taxonomy, split = ";")
                                   , FUN = function(x){
                                     tax_to_g <- paste(unlist(x)[1:6], collapse = ";")
                                     return(tax_to_g)
                                   })]

# Sum the abundance of these species within genera, and transpose to make a traditional community abundance matrix
genus_matrix <- transpose(df_g[, lapply(.SD, sum), by = c("Taxonomy_to_genus"), .SDcols = samplecols]
                          , make.names = "Taxonomy_to_genus")
genus_matrix[, c("Sample", "Treatment") := .(samplecols, Treatment)]
dfs_s_mat <- genus_matrix[,!c("Treatment", "Sample")]
dfs_s2 <- genus_matrix
matnames <- names(dfs_s_mat)

# Get the mean abundance of genera in eachg treatment
dfs_mean <- dfs_s2[, lapply(.SD, mean), .SDcols = matnames
                   , by = "Treatment"]

dfs_mean[,1:3]
# Put this in long format
dfs_mean <- melt(dfs_mean, id = "Treatment"
                 , value.name = "CPM", variable.name = "Taxonomy")
dfs_mean

# For each genus, calculate the fold change of the mean abundances, relative to in 10-year grassland  
dfs_diff <- dfs_mean[, .(Fold_change = (.SD$CPM +1e-6)/(.SD[Treatment == "10-Year\nGrassland",]$CPM + 1e-6)
                         , Treatment = Treatment
                         , CPM = CPM)
                     , by = c('Taxonomy')]
# Put this on a log2 scale
dfs_diff[, log2FC := log2(Fold_change)]


# View the log2 fold changes
ggplot(dfs_diff[log2FC > 1 | log2FC < -1,]
       , aes(x = log2FC, y = Taxonomy)) +
  geom_point() +
  theme(axis.text.y = element_blank())


#### Create a function to fit the models
spm <- genus_matrix
spm[, Treatment := factor(Treatment)]
spm[, Treatment := relevel(Treatment, ref = "10-Year\nGrassland")]


fitmod <- function(y){
  mod <- glm(car::logit(y) ~ Treatment, data = spm)
}

taxonomy_mods <- lapply(dfs_s_mat/1e6, fitmod) # Fit a logit glm to the species 
# abundance matrix, using the correct treatment factor as a reference

# Get model p values
pvals <- lapply(seq_along(taxonomy_mods)
                , FUN = function(i){
                  res <- drop1(taxonomy_mods[[i]], test = "F")
                  return(res$`Pr(>F)`[2])
                }
)

names(pvals) <- names(taxonomy_mods)
pvals <- data.table(Taxonomy = names(pvals), p = unlist(pvals))
pvals[, p_BY := p.adjust(p, method = "BY")]
pvals[, p_holm := p.adjust(p, method = "holm")]
pvals[, p_bonferroni := p.adjust(p, method = "bonferroni")]

# How many species abundances were significantly affected by 
# experimental treatment?
pvals[which(p < 0.05),]
pvals[which(p_BY < 0.05),]
pvals[which(p_holm < 0.05),]
pvals[which(p_bonferroni < 0.05),]

#l2fc2 <- data.table(l2FC = l2fcs, Taxonomy = names(taxonomy_mods))

# How were the significant differences distributed?
l2fcp <- merge(pvals, dfs_diff, by = "Taxonomy")
#l2fcp <- merge(pvals, l2fc2, by = "Taxonomy")
hist(l2fcp[p_bonferroni < 0.05, log2FC])
nrow(l2fcp[p_bonferroni < 0.05,])
nrow(l2fcp[p_bonferroni < 0.05 & log2FC > 0,])
nrow(l2fcp[p_bonferroni < 0.05 & log2FC < 0,])
nrow(l2fcp[p_bonferroni < 0.05 & log2FC > 1,])
nrow(l2fcp[p_bonferroni < 0.05 & log2FC < -1,])
nrow(l2fcp)

l2fcp[p_BY < 0.05 & log2FC > 1]
l2fcp[p_BY < 0.05 & log2FC < -1]
nrow(l2fcp[p < 0.05 & log2FC > 1])
nrow(l2fcp[p < 0.05 & log2FC < -1])

# Add metadata to the taxonomy data
#colstokeep <- c("Land_use", "Site", names(taxonomy_mods))
rf <- l2fcp

# Extract taxonomic information
rf[, Species := gsub(".+;", "", rf$Taxonomy)]
rf[, Domain := sub(";.*", "", rf$Taxonomy)]
rf[, Phylum := sub("^[a-zA-z]+;", "", rf$Taxonomy)]
rf[, Phylum := sub(";.*", "", rf$Phylum)]
rf[, Family := sub("^[a-zA-z]+;[a-zA-z]+;[a-zA-z]+;[a-zA-z]+;", "", rf$Taxonomy)]
rf[, Family := sub(";.*", "", rf$Family)]
rf[, Genus := sub("^[a-zA-z]+;[a-zA-z]+;[a-zA-z]+;[a-zA-z]+;[a-zA-z]+;", "", rf$Taxonomy)]
rf[, Genus := sub(";.*", "", rf$Genus)]

rf$Genus

# Plot species abundances
ggplot(rf[p_holm < 0.05 & Species != "UNASSIGNED",]
       , aes(x = CPM, y = Species, colour = Treatment)) +
  geom_point() + 
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  facet_grid(Domain~., scales = "free"
             , space = "free")

# Ensure there are no wayward hyphens in the
# taxonomic ranks (thinking of D-T in particular)
checkRanks <- function(i){
  ranks <- strsplit(rf$Taxonomy[i], split = ";")
  length(ranks[[1]]) == 6
}

# Get the highest taxonomic classification for all 
# taxonomic units
# Ugly code, but works
# the line if(! 6 %in classified_ranks) says start at the genus
getHighestRank <- function(i){
  ranks <- strsplit(rfu$Taxonomy[i], split = ";")
  classified_ranks <- which(ranks[[1]] != "UNASSIGNED")
  if(! 6 %in% classified_ranks){
    highest_rank <- max(classified_ranks)
    final_taxonomy <- paste0("Unclassified ", ranks[[1]][highest_rank])
  } else {
    final_taxonomy <- ranks[[1]][6]#[max(classified_ranks)]
  }
  return(final_taxonomy)
}

# Create a unique species, by land use table
rfu <- unique(rf[,.(Species, Domain, Phylum, Taxonomy#, p_bonferroni
)])

# Any warward semicolons?
hyphenbitches <- which(sapply(seq_along(rfu$Taxonomy), checkRanks) != TRUE)
rfu$Taxonomy[hyphenbitches]

# Get the highest taxonomic rank of all analysable species
rfu[, Highest_classification := sapply(seq_along(Taxonomy), getHighestRank)]
str(l2fcp)

rfuc <- merge(rfu[, .(Taxonomy, Highest_classification)], rf, by = "Taxonomy")

# Plot the species with > +-1 log2 fold change, where bonferroni-adjusted
# p-values (from LMMs) say that there was a significant difference
# rfuc[, Highest_classification := factor]

rfuc2 <- rfuc[, .(Highest_classification
                  , Taxonomy_to_genus = Taxonomy
                  #, Species
                  , Family
                  , Genus
                  , p_bonferroni
                  , p_BY
                  , p_holm
                  , p
                  , l2FC = log2FC
                  , col_order1 = rank(log2FC)
                  , Treatment
)
, keyby = c('Domain', 'Phylum')]
rfuc2[, col_order := 1]
rfuc2[Domain == "Archaea", col_order := col_order1]
rfuc2[Domain == "Bacteria", col_order := col_order1 + max(rfuc2[Domain == "Archaea"]$col_order1)]
rfuc2[Domain == "Eukaryota", col_order := col_order1 + max(rfuc2[Domain == "Bacteria"]$col_order)]

# Order the phyla nicely
rfuc2[, Phylum := factor(Phylum, levels = unique(Phylum))]
rfuc2[p_bonferroni < 0.05 & (l2FC > 1 | l2FC < -1),]#[! grep("Unclassified",Species),]


rfuc2[, Dbl_or_hlv := any((l2FC > 1 | l2FC < -1)), by = Highest_classification]
rfuc2[, sig_BY := any(p_BY < 0.05), by = Highest_classification]
rfuc2[, sig_bonferroni := any(p_bonferroni < 0.05), by = Highest_classification]
rfuc2[, sig_holm := any(p_holm < 0.05), by = Highest_classification]

rfuc2[p < 0.05,]

rfuc2[, Order := sub("^[a-zA-z]+;[a-zA-z]+;[a-zA-z]+;", "", rf$Taxonomy)]
rfuc2[, Order := sub(";.*", "", rfuc2$Order)]
rfuc2[, Class := sub("^[a-zA-z]+;[a-zA-z]+;", "", rf$Taxonomy)]
rfuc2[, Class := sub(";.*", "", rfuc2$Class)]

fr_at <- rfuc2[p < 0.05 & Dbl_or_hlv == TRUE & Treatment != "10-Year\nGrassland",]

fr_at[Family == "UNASSIGNED", Family := "No rank"]
fr_at[Order == "UNASSIGNED", Order := "No rank"]
fr_at[Class == "UNASSIGNED", Class := "No rank"]
fr_at[Family == "Archaea" , Family := "No rank"]
fr_at[Class == "Actinobacteria" , Class := "Actinomycetia"]
fr_at[Family == "Bacteria" , Family := "No rank"]

library(magrittr)

phy_class <- unique(fr_at[, Phylum, by = Class]) %>% setkey(Phylum)
phy_ord <- unique(fr_at[, Phylum, by = c("Phylum", "Class", "Order")]) %>% setkey(Phylum)
phy_fam <- unique(fr_at[, Phylum, by = c("Phylum", "Class", "Order", "Family")]) %>% setkey(Phylum)


strip_cols <- factor(c(
  as.character(unique(fr_at$Phylum))
  , as.character(phy_class$Phylum)
  , as.character(phy_ord$Phylum)
  , as.character(phy_fam$Phylum)
))

library(MetBrewer)

strip_cols <- met.brewer("Signac", length(levels(strip_cols)))[strip_cols]


## Create the plot
# TODO: Italicise
p_sp <- ggplot(fr_at
               , aes(x = l2FC, y = Highest_classification #reorder(Highest_classification, col_order)
                     , colour = Treatment
               )) +
  geom_point() + 
  scale_x_continuous(name = expression("Log"[2]*"(Fold change) Relative to 10-year grassland")) +
  scale_y_discrete(name = "Genus or highest classification") +
  scale_colour_manual(values = BlackoutPalette) +
  facet_nested(Phylum+Class+Order+Family~., scales = "free"
               , space = "free"
               , strip = strip_nested(size = "variable"
                                      , background_y = elem_list_rect(fill = strip_cols)
                                      , text_y = elem_list_text(angle = 0
                                                                , size = 8)
               )
  ) +
  geom_vline(xintercept = 0) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size = 10)
        , axis.text.y = element_text(size = 6)
        , panel.border = element_rect(colour = "black", linewidth = 1.1)
        , legend.text = element_text(size = 8)
        , legend.position = "bottom"
  ) +
  geom_vline(xintercept = c(-1,1)
             , linetype = "dashed")

p_sp

# Save a long figure for supplementary information
ggsave("Figure_S1_Species_changes.pdf"
       , plot = p_sp
       , height = 800 
       , width = 350
       , units = "mm")

# now zoom in on particular groups

Responsive_genera <- with(rfuc2[p < 0.05 & Dbl_or_hlv == TRUE & Treatment != "10-Year\nGrassland",]
                          , 
                          unique(data.table(Genus, Family, Phylum, Domain, Taxonomy_to_genus))
)

Responsive_genera

fwrite(Responsive_genera, file = "Responsive_genera.csv.gz")

#### (5) Zoom in on the taxonomy of particularly responsive phyla ####

rfuc2[, Order := sub("^[a-zA-z]+;[a-zA-z]+;[a-zA-z]+;", "", rf$Taxonomy)]
rfuc2[, Order := sub(";.*", "", rfuc2$Order)]
rfuc2[, Class := sub("^[a-zA-z]+;[a-zA-z]+;", "", rf$Taxonomy)]
rfuc2[, Class := sub(";.*", "", rfuc2$Class)]

fr_at <- rfuc2[p < 0.05 & Dbl_or_hlv == TRUE & Treatment != "10-Year\nGrassland",]
#fr_at_s <- fr_at[Phylum %in% c("Thaumarchaeota", "Firmicutes"),]
#fr_at_s <- fr_at[Phylum %in% c("Thaumarchaeota", "Firmicutes", "Actinobacteria"
#                               , "Proteobacteria"),]
#fr_at_s <- fr_at[grep("Thaumarchaeota|Firmicutes", Taxonomy_to_genus),]

fr_at_s <- fr_at[grep("Thaumarch|Firmicutes|Actinobacteria|Proteobacteria"
                      , Taxonomy_to_genus),]

fr_at_s[Family == "UNASSIGNED", Family := "No rank"]
fr_at_s[Order == "UNASSIGNED", Order := "No rank"]
fr_at_s[Class == "UNASSIGNED", Class := "No rank"]
fr_at_s[Family == "Archaea" , Family := "No rank"]
fr_at_s[Class == "Actinobacteria" , Class := "Actinomycetia"]
fr_at_s[Family == "Bacteria" , Family := "No rank"]

# Taxonomy according to NCBI
#
# https://elifesciences.org/articles/14589#:~:text=Robust%20phylum%2Dlevel%20phylogeny%20of,outgroup%20(Materials%20and%20methods).
# https://reader.elsevier.com/reader/sd/pii/S0723202018303291?token=018AE16BE349DD319B1A3CDBC8A9C7C68F09261FF564CE315E527D465FB1D6A253B823630F7DFE07870B430BFA4C09C6&originRegion=eu-west-1&originCreation=20220511145554
# https://www.nature.com/articles/nbt.4229?ref=https://githubhelp.com

## Firmicutes:
## Class: Bacilli, Clostridia, Erysipelotrichia, Tissierellia
## Order: Bacillales, Lactobacillales | Halanaerobiales, Clostridiales | Erysipelotrichales, Tissierellales
## Family in Lactobacillales: (far from Bacillales) Enterococcaceae, Carnobacteriaceae, No Rank
## Family in Bacillales:  (close to lactobacillales) Staphylococcaceae, Planococcaceae, Bacillaceae, Alicyclobacillaceae (close to Clostridiales)
## Family in Clostridiales: Leave it alone.

## Thaumarchaeota:
# seems close to fine
# https://www.sciencedirect.com/science/article/pii/S1369527411000555

## Actinobacteria
# https://www.frontiersin.org/articles/10.3389/fmicb.2018.02007/full
# https://www.researchgate.net/publication/236265596_Phylogenetic_analyses_of_phylum_Actinobacteria_based_on_whole_genome_sequences/figures
# https://en.wikipedia.org/wiki/Actinomycetia
# https://www.nature.com/articles/ismej2017156/figures/1
# https://en.wikipedia.org/wiki/Micrococcales

## Proteobacteria
# Classes:
# https://www.nature.com/articles/s41396-021-00995-x
# https://www.microbiologyresearch.org/docserver/fulltext/ijsem/63/8/2901_ijs049270.pdf?expires=1652371208&id=id&accname=guest&checksum=12C0FDC824BE2DC750FBA4BC46FBDC53
# alphaproteobacteria orders
# https://elifesciences.org/articles/42535#s2
# Gammaproteobacteria orders
# https://journals.asm.org/doi/10.1128/JB.01480-09#F3
# https://www.nature.com/articles/s41396-020-0588-4



### Correcting wrongly unranked genera

# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?lvl=0&id=165813
fr_at_s[Genus == "Sporanaerobacter", Family := "Tissierellaceae"]

# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1323375
fr_at_s[Genus == "Anoxybacter",] #Family := "Tissierellaceae"]

# https://www.frontiersin.org/articles/10.3389/fmicb.2018.02007/full
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2805586
fr_at_s[Genus == "Lawsonella", Family := "Lawsonellaceae"]

# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1410606
fr_at_s[Genus == "Candidatus Nitrosopelagicus",] #Family := "Tissierellaceae"]

# https://www.kegg.jp/kegg-bin/show_organism?org=tah
fr_at_s[Genus == "Candidatus Nitrosotenuis",] #Family := "Tissierellaceae"]

#https://www.frontiersin.org/articles/10.3389/fmicb.2020.608832/full
# https://www.frontiersin.org/articles/10.3389/fmicb.2018.00028/full#h4
fr_at_s[Genus == "Candidatus Nitroscaldus", Family := "Candidatus Nitrosocaldaceae"]

# ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1884634
fr_at_s[Highest_classification == "Candidatus Nanopelagicus", Family := "Candidatus Nanopelagicaceae"]

# https://link.springer.com/article/10.1007/s10096-022-04413-8
fr_at_s[Highest_classification == "Phytobacter", Family := "Enterobacteriaceae"]


unique(fr_at_s$Class)
# Replacing words with symbols for brevity
fr_at_s[Class == "Epsilonproteobacteria", Class := "ε-Proteobacteria"]
fr_at_s[Class == "Deltaproteobacteria", Class := "δ-Proteobacteria"]
fr_at_s[Class == "Alphaproteobacteria", Class := "α-Proteobacteria"]
fr_at_s[Class == "Gammaproteobacteria", Class := "γ-Proteobacteria"]
Class_levels <- c("Nitrososphaeria"
                  
                  # Actinobacteria
                  , "Actinomycetia"
                  
                  
                  # Firmicutes
                  , "Bacilli" , "Clostridia", "Erysipelotrichia", "Tissierellia"
                  
                  # Proteobacteria
                  , "ε-Proteobacteria", "δ-Proteobacteria", "α-Proteobacteria"
                  , "γ-Proteobacteria"
                  
                  
                  , "No rank") 

unique(fr_at_s$Order)
unique(fr_at_s[Class == "Epsilonproteobacteria",]$Order)
unique(fr_at_s[Class == "Deltaproteobacteria",]$Order)
unique(fr_at_s[Class == "Alphaproteobacteria",]$Order)
unique(fr_at_s[Class == "Gammaproteobacteria",]$Order)

fr_at_s[Class == "Epsilonproteobacteria", Class := "ε-Proteobacteria"]
fr_at_s[Class == "Deltaproteobacteria", Class := "δ-Proteobacteria"]
fr_at_s[Class == "Alphaproteobacteria", Class := "α-Proteobacteria"]
fr_at_s[Class == "Gammaproteobacteria", Class := "γ-Proteobacteria"]
unique(fr_at_s$Class)

Class_levels

Order_levels <- c("Nitrosopumilales", "Candidatus Nitrosocaldales", "Nitrososphaerales"
                  
                  # Actinobacteria
                  , "Frankiales", "Corynebacteriales", "Streptosporangiales"
                  , "Candidatus Nanopelagicales", "Micrococcales"
                  
                  # Firmicutes
                  , "Lactobacillales", "Bacillales", "Halanaerobiales", "Clostridiales"
                  , "Erysipelotrichales", "Tissierellales" 
                  
                  # Proteobacteria
                  , "Campylobacterales"
                  , "Desulfuromonadales"
                  , "Rhodospirillales", "Rhizobiales", "Rickettsiales"
                  , "Nevskiales", "Chromatiales", "Legionellales", "Oceanospirillales"
                  , "Alteromonadales", "Vibrionales", "Enterobacterales"
                  
                  ,"No rank")

Family_levels <- c("Candidatus Nitrosocaldaceae", "Nitrososphaeraceae", "Nitrosopumilaceae"
                   
                   # Actinobacteria
                   , "Frankiaceae", "Lawsonellaceae", "Nocardiopsaceae" 
                   , "Candidatus Nanopelagicaceae"
                   , "Microbacteriaceae", "Cellulomonadaceae", "Dermacoccaceae"
                   
                   # Firmicutes
                   , "Alicyclobacillaceae", "Bacillaceae", "Planococcaceae"
                   , "Staphylococcaceae", "Carnobacteriaceae" , "Enterococcaceae"      
                   , "Lachnospiraceae", "Peptococcaceae", "Peptostreptococcaceae"
                   , "Halanaerobiaceae"     , "Erysipelotrichaceae"
                   , "Tissierellaceae"
                   
                   # Proteobacteria
                   #, unique(fr_at_s[Phylum == "Proteobacteria",]$Family)
                   , unique(fr_at_s[Phylum == "Proteobacteria",]$Family)
                   
                   , "No rank")

fr_at_s[Phylum == "Proteobacteria",]$Family
Family_levels <- Family_levels[!duplicated(Family_levels)]

factor(fr_at_s$Class, levels = Class_levels)
factor(fr_at_s$Order, levels = Order_levels)
factor(fr_at_s$Family, levels = Family_levels)

fr_at_s[, Class := factor(Class, levels = Class_levels)]
fr_at_s[, Order := factor(Order, levels = Order_levels)]
fr_at_s[, Family := factor(Family, levels = Family_levels)]

#fr_at_s[Phylum == "Nitrososphaerota",]
fr_at_s2 <- fr_at_s[Phylum != "Actinobacteria",]
#fr_at_s2[Phylum == "Nitrososphaerota", Family]

library(magrittr)

phy_class <- unique(fr_at_s2[, Phylum, by = Class]) %>% setkey(Phylum)
phy_ord <- unique(fr_at_s2[, Phylum, by = c("Phylum", "Class", "Order")]) %>% setkey(Phylum)
phy_fam <- unique(fr_at_s2[, Phylum, by = c("Phylum", "Class", "Order", "Family")]) %>% setkey(Phylum)


strip_cols <- factor(c(
  as.character(unique(fr_at_s2$Phylum))
  , as.character(phy_class$Phylum)
  , as.character(phy_ord$Phylum)
  , as.character(phy_fam$Phylum)
))

#library(MetBrewer)

#strip_cols <- met.brewer("Signac", length(levels(strip_cols)))[strip_cols]
# X, P, X
strip_cols <- c("#5f9ea0","#81a870","#f6a192")[strip_cols]

## Create the plot
# Italicise genus names
fr_at_s2[grep("Unclassified", Highest_classification), Species_assignment := 0]
fr_at_s2[!grep("Unclassified", Highest_classification), Species_assignment := 1]
fr_at_s2[Species_assignment == 0, Highest_classification := 
           sub("Unclassified ", "Unclassified *", Highest_classification)]
fr_at_s2[Species_assignment == 1, Highest_classification := 
           paste0("*", Highest_classification)]
fr_at_s2[, Highest_classification := paste0(Highest_classification, "*")]
fr_at_s2[, Highest_classification]

library(ggtext)

# Make the plot

p_sp2 <- ggplot(fr_at_s2[l2FC >= 1 | l2FC <= -1,]
                , aes(x = l2FC, y = Highest_classification #reorder(Highest_classification, col_order)
                      , colour = Treatment
                )) +
  geom_point(size = 0.7) + 
  scale_x_continuous(name = expression("Log"[2]*"(Fold change) Relative to 10-year grassland")) +
  scale_y_discrete(name = "Genus or highest classification") +
  scale_colour_manual(values = BlackoutPalette) +
  facet_nested(Phylum+Class+Order+Family~., scales = "free"
               , space = "free"
               , strip = strip_nested(size = "variable"
                                      , background_y = elem_list_rect(fill = strip_cols)
                                      , text_y = elem_list_text(angle = 0
                                                                , size = 5)
               )
  ) +
  geom_vline(xintercept = 0) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size = 9)
        , axis.text.y = ggtext::element_markdown(size  = 7)#element_text(size = 7)
        , axis.text.x = element_text(size = 5)
        , panel.border = element_rect(colour = "black", linewidth = 1.1)
        , legend.text = element_text(size = 8)
        , legend.position = "bottom"
  ) +
  geom_vline(xintercept = c(-1,1)
             , linetype = "dashed")

p_sp2

# Save the output
tiff(filename = 'Figure_3_Genus_fold_changes.tif'
     , width = 15.5, height = 23
     , units = "cm"
     , res = 600)
p_sp2
dev.off()
