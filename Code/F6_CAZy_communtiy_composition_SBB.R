rm(list = ls()) # Clean the environment

#### (0) Setup ####

# Import useful packages
library(data.table)
library(ggplot2)
library(scales)
library(vegan)
source("Functions/stat_chull.R")
source("Functions/univariate_kw_dunn_plot.R");library(FSA)

mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

cf <- fread("Data/Henfaes_CAZy_Family_SP_TPM_of_all_DNA.csv") # Total abundance of each CAZy family in a sample

cf

species_matrix <- cf[,! c("Sample", "CAZymes", "Treatment")]

Sample <- cf$Sample
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

dt <- data.table(Sample = Sample, Treatment = Treatment)

### And environmental parameters

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

# And fibre analysis data
ff <- fread("Data/Fibre_Analysis_Data_Processed_2021-04-29.csv"
            , select = c("Sample", "Cellulose", "Hemicellulose", "Lignin", "Labile"))
ff

mf <- merge(mf, ff, by = "Sample")

# scale everything

mcols <- names(mf)[!names(mf) %in% c("Sample", "Treatment")]
mf[, (mcols) := lapply(.SD, scale), .SDcols = mcols]
mf <- merge(dt, mf, by = "Sample")

# Aggregate the scaled lignocellulose breakdown product types, and scale the outcome
mf[, Hemicellulose_breakdown_products := scale(fucose + xylose + `3,6-anhydro-D-galactose`)]
mf[, Lignin_breakdown_products := (`4-hydroxybenzoic acid` + `vanillic acid` + `benzoic acid`)]

#### (0.1) Keep only the families we care about ####

# Subset the CAZy families being analysed to only those with a high proportion of lignocellulolytic families

# Sum abundance of common cellulase CAZy families
common_cellulases <- c("GH5", "GH6", "GH7", "GH8", "GH9", "GH12"
                       , "GH44", "GH45", "GH48")

# GH8 (11%, 5%), GH10 (96%, 0.5%), GH11 (100%, 0.3%), GH30 (37%, 11%)
common_xylanases <- c("GH10", "GH11", "GH8", "GH30")

# https://www.pnas.org/doi/10.1073/pnas.2008888118#sec-1
LPMOs <- c("AA9", "AA10", "AA11", "AA13", "AA14", "AA15")

Auxiliary_cols <- names(cf)[grep("AA", names(cf))]

CBM_cols <- names(cf)[grep("CBM", names(cf))]

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
                                   names(cf)[grep(r, names(cf))]
                                 }))
                                 names(out) <- NULL 
                                 return(out)
                               })
names(active_families_list) <- c("Cellulases", "Xylanases"
                                 , "LPMOs", "Auxiliary Activities"
                                 , "CBMs"
                                 , "Lignocellulases")

unlist(active_families_list)[anyDuplicated(unlist(active_families_list))]
# GH8 is both a xylanase and a cellulase

cellulases_matrix <- as.data.frame(species_matrix)[, names(species_matrix) %in% unlist(active_families_list$Cellulases)]
xylanase_matrix <- as.data.frame(species_matrix)[, names(species_matrix) %in% unlist(active_families_list$Xylanases)]
aa_matrix <- as.data.frame(species_matrix)[, names(species_matrix) %in% unlist(active_families_list$`Auxiliary Activities`)]


#### (1.1) Abundance, richness and diversity of cellulase genes ####

dt[, Cellulase_abundance := rowSums(cellulases_matrix)] 
Cellulase_abundance_results <- univariate_kw_dunn_plot("Cellulase_abundance"
                                                      , y_lab = "Cellulase CAZy family abundance"
                                                      , Letter_height = 0.07
                                                      , p_type = "unadjusted"
                                                      #, p_type = "adjusted"
)
Cellulase_abundance_results$Kruskal_wallis
Cellulase_abundance_results$Dunn
Cellulase_abundance_results$Plot 
Cellulase_abundance_results$Effect

dt[, Cellulase_richness := specnumber(cellulases_matrix)] 
Cellulase_richness_results <- univariate_kw_dunn_plot("Cellulase_richness"
                                                    , y_lab = "Cellulase CAZy family richness"
                                                    , Letter_height = 30
                                                    , p_type = "unadjusted"
                                                    #, p_type = "adjusted"
                                                    )
Cellulase_richness_results$Kruskal_wallis
Cellulase_richness_results$Dunn
Cellulase_richness_results$Plot 
Cellulase_richness_results$Effect

ggsave("Figures_Supplementary/CAZy_richness_cellulases.pdf"
       , plot = Cellulase_richness_results$Plot)

# 1-year grassland >  10-year bare 
# processing sloughed root cells?

# Simpson
dt[, Cellulase_simpson := diversity(cellulases_matrix, index = "simpson")] 
Cellulase_Simpson_results <- univariate_kw_dunn_plot("Cellulase_simpson"
                                                     , y_lab = "Cellulase diversity (D)"
                                                     , Letter_height = 1
                                                     , p_type = "unadjusted")
Cellulase_Simpson_results$Kruskal_wallis
Cellulase_Simpson_results$Dunn
Cellulase_Simpson_results$Plot
Cellulase_Simpson_results$Effect

# No differesnce

# Shannon
dt[, Cellulase_shannon := diversity(cellulases_matrix, index = "shannon")] 
Cellulase_shannon_results <- univariate_kw_dunn_plot("Cellulase_shannon"
                                                     , y_lab = "Cellulase diversity (H')"
                                                     , Letter_height = 3
                                                     , p_type = "unadjusted")
Cellulase_shannon_results$Kruskal_wallis
Cellulase_shannon_results$Dunn
Cellulase_shannon_results$Plot
Cellulase_shannon_results$Effect

# 1-year grass > 10-year-bare
# processing sloughed root cells?

# Check out how the community of cellulase genes are different
pca1 <- prcomp(cellulases_matrix, scale. = TRUE, center = TRUE)
plot(pca1$x, col = BlackoutPalette[dt$Treatment], pch = 16)
text(pca1$rotation[, 1:2]# * pca1$scale +pca1$center
     , labels = row.names(pca1$rotation[, 1:2])
     )#, col = BlackoutPalette[dt$Treatment], pch = 16)

cel_f <- melt(data.table(cellulases_matrix, Treatment = dt$Treatment), id = "Treatment")


# Check out the community composition
nmds1 <- metaMDS(cellulases_matrix)
# Shepards test/goodness of fit
goodness(nmds1) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds1) # Produces a Shepards diagram
nmds_obj_cel <- data.table(nmds1$points)
nmds_obj_cel_species <- data.table(nmds1$species, CAZyme = rownames(nmds1$species))
range(nmds_obj_cel_species$MDS1)
range(nmds_obj_cel_species$MDS2)

# Plot the community composition of cellulase genes, and how they differ between
# treatments
Cellulase_community_comp <- ggplot(nmds_obj_cel
                         , aes(x = MDS1
                               , y = MDS2
                               , colour = Treatment
                               , fill = Treatment
                               , group = Treatment)) +
  geom_text(data = nmds_obj_cel_species, mapping = aes(MDS1, MDS2, label = CAZyme
                                                   , colour = NULL, fill = NULL, group = NULL)
            , size = 2) + 
  stat_chull(alpha = 0.7, colour = NA) +
  geom_point(size = 3) +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  theme_classic() +
  theme(text = element_text(size = 16)
        , axis.title = element_text(size = 15)
        , legend.position = "bottom") +
  annotate("text", x = 0.2, y = 0.25
            , label = paste0("Stress = ", round(nmds1$stress, 2))
            , size = 4) +
  ggtitle("Cellulases")
  
Cellulase_community_comp


# Order the CAZy genes nicely
source("Functions/Get_CAZy_family_ranks.R")

gene_levels <- get_gene_order(gene_families = as.character(unique(cel_f$variable)), class_of_interest = "GH")
cel_f[, variable := factor(variable, levels = rev(gene_levels))]
# This is an important figure, shows specifically where the differences in 
# cellulases are

# sorting the CAZy families by their GH family and then subfamilies would be helpful!
Cellulase_community_comp_family <- ggplot(cel_f, aes(x = value, y = variable, colour = Treatment)) + 
  stat_summary(fun.data = "mean_se"
               , position = position_dodge(0.1)
               , size = 0.3) + 
  scale_colour_manual(breaks = levels(dt$Treatment), values = BlackoutPalette) +
  xlab("Cellulase abundance (CPM)") +
  ylab("CAZy family") +
  theme_bw() +
  theme(panel.grid = element_blank()) 

Cellulase_community_comp_family

PermMANOVA_cellulase <- adonis2(cellulases_matrix ~ Treatment
                      , by = "margin" # Assess marginal significance of terms
                      , parallel = 3 # Make it faster
                      , permutations = 60000)

PermMANOVA_cellulase

# Environmental drivers of cellulase gene community
# Surprising that cellulose quantity isn't on there, but suppose it 
# may be more to do with quality than quantity 
PermMANOVA_cellulase_env <- adonis2(cellulases_matrix ~ 
                                     Percentage_Carbon 
                                   #+ Percentage_Nitrogen # (1) F = 1.4001, p = 0.221
                                   #+ PC # (3) F = 1.4091, p = 0.213
                                   + NC # 
                                   #+ Phosphorous # (2) F = 0.9485 p = 0.45934
                                   + Total_cations #
                                   # + Cellulose # (4) F = 1.4049, p = 0.203
                                   , data = mf
                                   , by = "margin" # Assess marginal significance of terms
                                   , parallel = 3 # Make it faster
                                   , permutations = 60000)

PermMANOVA_cellulase_env

PermMANOVA_cellulase_glucose <- adonis2(cellulases_matrix ~ glucose
                                         , data = mf
                                         , by = NULL
                                         , parallel = 3 # Make it faster
                                         , permutations = 60000)

PermMANOVA_cellulase_cellulose <- adonis2(cellulases_matrix ~ Cellulose
                                             , data = mf
                                             , by = NULL
                                             , parallel = 3 # Make it faster
                                             , permutations = 60000)

PermMANOVA_cellulase_glucose
PermMANOVA_cellulase_cellulose

#### (1.2) Abundance, richness and diversity of xylanase genes ####

# Abundance
dt[, Xylanlase_abundance := rowSums(xylanase_matrix)] 
Xylanlase_abundance_results <- univariate_kw_dunn_plot("Xylanlase_abundance"
                                                       , y_lab = "Xylanlase CAZy family abundance"
                                                       , Letter_height = 0.07
                                                       , p_type = "unadjusted"
                                                       #, p_type = "adjusted"
)
Xylanlase_abundance_results$Kruskal_wallis
Xylanlase_abundance_results$Dunn
Xylanlase_abundance_results$Plot 
Xylanlase_abundance_results$Effect

# Richness
dt[, Xylanase_richness := specnumber(xylanase_matrix)] 
Xylanase_richness_results <- univariate_kw_dunn_plot(y = "Xylanase_richness"
                                                      , y_lab = "Xylanase CAZy family richness"
                                                      , Letter_height = 11
                                                     , p_type = "unadjusted"
                                                     )
Xylanase_richness_results$Kruskal_wallis
Xylanase_richness_results$Dunn
Xylanase_richness_results$Plot
Xylanase_richness_results$Effect


# Differences
# 1-bare > 10-bare & 10-grass

# Simpson
dt[, Xylanase_simpson := diversity(xylanase_matrix, index = "simpson")] 
Xylanase_Simpson_results <- univariate_kw_dunn_plot("Xylanase_simpson"
                                                    , y_lab = "Xylanase diversity (D)"
                                                    , Letter_height = 1
                                                    , p_type = "unadjusted"
)
Xylanase_Simpson_results$Kruskal_wallis
Xylanase_Simpson_results$Dunn
Xylanase_Simpson_results$Plot # No differences
Xylanase_Simpson_results$Effect

dt[, .(Median = median(Xylanase_simpson), IQR = IQR(Xylanase_simpson))]

# Shannon
dt[, Xylanase_shannon := diversity(xylanase_matrix, index = "shannon")] 
Xylanase_Shannon_results <- univariate_kw_dunn_plot("Xylanase_shannon"
                                                    , y_lab = "Xylanase diversity (H')"
                                                    , Letter_height = 1
                                                    , p_type = "unadjusted")
Xylanase_Shannon_results$Kruskal_wallis
Xylanase_Shannon_results$Dunn
Xylanase_Shannon_results$Plot
Xylanase_Shannon_results$Effect

# Check out the community composition
nmds1 <- metaMDS(xylanase_matrix)
# Shepards test/goodness of fit
goodness(nmds1) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds1) # Produces a Shepards diagram
nmds_obj_xyl <- data.table(nmds1$points)
nmds_obj_xyl_species <- data.table(nmds1$species, CAZyme = rownames(nmds1$species))
range(nmds_obj_xyl_species$MDS1)
range(nmds_obj_xyl_species$MDS2)

# Plot the community composition of xylanase genes, and how they differ between
# treatments
Xylanase_community_comp <- ggplot(nmds_obj_xyl
                                   , aes(x = MDS1
                                         , y = MDS2
                                         , colour = Treatment
                                         , fill = Treatment
                                         , group = Treatment)) +
  geom_text(data = nmds_obj_xyl_species, mapping = aes(MDS1, MDS2, label = CAZyme
                                                   , colour = NULL, fill = NULL, group = NULL)
            , size = 2) + 
  stat_chull(alpha = 0.7, colour = NA) +
  geom_point(size = 3) +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  theme_classic() +
  theme(text = element_text(size = 16)
        , legend.position = "bottom") +
  annotate("text", x = -0.25, y = 0.3
           , label = paste0("Stress = ", round(nmds1$stress, 2))
           , size = 4#6
           ) +
  ggtitle("Xylanases")

Xylanase_community_comp

PermMANOVA_xylanase <- adonis2(xylanase_matrix ~ Treatment
                                , by = "margin" # Assess marginal significance of terms
                                , parallel = 3 # Make it faster
                                , permutations = 60000)

PermMANOVA_xylanase

# Environmental drivers of hemicellulase gene community
# Surprising that hemicellulose quantity isn't on there, but suppose it 
# may be more to do with quality than quantity 
PermMANOVA_xylanase_env <- adonis2(xylanase_matrix ~ 
                         Percentage_Carbon 
                         + Percentage_Nitrogen 
                         # + PC # (2) F(1,13) = -0.2283 p = 0.998
                         + NC # 
                         #+ Phosphorous # (1) F = -0.0400, p = 0.99035
                         + Total_cations #
                         #+ Hemicellulose # # (3) F = 0.879, p = 0.482
                       , data = mf
                       , by = "margin" # Assess marginal significance of terms
                       , parallel = 3 # Make it faster
                       , permutations = 60000)

PermMANOVA_xylanase_env

hist(mf$Hemicellulose_breakdown_products)

# What about individual correlations with breakdown products?
PermMANOVA_xylanase_breakdown <- adonis2(xylanase_matrix ~  Hemicellulose_breakdown_products
                                           , data = mf
                                           , by = NULL
                                           , parallel = 3 # Make it faster
                                           , permutations = 60000)

PermMANOVA_xylanase_fucose <- adonis2(xylanase_matrix ~  fucose
                                           , data = mf
                                           , by = NULL
                                           , parallel = 3 # Make it faster
                                           , permutations = 60000)

PermMANOVA_xylanase_xylose <- adonis2(xylanase_matrix ~  xylose
                                           , data = mf
                                           , by = NULL
                                           , parallel = 3 # Make it faster
                                           , permutations = 60000)


PermMANOVA_xylanase_galactose <- adonis2(xylanase_matrix ~ `3,6-anhydro-D-galactose`
                                       , data = mf
                                       , by = NULL
                                       , parallel = 3 # Make it faster
                                       , permutations = 60000)

PermMANOVA_xylanase_hemicellulose <- adonis2(xylanase_matrix ~ Hemicellulose
                                         , data = mf
                                         , by = NULL
                                         , parallel = 3 # Make it faster
                                         , permutations = 60000)

PermMANOVA_xylanase_breakdown
PermMANOVA_xylanase_fucose
PermMANOVA_xylanase_xylose
PermMANOVA_xylanase_galactose
PermMANOVA_xylanase_hemicellulose


source("Functions/AICc_adonis2.R")

#### (1.3) Richness and diversity of AA genes ####

# Abundance
dt[, AA_abundance := rowSums(aa_matrix)] 
AA_abundance_results <- univariate_kw_dunn_plot("AA_abundance"
                                                       , y_lab = "AA CAZy family abundance"
                                                         , Letter_height = 0.08
                                                       , p_type = "unadjusted"
                                                       #, p_type = "adjusted"
)
AA_abundance_results$Kruskal_wallis
AA_abundance_results$Dunn
AA_abundance_results$Plot 
AA_abundance_results$Effect

# Richness
dt[, AA_richness := specnumber(aa_matrix)] 
AA_richness_results <- univariate_kw_dunn_plot("AA_richness"
                                                     , y_lab = "AA family richness"
                                                     , Letter_height = 9
                                               , p_type = "unadjusted"
                                                     )
AA_richness_results$Kruskal_wallis
AA_richness_results$Dunn
AA_richness_results$Plot
AA_richness_results$Effect

# Differences
# 10-bare < 1-bare & 1-grass

# Simpson
dt[, AA_simpson := diversity(aa_matrix, index = "simpson")] 
AA_Simpson_results <- univariate_kw_dunn_plot("AA_simpson"
                                               , y_lab = "AA diversity (D)"
                                               , Letter_height = 0.8
                                              , p_type = "unadjusted"
)
AA_Simpson_results$Kruskal_wallis
AA_Simpson_results$Dunn
AA_Simpson_results$Plot
AA_richness_results$Effect

# 10 year bare different from: 1 year bare, 1 year grassland
# 10-year grassland different from: nothing
# 1 year bare different from: nothing else

# 10-year bare: a
# 1 year bare: b
# 1-year grassland: b
# 10-year grassland: ab


# Shannon
dt[, AA_shannon := diversity(aa_matrix, index = "shannon")] 
AA_Shannon_results <- univariate_kw_dunn_plot("AA_shannon"
                                              , y_lab = "AA diversity (H')"
                                              , Letter_height = 01.6
                                              , p_type = "unadjusted"
)
AA_Shannon_results$Kruskal_wallis
AA_Shannon_results$Dunn
AA_Shannon_results$Plot
AA_Shannon_results$Effect


# Check out the community composition
nmds1 <- metaMDS(aa_matrix)
# Shepards test/goodness of fit
goodness(nmds1) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds1) # Produces a Shepards diagram
nmds_obj_aa <- data.table(nmds1$points)
nmds_obj_aa_species <- data.table(nmds1$species, CAZyme = rownames(nmds1$species))
range(nmds_obj_aa_species$MDS1)
range(nmds_obj_aa_species$MDS2)

# Plot the community composition of xylanase genes, and how they differ between
# treatments
AA_community_comp <- ggplot(nmds_obj_aa
                                  , aes(x = MDS1
                                        , y = MDS2
                                        , colour = Treatment
                                        , fill = Treatment
                                        , group = Treatment)) +
  geom_text(data = nmds_obj_aa_species, mapping = aes(MDS1, MDS2, label = CAZyme
                                                   , colour = NULL, fill = NULL, group = NULL)
            , size = 2) + 
  stat_chull(alpha = 0.7, colour = NA) +
  geom_point(size = 3) +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  theme_classic() +
  theme(text = element_text(size = 16)
        , legend.position = "bottom") +
  annotate("text", x = -0.25, y = 0.3
           , label = paste0("Stress = ", round(nmds1$stress, 2))
           , size = 4) +
  ggtitle("Auxiliary activities")

AA_community_comp

PermMANOVA_AA <- adonis2(aa_matrix ~ Treatment
                               , by = "margin" # Assess marginal significance of terms
                               , parallel = 3 # Make it faster
                               , permutations = 60000)

PermMANOVA_AA

# Environmental drivers of cellulase gene community
# Surprising that hemicellulose quantity isn't on there, but suppose it 
# may be more to do with quality than quantity 
PermMANOVA_aa_env <- adonis2(aa_matrix ~ 
                                    #  Percentage_Carbon #(1) F = 0.2650, p = 0.8572
                                    # + Percentage_Nitrogen # (4) F = 0.7519, p = 0.46158
                                    # + PC # (5) F = 0.6253,  p =  0.5445
                                    # NC #(6) F = 2.8987, p = 0.0750654
                                    #+ Phosphorous # (3) F = 0.6516, p = 0.524
                                    + Total_cations #
                                    #+ Lignin # (2) F = 0.1145 p = 0.9722
                                    , data = mf
                                    , by = "margin" # Assess marginal significance of terms
                                    , parallel = 3 # Make it faster
                                    , permutations = 60000)

PermMANOVA_aa_env

PermMANOVA_aa_BA <- adonis2(aa_matrix ~ `benzoic acid`
                                        , data = mf
                                        , by = NULL
                                        , parallel = 3 # Make it faster
                                        , permutations = 60000)

PermMANOVA_aa_HBA <- adonis2(aa_matrix ~ `4-hydroxybenzoic acid`
                                , data = mf
                                , by = NULL
                                , parallel = 3 # Make it faster
                                , permutations = 60000)

PermMANOVA_aa_VA <- adonis2(aa_matrix ~ `vanillic acid`
                                , data = mf
                                , by = NULL
                                , parallel = 3 # Make it faster
                                , permutations = 60000)

PermMANOVA_aa_breakdown <- adonis2(aa_matrix ~ Lignin_breakdown_products
                            , data = mf
                            , by = NULL
                            , parallel = 3 # Make it faster
                            , permutations = 60000)

PermMANOVA_aa_lignin <- adonis2(aa_matrix ~ Lignin
                                          , data = mf
                                          , by = NULL
                                          , parallel = 3 # Make it faster
                                          , permutations = 60000)

PermMANOVA_aa_BA 
PermMANOVA_aa_HBA 
PermMANOVA_aa_VA 
PermMANOVA_aa_breakdown 
PermMANOVA_aa_lignin 

#### (2) NMDS of all lignocellulase genes and correlates ####

# Decided not to do the below, but interesing anyhow
## Scale the CAZy abundance data to help the NMDS
## Standardise does not work with bray curtis.
## Used frequency which divides each value by the total of that
## column (% abundance for that species), then multiplies by N non-zero
## cases. The result is that the mean of the species column has a value of 1
#decostand(species_matrix, method = "frequency")

species_matrix_rel <- as.data.frame(species_matrix)[, names(species_matrix) %in% unlist(active_families_list$Lignocellulases)]

# Bray curtis and NMDS mean that these don't need to be scaled

nmds1 <- metaMDS(species_matrix_rel)
# Shepards test/goodness of fit
goodness(nmds1) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds1) # Produces a Shepards diagram
#nmds_obj <- data.table(nmds1$points)
nmds_obj <- data.table(scores(nmds1, display = "sites"))
names(df)

Community_comp <- ggplot(nmds_obj
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
anova(Group_disp) # nodifferences between groups
permutest(Group_disp) # no differences between groups
boxplot(Group_disp)

# The design is nearly balanced and the dispersion between groups is okay, so 
# PermMANOVA should be okay :)

PermMANOVA <- adonis2(species_matrix_rel ~ Treatment
                      , by = "margin" # Assess marginal significance of terms
                      , parallel = 3 # Make it faster
                      , permutations = 60000)

PermMANOVA

#### (3) Which environmental parameters had an effect? ####

fit <- envfit(nmds1 ~ Percentage_Carbon 
              + Percentage_Nitrogen 
              # + `benzoic acid` # (1) p = 0.99
              # + #+ Nitrate # (2) p = 0.96
              # + PC # (3) p = 0.89
              # + glucose # (4) p = 0.69
              # + NC # (5) p = 0.61
              # + `3,6-anhydro-D-galactose` # (6) p = 0.6
              # + Phosphorous # (7) p = 0.48
              # + `vanillic acid` # (8) p = 0.49
              # + Soil_moisture_DW # (9) p = 0.40
              # + Ammonium (10) p = 0.35
              + Total_cations 
              # + `4-hydroxybenzoic acid`  # p = 0.24
              + xylose 
              + fucose 
              
              , data = mf)

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

# Hypothesise that CAZyme community is related to soil nutrients (C, N, P), and 
# amounts of these relative to C, pH (via microbial community composition will 
# also have an effect.

PermMANOVA2 <- adonis2(species_matrix_rel ~ 
                       Percentage_Carbon 
                       #+ Percentage_Nitrogen # (2) F = 0.8103, p = 0.546591
                       #+ PC # (1) F = 0.6930, p = 0.6531 # from PermMANOVA p = 0.35 2022-09-20
                       + NC # () p = 0.21 PermMANOVA
                       # + Phosphorous # (3)  F = 1.162, p = 0.313 PermMANOVA
                       + Total_cations
                       , data = mf
                       , by = "margin" # Assess marginal significance of terms
                       , parallel = 3 # Make it faster
                       , permutations = 60000)

PermMANOVA2 # Include these in the envfit. 

# Plot correlations onto the NMDS plot. Add in variables of interest which 
# likely result from the changed community composition 
fit <- envfit(nmds1 ~ 
                Percentage_Carbon 
              + NC
              + Total_cations
              + Hemicellulose_breakdown_products
              + Lignin_breakdown_products
              + glucose
              , data = mf
              )

fit

# Just check that the breakdown products aren't correlated
PermMANOVA3 <- adonis2(species_matrix_rel ~ 
                         Percentage_Carbon 
                       #+ Percentage_Nitrogen # (2) 0.49
                       # + PC # (1) 0.46
                       + NC # 
                       # + Phosphorous # (3) 0.55
                       + Total_cations #
                       # + Hemicellulose_breakdown_products (6) 0.38
                       # + Lignin_breakdown_products (5) 0.21
                       #+ glucose (4) 0.136
                       , data = mf
                       , by = "margin" # Assess marginal significance of terms
                       , parallel = 3 # Make it faster
                       , permutations = 60000)

PermMANOVA3 # Nope. The original was best. Include these in the envfit. 

# MAke these able to be used with ggplot
en_coord_cont = as.data.frame(scores(fit, "vectors")) * ordiArrowMul(fit) #* 0.8
en_coord_cont$Treatment <- "10-Year\nBare"
setDT(en_coord_cont)
en_coord_cont[, Variable := c("Carbon", "N:C", "Total cations", "Hemicellulose breakdown\nproducts"
                              , "Lignin breakdown\nproducts", "Glucose")]
en_coord_cont

# Plot the output
Community_comp <- Community_comp +
  geom_segment(data = en_coord_cont
               , aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2)
               , size = 1, alpha = 0.5, colour = "grey30"
               , arrow = arrow(ends = "last", length = unit(0.1, "inches"))
  ) +
  geom_text(data = en_coord_cont, aes(x = NMDS1*1.05, y = NMDS2*1.05), 
            label = en_coord_cont$Variable
            , colour = "navy", fontface = "bold", size = 2) +
  annotate("text", -0.2, y = 0.1
           , label = paste0("Stress = ", round(nmds1$stress, 2))
           , size = 4) 
  #geom_text(aes(x = -0.2, y = 0.1
  #              , label = paste0("Stress = ", round(nmds1$stress, 2)))
  #          , colour = "black")

Community_comp


# Are lignocellulases correlated with particular breakdown products?
# TODO

PermMANOVA_glucose <- adonis2(species_matrix_rel ~ glucose
                            , data = mf
                            , by = NULL
                            , parallel = 3 # Make it faster
                            , permutations = 60000)

PermMANOVA_hemicellulose_breakdown <- adonis2(species_matrix_rel ~ Hemicellulose_breakdown_products
                             , data = mf
                             , by = NULL
                             , parallel = 3 # Make it faster
                             , permutations = 60000)

PermMANOVA_Lignin_breakdown <- adonis2(species_matrix_rel ~ Lignin_breakdown_products
                                   , data = mf
                                   , by = NULL
                                   , parallel = 3 # Make it faster
                                   , permutations = 60000)

PermMANOVA_Cellulose <- adonis2(species_matrix_rel ~ Cellulose
                              , data = mf
                              , by = NULL
                              , parallel = 3 # Make it faster
                              , permutations = 60000)

PermMANOVA_Hemicellulose <- adonis2(species_matrix_rel ~ Hemicellulose
                                              , data = mf
                                              , by = NULL
                                              , parallel = 3 # Make it faster
                                              , permutations = 60000)

PermMANOVA_Lignin <- adonis2(species_matrix_rel ~ Lignin
                                , data = mf
                                , by = NULL
                                , parallel = 3 # Make it faster
                                , permutations = 60000)


PermMANOVA_glucose 

PermMANOVA_hemicellulose_breakdown 

PermMANOVA_Lignin_breakdown

PermMANOVA_Cellulose 

PermMANOVA_Hemicellulose

PermMANOVA_Lignin 


#### (4) How do community composition scores correlate with metabolite abundance? ####


nmds_obj <- data.table(scores(nmds1, display = "sites"))
nmds_obj[, Sample := dt$Sample]
nmds_obj[, Treatment := dt$Treatment]
nmds_obj$Sample == mf$Sample
nmds_obj[, Hemicellulose := mf$Hemicellulose_breakdown_products]
nmds_obj[, Lignin := mf$Lignin_breakdown_products]
nmds_obj[, Glucose := mf$glucose]

drop1(lm(Hemicellulose ~ NMDS1 + NMDS2# + MDS1:MDS2
         , data = nmds_obj), test = "F")
drop1(lm(Lignin ~ NMDS1 + NMDS2# + MDS1:MDS2
         , data = nmds_obj), test = "F")
drop1(lm(Glucose ~ NMDS1 + NMDS2 #+ MDS1:MDS2
         , data = nmds_obj), test = "F")

mds_x_plot <- melt(nmds_obj, id = c("NMDS1", "NMDS2", "Sample", "Treatment")
                   , variable.name = "Response", value.name = "y")
mds_x_plot <- melt(mds_x_plot, id = c("Sample", "Treatment", "Response", "y"), variable.name = "Predictor"
     , value.name = "x")

mod <- lm(Hemicellulose ~ NMDS1 + NMDS2, data = nmds_obj)
source("Functions/plotGlm.R")
plotGlm(mod)


source("Functions/get_prediction_frame.R")
newdata <- rbind(get_prediction_frame(data = nmds_obj, linear_predictors = "NMDS1"
                     , fixed_predictors = "NMDS2")
                 , get_prediction_frame(data = nmds_obj, linear_predictors = "NMDS2"
                     , fixed_predictors = "NMDS1"))
preds <- predict(mod, newdata = newdata, se.fit = TRUE)
newdata$fit <- preds$fit
newdata$lwr <- preds$fit - (1.96*preds$se.fit)
newdata$upr <- preds$fit + (1.96*preds$se.fit)
newdata$x <- NA
newdata$x[1:100] <- newdata$NMDS1[1:100]
newdata$x[101:200] <- newdata$NMDS2[101:200]
newdata$Predictor <- rep(c("NMDS1", "NMDS2"), each = 100)
setDT(newdata)
newdata[Predictor == "NMDS1",]
newdata[Predictor == "NMDS2", range(fit)]
newdata[Predictor == "NMDS2", range(x)]
newdata[, y := 1]
newdata[, Response := "Hemicellulose"]

mds_x_plot$Treatment

p <- ggplot(mds_x_plot, aes(x = x , y = y)) +
  geom_point(aes(colour = Treatment), size = 3) +
  facet_grid(Response ~ Predictor) +
  cowplot::theme_cowplot() +
  geom_line(data=newdata#[Predictor == "NMDS1",]
            , aes(y=fit, x =x), colour="blue", size = 1.2) +
  geom_ribbon(data = newdata, aes(ymin=lwr, ymax=upr), colour=NA, alpha=0.3) +
  xlab("Lignicellulolytic NMDS value") +
  ylab("Breakdown product abundance (standardized)") +
  theme(panel.border = element_rect(fill = NULL, colour = "black", size = 1.2)
        , panel.spacing = unit(2, "lines")
        , legend.position = "bottom") +
  scale_colour_manual(values = BlackoutPalette, breaks = levels(mds_x_plot$Treatment))

p

#### 4.1 Plot all of the gene community shifts together ####

library(patchwork)

Cellulase_community_comp <- ggplot(nmds_obj_cel
                                   , aes(x = MDS1
                                         , y = MDS2
                                         , colour = Treatment
                                         , fill = Treatment
                                         , group = Treatment)) +
  geom_text(data = nmds_obj_cel_species, mapping = aes(MDS1, MDS2, label = CAZyme
                                                       , colour = NULL, fill = NULL, group = NULL)
            , size = 2) + 
  stat_chull(alpha = 0.7, colour = NA) +
  geom_point(size = 3) +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  theme_classic() +
  theme(text = element_text(size = 16)
        , axis.title = element_text(size = 14)
        , legend.position = "bottom") +
  annotate("text", x = 0.2, y = 0.25
           , label = paste0("Stress = ", round(nmds1$stress, 2))
           , size = 4) +
  ggtitle("Cellulases")

Xylanase_community_comp <- ggplot(nmds_obj_xyl
                                  , aes(x = MDS1
                                        , y = MDS2
                                        , colour = Treatment
                                        , fill = Treatment
                                        , group = Treatment)) +
  geom_text(data = nmds_obj_xyl_species, mapping = aes(MDS1, MDS2, label = CAZyme
                                                       , colour = NULL, fill = NULL, group = NULL)
            , size = 2) + 
  stat_chull(alpha = 0.7, colour = NA) +
  geom_point(size = 3) +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  theme_classic() +
  theme(text = element_text(size = 16)
        , legend.position = "bottom") +
  annotate("text", x = -0.2, y = 0.3
           , label = paste0("Stress = ", round(nmds1$stress, 2))
           , size = 4#6
  ) +
  ggtitle("Xylanases")

AA_community_comp <- ggplot(nmds_obj_aa
                            , aes(x = MDS1
                                  , y = MDS2
                                  , colour = Treatment
                                  , fill = Treatment
                                  , group = Treatment)) +
  geom_text(data = nmds_obj_aa_species, mapping = aes(MDS1, MDS2, label = CAZyme
                                                      , colour = NULL, fill = NULL, group = NULL)
            , size = 2) + 
  stat_chull(alpha = 0.7, colour = NA) +
  geom_point(size = 3) +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  scale_fill_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  theme_classic() +
  theme(text = element_text(size = 16)
        , legend.position = "bottom") +
  annotate("text", x = -0.25, y = 0.3
           , label = paste0("Stress = ", round(nmds1$stress, 2))
           , size = 4) +
  ggtitle("Auxiliary activities")

Lignocellulase_community_comp <- ggplot(nmds_obj
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
        , legend.position = "bottom") +
  geom_segment(data = en_coord_cont
               , aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2)
               , size = 1, alpha = 0.5, colour = "grey30"
                 , arrow = arrow(ends = "last", length = unit(0.1, "inches"))
  ) +
  geom_text(data = en_coord_cont, aes(x = NMDS1*1.05, y = NMDS2*1.05), 
            label = en_coord_cont$Variable
            , colour = "navy", fontface = "bold", size = 3) +
  annotate("text", x = -0.15, y = 0.1
           , label = paste0("Stress = ", round(nmds1$stress, 2))
           , size = 4) + ggtitle("Lignocellulases")



Cellulase_community_comp
Xylanase_community_comp
AA_community_comp
Lignocellulase_community_comp 

p1 <- (Cellulase_community_comp + theme(legend.position = "none")| 
         Xylanase_community_comp + theme(legend.position = "none")) /
  (AA_community_comp + theme(legend.position = "bottom") + 
     guides(colour = guide_legend(title.position = "top", nrow = 2)
            , fill = guide_legend(title.position = "top", nrow = 2)
            )| 
     Lignocellulase_community_comp + theme(legend.position = "none")) &
  theme(axis.title = element_text(size = 14)
        , title = element_text(size = 14)
        , legend.text = element_text(size = 8)
        , legend.title = element_text(size = 12)) & 
  plot_layout(guides = 'collect') +  
  theme(legend.position = 'bottom')
  
p1

tiff(filename = 'Figures_2023_SBB/Figure_6_Lignocellulase_community_composition.tif'
     , width = 15.5, height = 21
     , units = "cm"
     , res = 600)
p1
dev.off()

#### (5) Which CAZymes were more or less abundant? ####

# Overall abundance
dt[, LC_abundance := rowSums(species_matrix_rel)] 
LC_abundance_results <- univariate_kw_dunn_plot("LC_abundance"
                                                , y_lab = "LC CAZy family abundance"
                                                , Letter_height = 0.08
                                                , p_type = "unadjusted"
                                                #, p_type = "adjusted"
)
LC_abundance_results$Kruskal_wallis
LC_abundance_results$Dunn
LC_abundance_results$Plot 
LC_abundance_results$Effect

# Set the data up for logit-transformed linear modlling
species_matrix_logit <- species_matrix_rel
setDT(species_matrix_logit)
species_matrix_logit[, colnames(species_matrix_logit) := lapply(.SD,FUN = function(x){x <- x/1e6})]
species_matrix_logit[, colnames(species_matrix_logit) := lapply(.SD, car::logit)]

# Make a bunch of lms on logit-transformed data
modlist <- lapply(seq_along(species_matrix_logit), FUN = function(i){
  mod <- lm(species_matrix_logit[[i]] ~ Treatment)
  })
names(modlist) <-  colnames(species_matrix_logit)


# Calculate fold change
library(magrittr)
fcs <- lapply(species_matrix_rel, FUN = function(x){
  grassland_mean <- mean(x[Treatment == '10-Year\nGrassland'])
  FC <- x / grassland_mean
})
names(fcs) <- colnames(species_matrix_rel)
setDT(fcs)

fcs[, Treatment := Treatment]
fc_means <- fcs[, .(FC = sapply(.SD, mean)
                    , FC_ci = sapply(.SD, ci)
                    )
                , by = Treatment]
fc_means[, CAZyme := rep(colnames(species_matrix_rel), 4)]
fc_means[, upr := FC + FC_ci]
fc_means[, lwr := FC - FC_ci]
fc_means[, log2_upr := log2(upr)]
fc_means[, log2_lwr := log2(lwr)]
fc_means[, log2_FC := log2(FC)]


# Get the p-values of each CAZyme's fold change
p_tab <- data.table(CAZyme = names(modlist)
           , p = lapply(modlist, function(x){
             res <- drop1(x, test = "F")
             res[2,6]
           }))
l2fcs <- merge(p_tab, fc_means, by = "CAZyme")
l2fcs[, isInf := ifelse(is.infinite(log2_FC), yes = 'Not present in\n10-Year Grassland', no = "Present in\n10-Year Grassland")]

# Which CAZy families met the criteria for being interesting?
sigCAZ <- l2fcs[(log2_FC > 1 | log2_FC < -1) 
                & p < 0.05, CAZyme]


# Plot the results
p_fc <- ggplot(l2fcs[Treatment != "10-Year\nGrassland" & CAZyme %in% sigCAZ
             & log2_FC != -Inf,]
       , aes(x = log2_FC, y = CAZyme, colour = Treatment
             , shape = isInf)) +
  geom_point(aes(fill = Treatment), size  = 3) +
  geom_point(data = l2fcs[Treatment != "10-Year\nGrassland" & CAZyme %in% sigCAZ
                         & log2_FC == Inf,]
             , size  = 3, colour = "black") +
  
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_vline(xintercept = 0) +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  scale_x_continuous(name = expression("Log"[2]*"(Fold change) Relative to 10-year grassland")
                     , expand = c(0,0)
                     , limits = c(-5, 5)) +
  ylab("CAZy family") +
  scale_shape_manual(name = "", values = c(22, 16), guide = "none") +
  #guides(shape = NULL) +
  cowplot::theme_cowplot() +
  theme(text = element_text(size = 16)
        , panel.border = element_rect(colour = "black", size = 1.1)
        , legend.position = "bottom"
        , plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "lines")
  )

p_fc

### Find more stats about what proportion of gene families responded


active_families_list$Cellulases
# Significant fold change in 5 / 19 GH5 containing CAZy combinations (26%)
active_families_list$Cellulases[grep("GH5", active_families_list$Cellulases)]
length(grep("GH5", active_families_list$Cellulases))
(5/19)*100

length(grep("GH45", active_families_list$Cellulases))

# 1/2 GH11
length(grep("GH11", active_families_list$Cellulases))
active_families_list$Xylanases[grep("GH11", active_families_list$Xylanases)]

length(grep("GH9", active_families_list$Xylanases))
active_families_list$Xylanases[grep("GH11", active_families_list$Xylanases)]


# 1 / 2
active_families_list$LPMOs


active_families_list 
active_families_list 







