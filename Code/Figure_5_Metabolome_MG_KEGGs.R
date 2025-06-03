#### Introduction ####

# Here we make the 2-panel figure of metabolome NMDS and the fermentation means 
# result

#### (0) Setup ####

# Set things up
rm(list = ls())

library(data.table)
library(magrittr)
library(vegan)
library(ggplot2)


# Allow everyone to run everyone else's scripts without rewriting file paths :)
if(file.exists("C:/Users/lcl24grb/OneDrive - Bangor University/Blackout_David_Nov_24/")){
  bl_path <- "C:/Users/lcl24grb/OneDrive - Bangor University/Blackout_David_Nov_24/"
}
if(file.exists("C:/Users/rbg23rfv/OneDrive - Bangor University/Blackout_David_Nov_24/")){
  bl_path <- "C:/Users/rbg23rfv/OneDrive - Bangor University/Blackout_David_Nov_24/"
}
if(file.exists("C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/RG_script_KEGG.R")){
  bl_path <- "C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/"
}

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

sample_cols <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7'
                 , 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7')

Treatment <- factor(c(rep("b_old", 3), rep("b_new", 4)
                      , rep("c_old", 3), rep("c_new", 4))
                    , levels = c("b_old", "b_new", "c_new", "c_old"))

# Watch this one. ordered as for levels
Treatment_labs <- c("10-Year\nBare", "1-Year\nBare", "1-Year\nGrassland", "10-Year\nGrassland")


#### (1) Fermentation analysis ####


##### Data import #####
# Foam DB
FOAM <- fread("C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/Fermentation/FOAM.tsv")
# This is where information about the function of different KOs comes from 

# Metacerberus annotation
mc <- fread("C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/Fermentation/henfaes_contig_HMMER_top_5.tsv")
# This is where we store contig associations with KEGG IDs

# Get the fermentation KEGGs, and ensure we only unique KEGG IDs to prevent double counting
KEGGs <- FOAM[grep("Fermentation", L1), ID]
KEGGs <- unique(KEGGs)

# Get the contigs associated with these KEGGs 
contigs_to_keep <- mc[ID %in% KEGGs, `Target Name`]

# Read in the data with contig labels and relative abundances
k <- fread("C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/TPM_contigs/Blackout_contigs_TPM.csv")
# This is where we analyse the relative abundance mapping to each sample

head(k)

##### Data cleaning #####

# Subset to contigs containing a fermentation FOAM KEGG.
k <- k[Contig %in% sub("_[0-9]+$", "", contigs_to_keep),];gc()
k <- k[n1 == 0,] # Make sure nothing mapped to in the negative is included


# Get total relative abundance of fermentation genes
ferm_vals <- k[, lapply(.SD, sum), .SDcols = sample_cols]

# Make this plottable
fv <- data.table(Sample = sample_cols
                 , Treatment
                 , CPM = unlist(ferm_vals))

#### L1 analysis ####

# View the outcome
ggplot(fv, aes(x = Treatment, y = CPM)) +
  geom_point() +
  stat_summary()

# What do conservative stats say?
kruskal.test(fv$CPM, g = Treatment)

# So there is no signal in the L1 fermentation KEGG category

#### L2 analysis ####

# Let's go a KEGG level down and see if more specific functions are affected

###### Function for subsetting the data ######

# This handy function can subset and sum the read conts associated with any K numbers

sum_KEGGs <- function(KEGG_ID){
  # Get the contigs associated with input KEGGs 
  contigs_to_keep <- mc[ID %in% KEGG_ID, `Target Name`]
  
  # Subset to contigs containing a FOAM KEGG.
  k_sub <- k[Contig %in% sub("_[0-9]+$", "", contigs_to_keep),]
  
  # Sum the CPM values from these contigs
  kegg_vals <- k_sub[, lapply(.SD, sum), .SDcols = sample_cols]
  
  # Format better
  kv <- data.table(Sample = sample_cols
                   , Treatment
                   , CPM = unlist(kegg_vals))
  return(kv)
}

###### L2-wise KEGG data cleaning ######

# For each KO, get the relative abundance
names(KEGGs) <- KEGGs

MG <- lapply(KEGGs, sum_KEGGs) %>% rbindlist(., idcol = "KO")
MG

# Now make this into a wide table
MGw <- dcast(MG, Sample + Treatment ~ KO
             , value.var = "CPM")
MGw

#MGw <-

# Remove fully zero columns
absent_KEGGs <- MGw[ , lapply(.SD, sum), .SDcols = KEGGs]
absent_KEGGs <- names(absent_KEGGs)[absent_KEGGs == 0]
MGw[, (absent_KEGGs) := NULL]
# Keep non-zero columns
KEGGs_present <- KEGGs[!KEGGs %in% absent_KEGGs]

# Now put back into long format for further analysis
# N.B., this is not the same as MW as contigs have multiple annotations
# The sum of all values from MGl should therefore be greater than the sum of
# CPM in our original dataset, fv.
MGl <- melt(MGw, id = c("Sample", "Treatment"))
setnames(MGl, old = 'variable', new = 'ID')

if(all(MGl[, sum(value), by = Sample][, V1] > fv$CPM)){
  paste("All CPM values in MGl are greater than in fv: continue analysis")
  }

##### L2 community analysis ####

# Do an NMDS
nmds_mg <- metaMDS(MGw[, ..KEGGs_present])

# Check the NMDS
plot(nmds_mg, "sites")
ordispider(nmds_mg, Treatment, label = TRUE)

# There are obvious huge differences in composition

##### L2 gene family analysis #####

# What are the L2 categories?
unique(FOAM[grep("Fermentation", L1), L2])

# Get the KEGG subsets to analyse
L2_functions <- unique(FOAM[grep("Fermentation", L1), L2])
names(L2_functions) <- L2_functions
KEGG_subsets <- lapply(L2_functions, FUN = function(x){
  out <- FOAM[grepl("Fermentation", L1) & L2 == x, ID]
  return(out)
})

# Make a dataset for each of these functions
dfL2 <- lapply(KEGG_subsets, FUN = function(x){
  MGl[ID %in% x, .(CPM = sum(value)), by = .(Sample, Treatment)]
}) %>% rbindlist(., idcol = "L2")

# Do stats on each fermentation process, then make the factor 
# levels represent the level of statistical support that there are 
# treatment differences
kw_L2 <- dfL2[, kruskal.test(CPM, Treatment), by = L2
              ][order(p.value, decreasing = TRUE),]
sig_order <- kw_L2[, L2]
dfL2[, L2 := factor(L2, levels = rev(sig_order))]


# Visualise the relative abundances of each fermentation process
ggplot(dfL2, aes(x = Treatment, y = CPM, colour = Treatment)) +
  stat_summary(fun.data = 'mean_se', position = position_dodge(0.9)) +
  theme_classic() +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  facet_wrap(L2 ~ ., scales = 'free_y') +
  ggtitle("Fermentation L2 Categories")

kw_L2[p.value < 0.05,]

glm_L2 <- dfL2[, glm(CPM ~ Treatment) %>% drop1(., test = 'Chisq'), by = L2
]
glm_L2 <- glm_L2[seq(2, nrow(glm_L2), by = 2)
                 ][, p_adj := p.adjust(`Pr(>Chi)`, 'fdr')
                   ][order(p_adj, decreasing = TRUE),]
glm_L2[p_adj < 0.05,]


sig_order <- glm_L2[, L2]
dfL2[, L2 := factor(L2, levels = rev(sig_order))]

# Visualise the relative abundances of each fermentation process
ggplot(dfL2, aes(x = Treatment, y = CPM, colour = Treatment)) +
  stat_summary(fun.data = 'mean_se', position = position_dodge(0.9)) +
  theme_classic() +
  scale_colour_manual(breaks = levels(Treatment), values = BlackoutPalette) +
  facet_wrap(L2 ~ ., scales = 'free_y') +
  ggtitle("Fermentation L2 Categories")

# Add significance indicators to the plot
dfL2[, Significance := ifelse(L2 %in% glm_L2[p_adj < 0.05,]$L2
              , yes = "p < 0.05"
              , no = "N.S.")][, Significance := factor(Significance
                                                       , levels = c("p < 0.05"
                                                                    , "N.S."))]
# Sort the categories by abundance
abundance_order <- dfL2[, mean(CPM), by = L2][order(V1),][,L2]# L2]
dfL2[, L2 := factor(L2, levels = rev(abundance_order))]

dfL2[, dput(unique(L2))]

subber <- data.table(L2 = c("Pyruvate fermentation to butanol", "Acetoin biosynthesis", 
  "Mixed acid fermentation", "Heterolactic fermentation", "Glutamate degradation V (via hydroxyglutarate)", 
  "Pyruvate Fermentation", "Lysine fermentation to acetate and butyrate", 
  "Homolactic fermentation", "Bifidobacterium shunt", "Acetylene degradation", 
  "Succinate fermentation to butyrate", "Acetyl-coA fermentation to butyrate II (clostridium)"
)
, L2_short = c("Pyruvate to BuOH", "Acetoin synthesis", 
               "Mixed acid", "Heterolactic", "Glutamate V", 
               "Pyruvate", "Lysine", 
               "Homolactic", "B. shunt", "Acetylene", 
               "Succinate to C4:0", "Acetyl-coA to C4:0"
)
)

dfL2s <- merge.data.table(dfL2, subber, by = 'L2')
dfL2s[, L2_short := factor(L2_short, levels = subber$L2_short)]


# Make the final fermentation plot
p_ferm <- ggplot(dfL2s, aes(x = Treatment, y = CPM
                 , colour = L2_short
                 , group = L2_short
                 , linetype = Significance)) +
  stat_summary(fun.data = 'mean_cl_boot'
               , geom = "errorbar"
               , width = 0.1
               , linetype = 'solid'
              ) +
  stat_summary(fun = 'mean'
               , geom = 'path'
               #, position = position_dodge(0.9)
  ) +
  theme_classic() +
  scale_colour_discrete("Fermentation L2 Categories") + 
  scale_y_log10() +
  scale_x_discrete(labels = Treatment_labs) +
  theme(legend.position = 'bottom'
        , strip.text.y.right = element_text(angle = 0)
        , text = element_text(size = 16)
        , plot.margin = margin(0,0,0,0)
        #, panel.margin.x = unit(-2, "lines")
        ) +
  guides(colour = 'none') + 
  facet_wrap(L2_short ~ .
             , scales = 'free_y'
             , ncol = 1
             , strip.position = "right"
             )

p_ferm

#### (2) Metabolome ####

##### Setup #####

library(vegan)

library(devtools)
source_gist("https://gist.github.com/robiwangriff/79633a738b128e64226bf58381232da7")
library(labdsv)

source("Functions/stat_chull.R")

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

sample_cols <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7'
                 , 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7')
Treatment <- factor(c(rep("bare_old", 3), rep("bare_new", 4)
                      , rep("control_old", 3), rep("control_new", 4))
                    , levels = c("bare_old", "bare_new", "control_new", "control_old"))

Treatment_labs <- c("10-Year\nBare", "1-Year\nBare", "1-Year\nGrassland", "10-Year\nGrassland")

trt<-c("bare_old","bare_old","bare_old","bare_new",
       "bare_new","bare_new","bare_new",
       "control_old","control_old","control_old",
       "control_new","control_new","control_new","control_new")

##### Data import and cleaning #####

base_path <- "C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/"
metabolome_path <- paste0(base_path, "Metabolomics/")

metab<-read.csv(paste0(metabolome_path, "Metabolome_Blackout_Raw.csv"), row.names=1,check.names=F)
metab<-t(metab)
metab.tax<-data.frame(id=paste0("X",1:ncol(metab)),metab=colnames(metab),stringsAsFactors=FALSE)
rownames(metab.tax)<-paste0("X",1:ncol(metab))
colnames(metab)<-paste0("X",1:ncol(metab))

#read wider metabolite info (KEGG ids etc)
metab_tax2<-read.csv(paste0(metabolome_path, "Metabolome_Blackout_tax.csv"),stringsAsFactors=FALSE)
#get higher level class
library(omu)
metab_tax2<- assign_hierarchy(count_data =metab_tax2, keep_unknowns = TRUE, identifier = "KEGG")
metab_tax2[] <- lapply(metab_tax2, function(x) if(is.factor(x)) as.character(x) else x)
metab.tax<-merge(metab.tax,metab_tax2,by.x="metab", by.y="BinBase.name")
library(stringr)
# Sort by numeric part after 'X'
metab.tax<- metab.tax[order(as.numeric(str_extract(metab.tax$id, "\\d+"))), ]
rownames(metab.tax)<-metab.tax$id
###
##analyses

#nmds
#mod6<-metaMDS(decostand(metab,"log"))

# Make sure we are comparing like with like--maybe not necessary but just in case.
# Doesn't really change outcome
mod6<-metaMDS(decostand(decostand(metab, "total"),"log"))
plot(mod6,"species")
text(mod6,"species")
plot(mod6,"sites")
ordispider(mod6,trt,label=T)

#Group vars
#change NA to unclassified

metab.tax$Subclass_1 <- as.character(metab.tax$Subclass_1)        # Convert factor to character
metab.tax$Subclass_1[is.na(metab.tax$Subclass_1)] <- "unclassified"
metab.tax$Subclass_1<-factor(metab.tax$Subclass_1)

# David edit 2024-12-12. Removing "none" as a compound class
metab.tax_dt <- as.data.table(metab.tax)
metab.tax_dt[Subclass_1 == 'none', Subclass_1 := Class]
#metab.tax_dt[Subclass_1 == 'none',]

metab.tax_dt_df <- as.data.frame(metab.tax_dt, row.names = row.names(metab.tax))
agg1 <- aggreg8(metab, metab.tax_dt_df ,"Subclass_1")

agg1

##### NMDS #####

ef<-envfit(mod6,agg1)
plot(mod6,"sites")
ordispider(mod6,trt,label=T)
plot(ef,p.max=0.05)


# MAke these able to be used with ggplot
en_coord_cont_Metabolome = as.data.frame(vegan::scores(ef, "vectors")) * ordiArrowMul(ef) #* 0.8
en_coord_cont_Metabolome$Treatment <- "bare_old"#"10-Year\nBare"
setDT(en_coord_cont_Metabolome)
en_coord_cont_Metabolome[, Variable := row.names(vegan::scores(ef, "vectors"))
                         #c("Carbon", "N:C", "Total cations", "Hemicellulose breakdown\nproducts"
                         # , "Lignin breakdown\nproducts", "Glucose")
]
en_coord_cont_Metabolome[, p_value := ef$vectors$pvals]
en_coord_cont_Metabolome[, p_value_adj := p.adjust(p_value, method = "fdr")]
en_coord_cont_Metabolome[, Significance := "N.S.",]
# en_coord_cont_Metabolome[p_value_adj < 0.1, Significance := ifelse(p_value_adj < 0.05
#                                                                    , yes = "p < 0.05"
#                                                                    , no = "p < 0.1"),]
# en_coord_cont_Metabolome_s <- en_coord_cont_Metabolome[p_value_adj < 0.1,]
# en_coord_cont_Metabolome_s <- en_coord_cont_Metabolome_s[Variable != "unclassified",]

en_coord_cont_Metabolome[p_value < 0.1, Significance := ifelse(p_value < 0.05
                                                                   , yes = "p < 0.05"
                                                                   , no = "p < 0.1"),]
en_coord_cont_Metabolome_s <- en_coord_cont_Metabolome[p_value < 0.1,]
en_coord_cont_Metabolome_s <- en_coord_cont_Metabolome_s[Variable != "unclassified",]


library(ggplot2)
library(ggrepel)
library(data.table)
# Create the ordination plot with ggrepel to avoid label overlap # not finished this yet!
df <- data.table(Sample = sample_cols
                 , Treatment = trt  
                 , vegan::scores(mod6, "sites", choices = 1:2) 
)
setnames(df, old = c("NMDS1", "NMDS2"), new = c("MDS1", "MDS2"))

na.omit(metab.tax[metab.tax$Subclass_2 == 'none',])

##### Metabolite plotting #####

M_comp <- ggplot(df
                 , aes(x = MDS1
                       , y = MDS2
                       , colour = Treatment
                       , fill = Treatment
                       , group = Treatment)) +
  stat_chull(alpha = 0.7, colour = NA) +
  geom_point(size = 3) +
  scale_colour_manual(breaks = levels(Treatment)
                      , values = BlackoutPalette
                      , labels = Treatment_labs) +
  scale_fill_manual(breaks = levels(Treatment)
                    , values = BlackoutPalette
                    , labels = Treatment_labs) +
  theme_classic() +
  geom_segment(data = en_coord_cont_Metabolome_s
               , aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2
                     , linetype = Significance)
               , size = 1, alpha = 0.5, colour = "grey30"
               , arrow = arrow(ends = "last", length = unit(0.1, "inches"))
  ) +
  geom_text_repel(data = en_coord_cont_Metabolome_s, aes(x = NMDS1*1.05, y = NMDS2*1.05), 
                  label = en_coord_cont_Metabolome_s$Variable
                  , colour = "navy", fontface = "bold", size = 3) +
  annotate("text", -0.015, y = 0.022
           , label = paste0("Stress = ", round(mod6$stress, 2))
           , size = 4) +
  theme(text = element_text(size = 10#16
                            )
        , legend.position = "bottom"
        , legend.box = "vertical"
        , plot.margin = margin(0,0,0,0)
        ) +
  guides(colour = guide_legend(nrow = 2 , byrow = TRUE))


M_comp

#en_coord_cont_Metabolome_s[Variable != "Oligosaccharides", Variable]
Grassland_compounds <- en_coord_cont_Metabolome_s[NMDS1 > 0,]
Bare_compounds <- en_coord_cont_Metabolome_s[NMDS1 < 0,]

Grassland_compounds
Bare_compounds[abs(NMDS1) > 0.01,][order(NMDS1),]

en_coord_cont_Metabolome_s[abs(NMDS1) > 0.01,][order(NMDS1),]

#### (3) Metabolites (aggregated)  ####

MGw <- fread("C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/TPM_contigs/Blackout_KOs_by_sample_TPM.csv")

spc <- as.data.frame(MGw[,3:ncol(MGw)])

#####################
# now an env dataframe

env<-data.frame(treat=Treatment)

#read in foam ontology

# Foam DB
FOAM <- fread(paste0(bl_path,"Fermentation/FOAM.tsv"))
FOAM<-as.data.frame(FOAM)

#use aggreg8_mg a custom function to aggregate, handling multiple genes being in different pathways
source_gist("https://gist.github.com/robiwangriff/79633a738b128e64226bf58381232da7")

# aggs<-aggreg8_mg(spc=spc,ann=FOAM,level_col="L1",gene_col="ID")
# 
# #####plot
# #add treat column to agg df
# aggs$treat<-env$treat
# 
# #make long
# aggs.l<-melt(aggs,id.vars="treat",variable.name="L1",value.name="count")
# #
# ggplot(aggs.l,aes(x=treat,y=count,fill=treat))+
#   geom_boxplot()+
#   facet_wrap(~L1,scales="free")+
#   theme_minimal()

#show only L1 =fermentation
aggs<-aggreg8_mg(spc=spc,ann=FOAM,level_col="L2",gene_col="ID")

#####plot
#add treat column to agg df
aggs$treat<-env$treat

#make long
aggs.l<-melt(aggs,id.vars="treat",variable.name="L2",value.name="count")

sel<-unique(FOAM[FOAM$L1=="01_Fermentation",]$L2)

ggplot(aggs.l[aggs.l$L2%in%sel,],aes(x=treat,y=count,fill=treat))+
  geom_boxplot()+
  facet_wrap(~L2,scales="free")+
  theme_minimal()

#add to metabolome ordination

metab<-read.csv(paste0(bl_path, "Metabolomics/Metabolome_Blackout_Raw.csv"), row.names=1,check.names=F)
metab<-t(metab)

trt<-c("bare_old","bare_old","bare_old","bare_new",
       "bare_new","bare_new","bare_new",
       "control_old","control_old","control_old",
       "control_new","control_new","control_new","control_new")

#nmds
#mod6<-metaMDS(decostand(metab,"log"))
mod6<-metaMDS(metab)
#plot(mod6,"species")
#text(mod6,"species")
plot(mod6,"sites")
ordispider(mod6,trt,label=T)

#envfit
ef1<-envfit(mod6,aggreg8_mg(spc=spc,ann=FOAM,level_col="L1",gene_col="ID"))
ef2<-envfit(mod6,aggreg8_mg(spc=spc,ann=FOAM,level_col="L2",gene_col="ID"))
plot(ef2,p.max=0.05,cex=0.7)

setDT(FOAM)
FOAM[L2 %in% names(ef2$vectors$pvals[ef2$vectors$pvals < 0.05]), table(unique(L1))]

#### (4) Metabolites and metagenome correlations ####

#rm(list = ls())

##### Setup #####

library(data.table)
library(magrittr)
library(vegan)
library(ggplot2)
library(devtools)

BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

sample_cols <- c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7'
                 , 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7')

Treatment <- factor(c(rep("b_old", 3), rep("b_new", 4)
                      , rep("c_old", 3), rep("c_new", 4))
                    , levels = c("b_old", "b_new", "c_new", "c_old"))


# Allow everyone to run everyone else's scripts without rewriting file paths :)
if(file.exists("C:/Users/lcl24grb/OneDrive - Bangor University/Blackout_David_Nov_24/")){
  bl_path <- "C:/Users/lcl24grb/OneDrive - Bangor University/Blackout_David_Nov_24/"
}
if(file.exists("C:/Users/rbg23rfv/OneDrive - Bangor University/Blackout_David_Nov_24/")){
  bl_path <- "C:/Users/rbg23rfv/OneDrive - Bangor University/Blackout_David_Nov_24/"
}
if(file.exists("C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/RG_script_KEGG.R")){
  bl_path <- "C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/"
}


MGw <- fread(paste0(bl_path, "TPM_contigs/Blackout_KOs_by_sample_TPM.csv"))
spc<-as.data.frame(MGw[,3:ncol(MGw)])


# Foam DB
FOAM <- fread(paste0(bl_path,"Fermentation/FOAM.tsv"))
FOAM<-as.data.frame(FOAM)

#use aggreg8_mg a custom function to aggregate, handling multiple genes being in different pathways
source_gist("https://gist.github.com/robiwangriff/79633a738b128e64226bf58381232da7")

aggs <- aggreg8_mg(spc = spc # KOs by sample 
                   , ann = FOAM # FOAM annotations
                   , level_col = "L2" # Level to summarise at
                   , gene_col = "ID" # Gene annotation column from "ann"
)

#####plot
#add treat column to agg df
env <- data.frame(treat=Treatment)
aggs$treat <- env$treat

#make long
aggs.l <- melt(aggs,id.vars="treat",variable.name="L2",value.name="count")

sel <- unique(FOAM[FOAM$L1=="01_Fermentation",]$L2)

ggplot(aggs.l[aggs.l$L2%in%sel,], aes(x = treat, y = count, fill = treat))+
  geom_boxplot()+
  facet_wrap(~ L2, scales = "free")+
  theme_minimal()

metab<-read.csv(paste0(bl_path, "Metabolomics/Metabolome_Blackout_Raw.csv"), row.names = 1, check.names = F)
metab<-t(metab)
metab <- decostand(decostand(metab, "total"),"log")
trt<-c("bare_old","bare_old","bare_old","bare_new",
       "bare_new","bare_new","bare_new",
       "control_old","control_old","control_old",
       "control_new","control_new","control_new","control_new")

##### Community ecology #####

#nmds
#mod6<-metaMDS(decostand(metab,"log"))
mod6<-metaMDS(metab)
#plot(mod6,"species")
#text(mod6,"species")
plot(mod6,"sites")
ordispider(mod6,trt,label=T)

#envfit
ef1<-envfit(mod6,aggreg8_mg(spc=spc,ann=FOAM,level_col="L1",gene_col="ID"))
#ef2<-envfit(mod6,aggreg8_mg(spc=spc,ann=FOAM,level_col="L2",gene_col="ID"))
efd <- aggreg8_mg(spc=spc,ann=FOAM,level_col="L2",gene_col="ID")
efd <- efd[, colSums(efd) != 0]
ef2 <- envfit(mod6, efd)

#hist(colSums(efd))
#sum(colSums(efd) == 0)
plot(ef2,p.max=0.05,cex=0.7)



##### Plotting #####

df_mm <- vegan::scores(mod6, display = 'sites', choices = 1:2)
df_mm <- data.table(df_mm)
setnames(df_mm, c("MDS1", "MDS2"))
df_mm[, Sample := sample_cols]
df_mm[, Treatment := Treatment]
#df_mm[, Treatment_labs := Treatment_labs]

# MAke these able to be used with ggplot
en_coord_cont_mm = as.data.frame(vegan::scores(ef2, "vectors")) * ordiArrowMul(ef2) #* 0.8
en_coord_cont_mm$Treatment <- "b_old"#"bare_old"#"10-Year\nBare"
setDT(en_coord_cont_mm)
en_coord_cont_mm[, Variable := row.names(vegan::scores(ef2, "vectors"))
                         #c("Carbon", "N:C", "Total cations", "Hemicellulose breakdown\nproducts"
                         # , "Lignin breakdown\nproducts", "Glucose")
]
en_coord_cont_mm[, p_value := ef2$vectors$pvals]
en_coord_cont_mm[, p_value_adj := p.adjust(p_value, method = "fdr")]
en_coord_cont_mm[, Significance := "N.S.",]
en_coord_cont_mm[p_value < 0.1, Significance := ifelse(p_value < 0.05
                                                                   , yes = "p < 0.05"
                                                                   , no = "p < 0.1"),]
en_coord_cont_mm_s <- en_coord_cont_mm[p_value < 0.1,]
en_coord_cont_mm_s <- en_coord_cont_mm_s[Variable != "unclassified",]

unclear <- en_coord_cont_mm_s[grep("^From", Variable), Variable]
class(FOAM)
FOAMdt <- as.data.table(FOAM)
drop_in <- unique(FOAMdt[L2 %in% unclear, .(L2, paste(L1, tolower(L2)))])

en_coord_cont_mm_s[, Variable]

unique(FOAMdt[L2 %in% c('Cellulose', 'Xylose', 'Hexose'), .(L1, L2)])

# Ensure everything is intelligable
lapply(1:nrow(drop_in), FUN = function(i){
 en_coord_cont_mm_s[Variable == drop_in[i, L2], Variable := drop_in[i, V2] ]
})

en_coord_cont_mm_s[, Variable := sub("^[0-9]+_", "", Variable)]
en_coord_cont_mm_s[Variable == 'Cellulose', Variable := "Cellulose hydrolysis"]
en_coord_cont_mm_s[Variable == 'Xylose', Variable := "Polymer hydrolysis (xylose)"]
en_coord_cont_mm_s[Variable == 'Hexose', Variable := "Glycolysis (conversion of pentose to hexose)"]

en_coord_cont_mm_s[NMDS2 > 0, Variable]

#min(en_coord_cont_mm$p_value_adj)
#hist(en_coord_cont_mm$p_value_adj)
#hist(en_coord_cont_mm$p_value)

# ggplot2ify

# Metabolome scores of full set of metabolites measured

mm_comp <- ggplot(df_mm
                 , aes(x = MDS1
                       , y = MDS2
                       , colour = Treatment
                       , fill = Treatment
                       , group = Treatment)) +
  stat_chull(alpha = 0.7, colour = NA) +
  geom_point(size = 3) +
  scale_colour_manual(breaks = levels(Treatment)
                      , values = BlackoutPalette
                      , labels = Treatment_labs) +
  scale_fill_manual(breaks = levels(Treatment)
                    , values = BlackoutPalette
                    , labels = Treatment_labs) +
  theme_classic() +
  geom_segment(data = en_coord_cont_mm_s
               , aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2
                     , linetype = Significance, colour = NULL)
               , size = 1, alpha = 0.5, colour = "grey30"
               , arrow = arrow(ends = "last", length = unit(0.1, "inches"))
  ) +
  geom_text_repel(data = en_coord_cont_mm_s
                  , aes(x = NMDS1*1.05, y = NMDS2*1.05
                        , colour = NULL), 
                  label = en_coord_cont_mm_s$Variable
                  , colour = "navy", fontface = "bold"
                  , size = 1.8
                  , max.overlaps = 1000) +
  #annotate("text", -0.015, y = 0.022
  #         , label = paste0("Stress = ", round(mod6$stress, 2))
  #         , size = 4) +
  theme(text = element_text(size = 10#16
  )
  , legend.position = "bottom"
  , legend.box = "vertical"
  , plot.margin = margin(0,0,0,0)
  ) +
  guides(colour = guide_legend(nrow = 2 , byrow = TRUE))


mm_comp


Grassland_L2s <- en_coord_cont_mm[NMDS2 > 0 & NMDS1 > 0 & p_value < 0.05
                                  , Variable]
# What are grassland associated processes?
unique(FOAMdt[L2 %in% Grassland_L2s, .(L1, L2)])[, .N, by = .(L1)][ order(N),]
Grassland_L2s

#FOAMdt[grep("none", L2), .(L1, L2)]

# Methanogenesis, methylotrophy, N, amino acids, polymers, 

# What are bare associated processes?
Bare_L2s <- en_coord_cont_mm[NMDS2 < 0 & p_value < 0.05
                                  , Variable]
unique(FOAMdt[L2 %in% Bare_L2s, .(L1, L2)])
unique(FOAMdt[L2 %in% Bare_L2s, L2])

unique(FOAMdt[L2 %in% Bare_L2s, .(L1, L2)])[, .N, by = .(L1)][ order(N),]


#### (5) Plotting together ####

library(patchwork)

# tiff("C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/Figures/Metabolome_and fermentation.tiff"
#      , width = 7
#      , height = 7
#      , units = "in"
#      , res = 600
#      )
# M_comp | p_ferm * theme(text = element_text(size = 10)
#                         , plot.margin = unit(c(0, 0, 0, 0), "cm")
#                         ) +
#   patchwork::plot_annotation(tag_levels = "A") +
#   patchwork::plot_layout(widths = c(0.6, 0.5)
#                          #, margin = unit(0, "cm")
#                          )
# dev.off()

ggsave(filename = "Figures_2024/Figure_5_Metabolome_L1_KEGGs.tiff"
       , plot = M_comp
       , device = "tiff"
       , width = 6
       , height = 6
       , dpi = 600
)


combined_plot <- (M_comp | mm_comp & theme(text = element_text(size = 10)
                                 , plot.margin = unit(c(0, 0.5, 0, 0), "cm")
                                 , legend.text = element_text(size = 8)
                                 , legend.position = 'bottom'
                                  ) 
                  ##  | p_ferm + theme(text = element_text(size = 10)
                        # , plot.margin = unit(c(0, 0, 0, 0), "cm")
                        # , legend.text = element_text(size = 8)
                        #)
                  ) / guide_area() +
  patchwork::plot_layout(guides = 'collect'
                         , heights = c(1, 0.3)) +
  patchwork::plot_annotation(tag_levels = "A")

combined_plot

ggsave(filename = "C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/Figures/Metabolome_and_metagenome_KEGGs.tiff"
  , plot = combined_plot
  , device = "tiff"
  , width = 7
  , height = 6
  , dpi = 600
)

# ggsave(filename = "C:/Users/bspa44/OneDrive - Bangor University/Blackout_David_Nov_24/Figures/Metabolome_and_metagenome_KEGGs_test.tiff"
#        , plot = combined_plot
#        , device = "tiff"
#        , width = 7
#        , height = 6
#        , dpi = 600
# )
