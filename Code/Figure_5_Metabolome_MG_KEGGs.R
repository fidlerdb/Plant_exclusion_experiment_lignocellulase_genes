##### Setup #####
rm(list = ls())

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