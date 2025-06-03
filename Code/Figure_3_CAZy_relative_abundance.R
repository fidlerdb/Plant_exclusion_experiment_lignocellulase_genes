spc<-read.csv("Taxonomy/Multi-Method/Cleaned_Data/Henfaes_Species_Multi-Method_TPM_1kbpmin.csv")
tax<-spc[,1:8]
spc<-spc[,9:ncol(spc)]
spc<-t(spc)
trt<-c("bare_old","bare_old","bare_old","bare_new",
       "bare_new","bare_new","bare_new",
       "control_old","control_old","control_old",
       "control_new","control_new","control_new","control_new")

trt<- factor(c(rep("10-Year\nBare", 3), rep("1-Year\nBare", 4),
       rep("10-Year\nGrassland", 3), rep("1-Year\nGrassland", 4))
       , levels = c("10-Year\nBare", "1-Year\nBare"
                    , "1-Year\nGrassland", "10-Year\nGrassland")
       )

meta<-read.csv("Metabolomics/Cleaned Data/Henfaes_Taxonomy_Metabolome.csv",row.names=1)

aa <- read.csv("CAZy_2/outputData/Blackout_aa_matrix.csv")
cl <- read.csv("CAZy_2/outputData/Blackout_cellulase_matrix.csv")
xy <- read.csv("CAZy_2/outputData/Blackout_xylanase_matrix.csv")

enz <- do.call("cbind", list(aa, cl, xy))
enz <- enz[, !duplicated(names(enz))] # Remove duplicated columns


library(ggplot2); theme_set(theme_classic())
library(data.table)
df <- data.table(Sample = row.names(meta)
           , Treatment = trt
           , Cellulases = rowSums(cl)
           , Hemicellulases = rowSums(xy)
           , AAs = rowSums(aa))


BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')



p_caz <- ggplot(melt(df, id.vars = c("Sample", "Treatment"))
       , aes(x = Treatment, y = value, fill = Treatment)) + 
  geom_boxplot() +
  geom_point() +
  facet_wrap(~ variable) +
  scale_fill_manual(values = BlackoutPalette) +
  scale_y_continuous("Relative abundance (CPM)") +
  theme(legend.position = "none"
        , axis.text.x = element_text(size = 7))


df[Treatment == '1-Year\nGrassland', Time := 1]
df[Treatment == '10-Year\nGrassland', Time := 10]
df[Treatment == '1-Year\nBare', Time := -1]
df[Treatment == '10-Year\nBare', Time := -10]

# Model
mod_cel <- glm(Cellulases/1e6 ~ Time, data = df, family = glmmTMB::beta_family())
summary(mod_cel)

mod_hem <- glm(Hemicellulases/1e6 ~ Time, data = df, family = glmmTMB::beta_family())
summary(mod_hem)

mod_aa <- glm(AAs/1e6 ~ Time, data = df, family = glmmTMB::beta_family())
summary(mod_aa)


library(glmmTMB)

# Model
mod_cel <- glmmTMB(Cellulases/1e6 ~ Time, data = df, family = glmmTMB::beta_family())
summary(mod_cel)
drop1(mod_cel, test = 'Chisq')

mod_hem <- glmmTMB(Hemicellulases/1e6 ~ Time, data = df, family = glmmTMB::beta_family())
summary(mod_hem)
drop1(mod_hem, test = 'Chisq')

mod_aa <- glmmTMB(AAs/1e6 ~ Time, data = df, family = glmmTMB::beta_family())
summary(mod_aa)
drop1(mod_aa, test = 'Chisq')[2, 3]

ann_df <- data.table(variable = rep(c("Cellulases", "Hemicellulases", "AAs"), 2)
                     , value = rep(c(0.005, 0.001), each = 3)
                     , Intercept = c(mod_cel$sdr$par.fixed[1]
                                     , mod_hem$sdr$par.fixed[1]
                                     , mod_aa$sdr$par.fixed[1]
                                     , NA, NA, NA)
                     , Slope = c(mod_cel$sdr$par.fixed[2]
                                 , mod_hem$sdr$par.fixed[2]
                                 , mod_aa$sdr$par.fixed[2]
                                 , NA, NA, NA
                                 )
                     , Chi2 = c(NA, NA, NA
                                , drop1(mod_cel, test = 'Chisq')[2, 3]
                                , drop1(mod_hem, test = 'Chisq')[2, 3]
                                , drop1(mod_aa, test = 'Chisq')[2, 3])
                     , p = c(NA, NA, NA
                             , drop1(mod_cel, test = 'Chisq')[2, 4]
                             , drop1(mod_hem, test = 'Chisq')[2, 4]
                             , drop1(mod_aa, test = 'Chisq')[2, 4])
                     )
ann_df[1:3, Label := paste('logit(mu) ==', signif(Intercept, 3), signif(Slope, 2), ' %*% Time')]
ann_df[4:6, Label := paste("chi[1]^2 == ", signif(Chi2, 3), " * \",\"~~italic(p) == ", signif(p, 2))]


p_caz3 <- ggplot(melt(df, id.vars = c("Sample", "Treatment", "Time"))
                 , aes(x = Time #Treatment
                       , y = value/1e6
                       , fill = Treatment)) + 
  geom_smooth(aes(fill = NULL)
              #, se = FALSE
              , method = 'glm'
              , method.args=list(family=glmmTMB::beta_family)
  ) +
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values = BlackoutPalette) +
  theme(legend.position = "none"
        , panel.grid = element_blank()
        , axis.text.x = element_text(
          hjust = c(0, 1, 0, 1)
          , size = 6
        )
        , plot.margin = unit(c(0,0,0,0), "cm")) +
  facet_wrap(~ factor(variable, levels = c("Cellulases", "Hemicellulases", "AAs"))
             , nrow = 1) +
  scale_y_continuous("Relative abundance (CPM)"
                     #, expand = expansion(mult = c(0, .1))
                     #, limits = c(0, NA)
                     , labels = scales::label_percent(scale = 1e6
                                                      , suffix = "")
  ) +
  scale_x_continuous(name = 'Treatment'
                     , breaks = c(-10, -1, 1, 10)
                     , labels = c("10-Year\nBare"
                                  , "1-Year\nBare"
                                  , '1-Year\nGrassland'
                                  , '10-Year\nGrassland')
  ) +
  # equation lines
  geom_text(data = ann_df[1:3,]
            , aes(x = 0, fill = NULL
                  , label = Label)
            , size = 3
            , parse = TRUE
            ) +
  geom_text(data = ann_df[4:6,]
            , aes(x = 0, fill = NULL
                  , label = Label)
            , size = 3
            , parse = TRUE
  )
  p_caz3









# p_caz2 <- ggplot(melt(df, id.vars = c("Sample", "Treatment", "Time"))
#                 , aes(x = Time #Treatment
#                       , y = value/1e6
#                       , fill = Treatment)) + 
#   geom_smooth(aes(fill = NULL)
#               #, se = FALSE
#               , method = 'glm'
#               , method.args=list(family=glmmTMB::beta_family)
#   ) +
#   geom_boxplot() +
#   geom_point() +
#   #facet_wrap(~ variable) +
#   scale_fill_manual(values = BlackoutPalette) +
#   theme(legend.position = "none"
#       , panel.grid = element_blank()
#       , axis.text.x = element_text(
#         hjust = c(0, 1, 0, 1)
#         , size = 6
#       )
#       , plot.margin = unit(c(0,0,0,0), "cm")) +
#   facet_wrap(~ factor(variable, levels = c("Cellulases", "Hemicellulases", "AAs"))#, scales = "free_y"
#              , nrow = 1) +
#   scale_y_continuous("Relative abundance (CPM)"
#                      #, expand = expansion(mult = c(0, .1))
#                      #, limits = c(0, NA)
#                      , labels = scales::label_percent(scale = 1e6
#                                                       , suffix = "")
#   ) +
#   scale_x_continuous(name = 'Treatment'
#                      , breaks = c(-10, -1, 1, 10)
#                      , labels = c("10-Year\nBare"
#                                   , "1-Year\nBare"
#                                   , '1-Year\nGrassland'
#                                   , '10-Year\nGrassland')
#   ) +
#   annotate("text", x = 0, y = 0.02e-6, label = "logit(μ) == β[0] + β[1]*Time"# ~~~ Y %~% BetaBinomial(n, μ, φ)"
#   , parse = TRUE) +
#   #annotate("text", x = 0, y = 0.018e-6, label = "Y %~% BetaBinomial(n, μ, φ)"
#   #         , parse = TRUE) +
#   annotate("text", x = 0, y = 0.001e-6, label = "chi[1]^2 == 5.21 * \",\"~~italic(p) == 0.022"
#            , parse = TRUE)

p_caz2

kruskal.test(df$Cellulases, df$Treatment)
kruskal.test(df$Hemicellulases, df$Treatment)
kruskal.test(df$AAs, df$Treatment)
FSA::dunnTest(df$AAs ~ df$Treatment)


# Save a figure for SBB
tiff("Figures_2024/Review_updates/Figure_3_Lignocellulase_relative_abundance.tif"
     , width = 16.5, height = 17
     , units = "cm"
     , res = 600)
p_caz3
dev.off()












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

# tax_sum<-aggreg8(decostand(spc,"total"),tax,"Phylum")
# 
# 
# spc.s<-decostand(spc,"total")
# mod4<-metaMDS(spc.s)
# plot(mod4,"sites")
# ordispider(mod4,trt,label=T)
# ef<-envfit(mod4,tax_sum)
# plot(ef,p.max=0.005)
# 
# ef2<-envfit(mod4,meta[,46:ncol(meta)])
# plot(ef2,p.max=0.005,col="red")


#####enzymes
#enz<-read.csv("enz.csv") # I cbind-ed this in excel

aa <- read.csv("CAZy_2/outputData/Blackout_aa_matrix.csv")
cl <- read.csv("CAZy_2/outputData/Blackout_cellulase_matrix.csv")
xy <- read.csv("CAZy_2/outputData/Blackout_xylanase_matrix.csv")

enz <- do.call("cbind", list(aa, cl, xy))
enz <- enz[, !duplicated(names(enz))] # Remove duplicated columns

mod5<-metaMDS(enz)
plot(mod5,"sites")
ordispider(mod5,trt,label=T)
library(devtools)
source_gist("https://gist.github.com/robiwangriff/79633a738b128e64226bf58381232da7")

library(labdsv)
inds<-indies2df(indvals= indval(enz,trt),taxa=data.frame(row.names=names(enz)),return.alltax=TRUE)

write.csv(inds, "CAZy_2/outputData/Blackout_CAZyme_indvals.csv")

inds2<-indies2df(indvals = indval(enz, sub("_.*", "", trt))
                 ,taxa=data.frame(row.names=names(enz))
                 ,return.alltax=TRUE)

inds$Row.names %in% inds2$Row.names
inds2$Row.names %in% inds$Row.names

inds2$Row.names == inds$Row.names
inds2$group_name == sub("_.*", "", inds$group_name)
itest <- data.table(CAZyme = inds$Row.names, full = inds$group_name, simple = inds2$group_name)
itest[!is.na(full) & !is.na(simple),]
itest[!is.na(simple),]

enz
trt
inds3<-indies2df(indvals = indval(enz[trt!='bare_new'], sub("_.*", "", trt[trt!='bare_new']))
                 ,taxa=data.frame(row.names=names(enz[trt!='bare_new']))
                 ,return.alltax=TRUE)

inds3
#write.csv(inds, "CAZy_2/outputData/Blackout_CAZyme_indvals.csv")


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
spc
spc$CAZyme <- sub("\\.", "|", spc$CAZyme)

range(spc$NMDS1)



unique(blc$Treatment)

trtlevs <- c("bare_old",    "bare_new",  "control_new",  "control_old" )

Community_comp <- ggplot(blc, aes(x = NMDS1, y = NMDS2, colour = Treatment)) +
  geom_point(size = 3) +
  cowplot::theme_cowplot() + 
  geom_text(data = spc
            , aes(label = CAZyme
                  , colour = names(colors.spc))
            , size = 3) +
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
  scale_y_continuous("MDS2") +
  scale_x_continuous("MDS1", limits = c(-0.55, 0.55))

Community_comp

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
en_coord_cont = as.data.frame(vegan::scores(fit, "vectors")) * ordiArrowMul(fit) * 1.8
en_coord_cont$Treatment <- "10-Year\nBare"

row.names(en_coord_cont) <- sub("Percentage_", "", row.names(en_coord_cont))
row.names(en_coord_cont) <- sub("_", " ", row.names(en_coord_cont))
library(stringr)
row.names(en_coord_cont) <- stringr::str_to_sentence(row.names(en_coord_cont))
row.names(en_coord_cont)[row.names(en_coord_cont) == "Nc"] <- "N:C"
en_coord_cont$pval <- p.adjust(fit$vectors$pvals, method = "fdr") 
en_coord_cont$Linetype <- ifelse(en_coord_cont$pval < 0.05, yes = 'p < 0.05', no = 'p < 0.1')

Community_comp <- Community_comp +
  geom_segment(data = en_coord_cont
               , aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2
                     , linetype = Linetype)
               , size = 1, alpha = 0.5, colour = "grey30"
               , arrow = arrow(ends = "last", length = unit(0.1, "inches"))
  ) +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), 
            label = row.names(en_coord_cont)
            , colour = "navy", fontface = "bold") +
  annotate(geom = 'text', x = -0.45, y = 0.35
           , label = paste0("Stress = ", round(mod5$stress, 2))
           , colour = "black")

Community_comp

# Save a figure for SBB
tiff("Figures_2024/Figure_4_Lignocellulase_community_composition_Rob.tif"
     , width = 16.5, height = 17
     , units = "cm"
     , res = 600)
Community_comp
dev.off()