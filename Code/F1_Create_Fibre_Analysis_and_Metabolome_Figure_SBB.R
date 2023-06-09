rm(list = ls())

library(data.table)
library(ggplot2)
library(cowplot)
library(ggnewscale)
library(gridExtra)

#### Data import ####

# Colour palette for the treatment variable
BlackoutPalette <- c('#7b3294'
                     ,'#c2a5cf'
                     ,'#a6dba0'
                     ,'#008837')

# Function for 95% confidence intervals
mean_ci <- function(x){
  mean_se(x, mult = 1.96)
}

# Fibre analysis data
ff <- readRDS("Data/Fibre_Data_For_Figure_2021-05-05.rds")
ff

# Metabolite data
mf <- readRDS("Data/Metabolite_Data_For_Figure.rds")
mf

#### Data cleaning ####

## Fibre analysis data
ff

ff_d <- ff[[1]]
colsToKeep <- c("Sample", "Treatment","Proportion_Cellulose","Proportion_Hemicellulose", "Proportion_Lignin")
ff_d <- ff_d[,..colsToKeep]
setnames(ff_d, c("Sample", "Treatment","Cellulose","Hemicellulose", "Lignin"))
ff_d[, Measurement := "Fibre Analysis"]

ff_d <- melt(ff_d, id = c("Sample", "Treatment", "Measurement"), variable.name = "Polymer", value.name = "Abundance")
ff_d[, Abundance := Abundance * 100]

ff_d[,BinBase.name := NA] # For compatibility with the metabolite data

# Simple plot of the data in this format
ggplot(ff_d, aes(x = Treatment, y = Abundance)) + geom_point() + facet_grid(Polymer~Measurement
                                                                            #, nrow = 3
                                                                            , scales = "free_y"
                                                                          
                                                                            )
levels(ff_d$Treatment) <- sub(" ", "\n" , levels(ff_d$Treatment)) # Two lines for plotting

## Metabolite data
mf

mf_d <- rbindlist(mf[1:3], idcol = "Polymer")

# Sort out a polymer variable
mf_d$Polymer <- as.character(mf_d$Polymer)
mf_d[Polymer == '1', Polymer := "Cellulose"]
mf_d[Polymer == '2', Polymer := "Hemicellulose"]
mf_d[Polymer == '3', Polymer := "Lignin"]

# Sort out the treatment variable
mf_d$Treatment
levels(mf_d$Treatment) <- c("1-Year\nBare", "10-Year\nBare", "1-Year\nGrassland", "10-Year\nGrassland")
mf_d$Treatment <- relevel(mf_d$Treatment, ref = "10-Year\nBare")

setnames(mf_d, c("Polymer", "Chemical", "Sample", "RawAbundance", "Treatment",  "Abundance"))
# Add a Measurement column
mf_d[, Measurement := "Metabolome"]

mf_d[, Abundance := scale(RawAbundance, center = TRUE, scale = TRUE), by = Chemical]

mf_d[Chemical == "4-hydroxybenzoic acid", Chemical := "4-hydroxybenzoic\nacid"]
mf_d[Chemical == "3,6-anhydro-D-galactose", Chemical:= "3,6-anhydro-D-\ngalactose"]

#### Fibre analysis plot to add letters to  ####

pf <- ggplot(ff_d, aes(x = Treatment, y = Abundance
                       , fill = Treatment)) + 
  stat_summary(fun = "median", geom = "bar") +
  geom_point(position = position_jitter(0.15)) +
  ylab("Percentage of biomass") + 
  scale_fill_manual(values = BlackoutPalette
                     , breaks = levels(ff_d$Treatment)) +
  theme_bw() +
  theme(legend.position = "none"
        , panel.grid = element_blank()) +
  facet_grid(Polymer~Measurement, scales = "free_y")

#### Metabolite plot to add letters to ####

pm <- ggplot(mf_d, aesthetics = aes(x = Treatment, y = Abundance)) + 
  # Data
  geom_point(data = mf_d[Polymer == "Cellulose",]
             , aes(x = Treatment, y = Abundance, shape = Chemical)
             , position = position_jitter(0.05)
             , alpha = 0.5
  ) + 
  new_scale("shape") +
  geom_point(data = mf_d[Polymer == "Hemicellulose",]
             , aes(x = Treatment, y = Abundance, shape = Chemical)
             , position = position_jitter(0.05)
             , alpha = 0.5
  ) +
  new_scale("shape") +
  geom_point(data = mf_d[Polymer == "Lignin",]
             , aes(x = Treatment, y = Abundance, shape = Chemical)
             , position = position_jitter(0.05)
             , alpha = 0.5
  ) +
  
  # Mean values
  stat_summary(fun.data = 'mean_ci'
               , aes(x = Treatment, y = Abundance, colour = Treatment), size = 1.1) +
  scale_color_manual(values = BlackoutPalette
                     , breaks = levels(mf_d$Treatment)) +
  
  # Aesthetics
  facet_grid(Polymer ~ Measurement) +
  theme_bw() +
  theme(panel.grid = element_blank()
        , legend.text = element_text(size = 7)) +
  ylab("Standardized chemical abundance")  +
  theme(#legend.spacing.y = unit(0.5, 'cm')
        #,
        legend.text.align = 0
        , legend.box.just = 'left'
        , legend.box.margin = unit(c(0.1, 1, 0.1, 0), units = "cm")
    )  +
  ## important additional element
  guides(colour = "none"
         #, fill = guide_legend(byrow = TRUE)
         )

pm

pl <- get_legend(pm)
#pl <- pl + theme(text = element_text(size = 7))
pm <- pm + theme(legend.position = "none")

#### Plot legend ####

# Put the grobs in the correct order, also spaces them nicely
# Cellulose, hemicellulose lignin

pl2 <- grid.arrange(grobs = list(pl$grobs[[3]], pl$grobs[[1]], pl$grobs[[2]])
                    , ncol = 1
)
grid.arrange(pl2)

#### Add significance letters ####

format_Letters <- function(x, polymer){
  setnames(x, c("Treatment", "Letters", "Abundance"))
  x[,Polymer := polymer]
  x[, Treatment := sub(" ", "\n", Treatment)]
  return(x)
}

fib_cel_letters <- format_Letters(ff[[2]], "Cellulose")
fib_hem_letters <- format_Letters(ff[[3]], "Hemicellulose")
fib_lig_letters <- format_Letters(ff[[4]], "Lignin")
fib_lig_letters$Abundance <- 12

pf2 <- pf + geom_text(data = fib_cel_letters, aes(x = Treatment, y = Abundance, label = Letters), colour = 'black') +
  geom_text(data = fib_hem_letters, aes(x = Treatment, y = Abundance, label = Letters), colour = 'black') +
  geom_text(data = fib_lig_letters, aes(x = Treatment, y = Abundance, label = Letters), colour = 'black') #+
  # ylim(-1,100)

pm

mf[[4]] # Cel
mf[[5]] # Hem
mf[[6]] # Lig

format_letters_met <- function(x, polymer){
  setDT(x)
  x[,Polymer := polymer]
  return(x)
}

met_cel_letters <- format_letters_met(mf[[4]], "Cellulose")
met_hem_letters <- format_letters_met(mf[[5]], "Hemicellulose")
met_lig_letters <- format_letters_met(mf[[6]], "Lignin")

met_cel_letters$Letters <- c("b", "a", "ab", "ab") # Re-jig these so that they work the normal way round
met_lig_letters$Letters <- c("b", "a", "c", "ac") # Re-jig these so that they work the normal way round

met_cel_letters$NormAbund <- 3.2
met_hem_letters$NormAbund <- 3.2
met_lig_letters$NormAbund <- 3.2

pm2 <- pm + geom_text(data = met_cel_letters, aes(x = Treatment, y = NormAbund, label = Letters)) +
  geom_text(data = met_hem_letters, aes(x = Treatment, y = NormAbund, label = Letters)) + 
  geom_text(data = met_lig_letters, aes(x = Treatment, y = NormAbund, label = Letters))

#### Create the finalised plot ####

Figure2 <- plot_grid(pf2 + theme(axis.text.x = element_text(size = 6)
                                 , plot.margin = unit(c(0,0,0,0.25), units = "cm"))
                     , pm2 + theme(axis.text.x = element_text(size = 6)
                                   , plot.margin = unit(c(0,0,0,0.25), units = "cm"))
                     , pl2 
                     , ncol = 3, rel_widths = c(1,1,0.46)
                     )

tiff(filename = 'Figures_2023_SBB/Figure_1_Fibre_Analysis_and_Metabolomics.tif'
     , width = 15.5, height = 17
     , units = "cm"
     , res = 600)
Figure2 
dev.off()
