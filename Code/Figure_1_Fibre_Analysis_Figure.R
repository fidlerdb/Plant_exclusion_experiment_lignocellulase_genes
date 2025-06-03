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
ff <- readRDS("Fibre_Analysis/Output_Data/Fibre_Data_For_Figure_2021-05-05.rds")
ff

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

#### Linear time scale ####

ff_d

ff_d[Treatment == '1-Year\nGrassland', Time := 1]
ff_d[Treatment == '10-Year\nGrassland', Time := 10]
ff_d[Treatment == '1-Year\nBare', Time := -1]
ff_d[Treatment == '10-Year\nBare', Time := -10]

#### Models for line equations ####

library(glmmTMB)

# Model
mod_cel <- glmmTMB(Abundance/100 ~ Time
                   , data = ff_d[Polymer == 'Cellulose', ]
                   , family = glmmTMB::beta_family())

summary(mod_cel)
drop1(mod_cel, test = 'Chisq')

mod_hem <- glmmTMB(Abundance/100 ~ Time
                   , data = ff_d[Polymer == 'Hemicellulose', ]
                   , family = glmmTMB::beta_family())

summary(mod_cel)
drop1(mod_cel, test = 'Chisq')

mod_aa <- glmmTMB(Abundance/100 ~ Time
                   , data = ff_d[Polymer == 'Lignin', ]
                   , family = glmmTMB::beta_family())

summary(mod_aa)
drop1(mod_aa, test = 'Chisq')

ann_df <- data.table(Polymer = rep(c("Cellulose", "Hemicellulose", "Lignin"), 2)
                     , Abundance = c(40, 45, 8, 35, 40, 7) 
                       #rep(c(40, 35), each = 3)
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
ann_df[1:3, Label := paste('logit(mu) ==', signif(Intercept, 3), " + ",signif(Slope, 2), ' %*% Time')]
ann_df[4:6, Label := paste("chi[1]^2 == ", signif(Chi2, 3), " * \",\"~~italic(p) == ", signif(p, 2))]



#### Fibre analysis plot to add letters to  ####

pf <- ggplot(ff_d, aes(x = Time  #Treatment
                       , y = (Abundance * 10)/1e3
                       , fill = Treatment)) + 
  geom_smooth(aes(fill = NULL)
              #, se = FALSE
              , method = 'glm'
              , method.args=list(family=glmmTMB::beta_family)
              ) +
  geom_boxplot(outliers = FALSE) +
  #stat_summary(fun = "median", geom = "bar") +
  geom_point(position = position_jitter(0.15)
             , size = 0.1) +
  ylab("Polymer relative abundance (g/kg dry soil)" #"Grams per kg dry soil"
       ) + 
  scale_fill_manual(values = BlackoutPalette
                     , breaks = levels(ff_d$Treatment)) +
  theme_bw() +
  theme(legend.position = "none"
        , panel.grid = element_blank()
        , axis.text.x = element_text(
          hjust = c(0, 1, 0, 1)
          , size = 7.5
          )
        , plot.margin = unit(c(0,0,0,0), "cm")
        ) +
  facet_wrap(~ Polymer, scales = "free_y", nrow = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))
                     , limits = c(0, NA)
                     , labels = scales::label_percent(scale = 1e3
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

pf

# Save the plot
tiff(filename = 'Figures_2024/Review_updates/Figure_1_Fibre_Analysis.tif'
     , width = 3, height = 6
     , units = 'in'#"cm"
     , res = 600)
pf 
dev.off()

