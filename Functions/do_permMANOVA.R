do_permMANOVA <- function(community_matrix, Treatment){

species_matrix_rel <- decostand(community_matrix, method = "frequency")

nmds1 <- metaMDS(species_matrix_rel)
# Shepards test/goodness of fit
stressplot_out <- stressplot(nmds1) # Produces a Shepards diagram
nmds_obj <- data.table(nmds1$points)

# Create a plot
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
  theme(text = element_text(size = 16))# +

### Check for significant differences in  microbial community composition

# Check PermMNOVA assumptions 
Group_disp <- betadisper(vegdist(species_matrix_rel, method = "bray")
                         , group = Treatment
                         , type = "median")


# Do a permMANOVA
PermMANOVA <- adonis2(species_matrix_rel ~ Treatment
                      , by = "margin" # Assess marginal significance of terms
                      , parallel = 3 # Make it faster
                      , permutations = 60000)


return(list(Stressplot = stressplot_out
       , NMDS = Community_comp
       , Group_disp = Group_disp
       , PermMANOVA = PermMANOVA))
}
