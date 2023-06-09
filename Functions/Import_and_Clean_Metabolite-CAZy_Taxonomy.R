# Read in the taxonomy data--normalized
df <- fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Henfaes_dbCAN_Taxonomy_ReadAbundance.csv")

# Remove viral contigs
df <- df[df$SuperKingdom != "Viruses",]
# Remove contigs which have no phylum-level assignment
df <- df[!is.na(df$Phylum),]

# Read in the CAzyme data
cl <- fread("CAZy/04_CAZyme_Origin_Analysis/Cleaned_Data/Positively_Correlated_CAZymes-Metabolites.csv")

# Sort out the weird ones
cl[cl$CAZymes == "GH1.3.2.1..",]$CAZymes <- "GH1"
cl[cl$CAZymes == "GH10.3.2.1.8",]$CAZymes <- "GH10_endo-1,4-β-xylanase"
cl[cl$CAZymes == "GH3.3.2.1.21",]$CAZymes <- "GH3_β-Glucosidase"

# Create useful vectors for building a matrix of CAZymes vs Orders/Phyla
AllCazymes <- unique(cl$CAZymes)
AllPhyla <- unique(df$Phylum)
AllMetabolites <- unique(cl$Metabolite)

# All metabolite-cazyme combinations
Met_Caz <- paste(rep(AllMetabolites, each = length(AllCazymes)), 
                 rep(AllCazymes, 7)
                 , sep = "-"
)

# Create an empty matrix the correct size
CRM <- matrix(NA
              , nrow = length(AllPhyla)
              , ncol = length(Met_Caz))
CRM <- data.frame(CRM)
names(CRM) <- Met_Caz
rownames(CRM) <- AllPhyla