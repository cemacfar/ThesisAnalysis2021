#import the 16S dataset
setwd("~/Bioinformatics")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("DESeq2")
BiocManager::install("decontam")
#BiocManager::install("ShortRead")
BiocManager::install(version='devel')
BiocManager::install("microbiome")

install.packages("readxl")
install.packages("tidyverse")
install.packages("devtools")

library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

library(phyloseq) #needed moving forward for statistical testing
#library(DESeq2) #will be used to do differential abundance
library(decontam) #gets rid of the blanks
library(readxl) #read in the datasheet
library(tidyverse) #needed to help organize the howler metadata
#library(ShortRead) #needed to read in the Fasta files

#decontaminate the samples first
Howler_Metadata16s = read_excel("Howler Metadata.xlsx",
                                sheet = "16s")
asv_fasta_dada = readFasta("dadaspASVs.fa")
asv_tax_dada = read.table("dadaspASVs_taxonomy.txt",
                          header = T,
                          row.names = 1,
                          check.names = F,
                          sep = "\t")
asv_tab_dada = read.table("dadaspASVs_counts.txt",
                          header = T,
                          row.names = 1,
                          check.names = F,
                          sep = "\t")

colnames(asv_tab_dada)
#blanks are last samples, 74 to 84
#logical vector for last 11 samples
decontam_logic_dada = c(rep(TRUE, 73), rep(FALSE, 11))

contam_df_dada = isContaminant(t(asv_tab_dada), neg = decontam_logic_dada)

table(contam_df_dada$contaminant)
contam_asvs_dada = row.names(contam_df_dada[contam_df_dada$contaminant == TRUE,])

# making new count table
asv_tab_no_contam_dada = asv_tab_dada[!row.names(asv_tab_dada) %in% contam_asvs_dada, ]
# making new taxonomy table
asv_tax_no_contam_dada = asv_tax_dada[!row.names(asv_tax_dada) %in% contam_asvs_dada, ]




#phyloseq object
phyotu16 = otu_table(asv_tab_no_contam_dada, taxa_are_rows = TRUE)
phytax16 = tax_table(as.matrix(asv_tax_no_contam_dada))
phytest16 = phyloseq(phyotu16, phytax16)
phytest16

howmeta.df16 = Howler_Metadata16s %>% select(-sample) %>% as.data.frame
rownames(howmeta.df16) = Howler_Metadata16s$sample
howler.data16 = sample_data(howmeta.df16, errorIfNULL = F)

setdiff(colnames(phyotu16), rownames(howler.data16))
#no differences in names, should merge together 

phytestsamp16 = merge_phyloseq(phytest16, howler.data16)
phytestsamp16

#rename for ease moving forward
phy16sbl = phytestsamp16

#take out blank samples
phy16s = prune_samples(sample_names(phy16sbl)[1:73], phy16sbl)
sample_names(phy16s)


############################
#a/b tests bacteria location
############################
library(ggplot2)
#alpha diversity, richness esitmate by sample plot_richness() function on our phyloseq object
#CHAO1 is RICHNESS and SHANNON is DIVERSITY
Frozsamp = subset_samples(phy16s, storage=="Liquid Nitrogen")
Costasamp = subset_samples(phy16s, field_site != "Monkey River")

#by grouping rather than all samples
#visualize without tranformed data
library(viridis)
library(wesanderson)
plot_richness(Costasamp, x="field_site" ,color="storage", measures=c("Chao1", "Shannon"), title = "Howler Monkey 16s alpha diversity storage") + 
  theme(legend.title = element_blank()) +
  geom_boxplot()
plot_richness(Costasamp, x="field_site", measures=c("Chao1", "Shannon"), title = "Howler Monkey 16s alpha diversity") + 
  theme(legend.title = element_blank()) +
  geom_boxplot()
#final figures
sample_data(Costasamp)$field_site = factor(sample_data(Costasamp)$field_site, levels = c("Santa Rosa","Eladio's Ranch","Tamarindo"))
Allhabitat = plot_richness(Costasamp, 
              x="field_site", 
              color = "field_site",
              measures=c("Chao1","Shannon", "Simpson")) + 
  theme(legend.title = element_blank()) +
  geom_boxplot() +
  labs(x = "Habitat") +
  scale_color_manual(values = wes_palette("Darjeeling2"))

Fragm = plot_richness(Costasamp,
              x="fragmentation",
              color = "fragmentation",
              measures = c("Chao1", "Shannon", "Simpson")) +
  theme(legend.title = element_blank()) +
  geom_boxplot() +
  labs(x = "Habitat") +
  scale_color_manual(values = wes_palette("Darjeeling2"))
library(ggpubr)
figurealphacombo = ggarrange(Allhabitat, Fragm, nrow = 2, labels = "AUTO")

alphtest = estimate_richness(Costasamp)
#Shannon is diversity. Sig diff in diversity of bacteria, taking into account the richness and evenness
pairwise.wilcox.test(alphtest$Shannon, sample_data(Costasamp)$field_site) 
#CHAO1 is richness. Sig diff in richness of bacteria, none
#Simpson measures what shannon does but places less emphasis on rare species, they will have less impact
pairwise.wilcox.test(alphtest$Chao1, sample_data(Costasamp)$field_site)
#16s beta diversity clustering
pairwise.wilcox.test(alphtest$Chao1, sample_data(Costasamp)$fragmentation)
pairwise.wilcox.test(alphtest$Shannon, sample_data(Costasamp)$fragmentation)

#TRANSFORM NOW
#visualize with the transformed data, could get rid of the proportion function and transform diff

library(microbiome)
phy16log = transform(Costasamp, transform = "log10")
phy16hell = transform(Costasamp, transform = "hellinger")
phy16sprop2 = transform_sample_counts(Costasamp, function(otu) otu/sum(otu))

ord.nmds.bray161 = ordinate(phy16sprop2, method="NMDS", distance="bray")
plot_ordination(phy16sprop2, ord.nmds.bray161, color="field_site", title="16s Bray NMDS Costa Rica proportion")

ord.nmds.bray16log = ordinate(phy16log, method="NMDS", distance="bray")
plot_ordination(phy16log, ord.nmds.bray16log, color="field_site", title="16s Bray NMDS Costa Rica log trans")

ord.nmds.bray16hell = ordinate(phy16hell, method="NMDS", distance="bray")
plot_ordination(phy16hell, ord.nmds.bray16hell, color="field_site", title="16s Bray NMDS Costa Rica hellinger trans")





ord.nmds.bray16log = ordinate(phy16log, method="PCoA", distance="bray")
plot_ordination(phy16log, ord.nmds.bray16log, color="field_site") +
  scale_color_manual(values = wes_palette("Darjeeling2")) +
  labs(color = "Habitat") +
  stat_ellipse()

ord.nmds.bray16log = ordinate(phy16log, method="PCoA", distance="bray")
plot_ordination(phy16log, ord.nmds.bray16log, color="fragmentation") +
  scale_color_manual(values = wes_palette("Darjeeling2")) +
  labs(color = "Habitat")

#SIGNIFICANCE TESTING#

#alpha significance
#no transformation of the data here, just use the Costa Rica samples used in the visualization
alpha = estimate_richness(Costasamp)
      #Shannon is diversity. Sig diff in diversity of bacteria, taking into account the richness and evenness
pairwise.wilcox.test(alpha$Shannon, sample_data(Costasamp)$field_site) #santa and tam 0.03, santa and hela 0.0005, santa and tam 0.6
      #CHAO1 is richness. Sig diff in richness of bacteria, none
      #Simpson measures what shannon does but places less emphasis on rare species, they will have less impact
pairwise.wilcox.test(alpha$Simpson, sample_data(Costasamp)$field_site)

#beta significance
#use the transformed data here, have to deal with the low values and zeros somehow
#betadisper to check for homogen. dispers.
library(vegan)

#TESTING OTHER WAY NO DESEQ
phydist = phyloseq::distance(phy16log, "bray")
df = as(sample_data(phy16log), "data.frame")
groups = df[["field_site"]]
groups

mod = betadisper(phydist, groups)
anova(mod)
mod.HSD = TukeyHSD(mod)

plot(mod)
boxplot(mod)
plot(mod.HSD)

mod2 = betadisper(phydist, groups, bias.adjust = T)
anova(mod2)
mod.HSD2 = TukeyHSD(mod2)
plot(mod.HSD2)

adonis(phydist ~ groups, perm = 999)

df

pairwise.adonis(phydist, groups)
