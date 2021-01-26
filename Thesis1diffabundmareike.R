#new diff abund figure from Mareike's code

#subset data that is only in the ethanol samples
CoprologicalEthanol = subset_samples(transf, storage=="Ethanol")
#conglomerate down to the phylum level for ease
##ce.phylum = tax_glom(CoprologicalEthanol, taxrank = "Genus", NArm=F)
#convert to DESeq2 
ethanol.helminth = phyloseq_to_deseq2(CoprologicalEthanol, design = ~Helminth_Infection)
#function to avoid error with 0s in case of a a high prevalence of sparsely sampled OTUs 
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
#apply function to deseq object (replaces 0s with 1s, I think, mareike's function)
geoMeans = apply(counts(ethanol.helminth), 1, gm_mean)
ethanol.helminth = estimateSizeFactors(ethanol.helminth, geoMeans = geoMeans)

# run DESeq function to compare groups
ethanol.helminth = DESeq(ethanol.helminth, test="Wald", fitType="local") 
#results
res.ethhelm = results(ethanol.helminth, cooksCutoff = T)
alpha = 0.01
test = na.omit(res.ethhelm)
sigtab.eth = res.ethhelm[which(res.ethhelm$padj < alpha), ]
sigtab.eth = cbind(as(sigtab.eth, "data.frame"), 
                   as(tax_table(CoprologicalEthanol)[rownames(sigtab.eth), ], "matrix"))
dim(sigtab.eth)
mcols(res.ethhelm, use.names = TRUE)
# Phylum order - order from highest to lowest log2fold change for visual appeal
x = tapply(sigtab.eth$log2FoldChange, sigtab.eth$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab.eth$Phylum = factor(as.character(sigtab.eth$Order), levels=names(x))
# Genus order
x = tapply(sigtab.eth$log2FoldChange, sigtab.eth$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab.eth$Genus = factor(as.character(sigtab.eth$Genus), levels=names(x))

# Make the figure
p_DiffExp_helm <- ggplot(sigtab.eth, aes(x=Family, y=log2FoldChange, color=Order, size = baseMean)) + 
  geom_point() + 
  geom_hline(yintercept = 0, linetype="dotted", alpha=0.5) +
  theme(axis.text.x = element_text(angle = -50, hjust = 0, vjust=0.5)) +
  scale_color_viridis(discrete = T); p_DiffExp_helm

test = na.omit(res.ethhelm)
difftax162 = cbind(as(test, "data.frame"), as(tax_table(CoprologicalEthanol)[rownames(test), ], "matrix"))
difftax162

p_DiffExp_helm_ALL = ggplot(difftax162, aes(x=Class, y=log2FoldChange, color=Phylum)) +
  geom_jitter(size = 1, width = 0.2) +
  theme(axis.text.x = element_text(angle = -90)) +
  scale_colour_viridis(discrete = T); p_DiffExp_helm_ALL

library(ggpubr)
figure4combo = ggarrange(p_DiffExp_helm_ALL, p_DiffExp_helm, nrow = 2, labels = "AUTO")





