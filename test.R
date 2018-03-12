source("https://bioconductor.org/biocLite.R")  
library("tidyverse")
library("phyloseq")
library("ggplot2")
library("QUBIC")
library("metagenomeSeq")
library("dplyr")
library("spatstat") # "im" function 
library("pvclust")
library("colorspace") # get nice colors
library("dendextend")

load("data/mothur_phyloseq.RData")
abundantTaxa = filter_taxa(mothur, function(x) sum(x == 0) <= 4, TRUE)
basedOnGenus <- as.data.frame(tax_table(abundantTaxa)) %>% 
  filter(!str_detect(Genus, 'uncultured'), !str_detect(Genus, 'unclassified')) %>%
  select(Genus)

knownTaxa = subset_taxa(abundantTaxa, Genus %in% basedOnGenus$Genus)
basedOnphylums <- as.data.frame(tax_table(knownTaxa)) %>% 
  group_by(Phylum) %>% 
  count() %>% 
  filter( n > 5)

workingTaxa = subset_taxa(knownTaxa, Phylum %in% basedOnphylums$Phylum)
workingTaxa <- prune_taxa(taxa_sums(knownTaxa) > 5, knownTaxa)

#############################################################

## Perform correlations
D <-  phyloseq_to_metagenomeSeq(workingTaxa)
# agg <- aggregateByTaxonomy(D, lvl="Phylum", norm=FALSE, aggfun=colSums)
dim(D)
mat = MRcounts(D)
cors = correlationTest(D, norm=TRUE, log=FALSE)
ind  = correctIndices(nrow(mat))
cormat = as.matrix(dist(x = mat, diag = FALSE))
cormat[cormat > 0] = 0
cormat[upper.tri(cormat)][ind] = round(cors[, 1], 4)
#cormat[lower.tri(cormat)][ind] = round(cors[, 1], 4)
#cormat <- round(1000 * cormat) # (prints more nicely)

## Perform hclust
dev.off()
graphics.off()
#par("mar")
#par(mar=c(1,1,1,1))
par(mar = rep(0,4))
bthlcust <- pvclust(workingTaxa@otu_table, method.dist="correlation",nboot=1000, parallel=FALSE)
#plot(bthlcust)
#pvrect(bthlcust)

dend <- as.dendrogram(bthlcust)
# Color the branches based on the clusters:
dend <- color_branches(dend)
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend, hang_height=0.1)
# cut the dendrogram based on height:
OTUs = cut_lower_fun(dend = dend, h = 0.2, FUN = labels)
pt = prune_taxa(OTUs[[2]], workingTaxa)
plot_richness(pt, x="O2_uM")

heatmap(cormat)

# And plot:
# circlize_dendrogram(dend)


###############

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

estimate_richness(workingTaxa, split = TRUE, measures = NULL)

p = plot_richness(workingTaxa, x="Depth_m", color="Depth_m")
p + geom_point(size=5, alpha=0.7)

plot_richness(workingTaxa, x="O2_uM")
