source("https://bioconductor.org/biocLite.R")  
library("tidyverse")
library("phyloseq")
library("QUBIC")
library("metagenomeSeq")
library("spatstat") # "im" function 
library("pvclust")
library("colorspace") # get nice colors
library("dendextend")
library("dbscan")

load("data/mothur_phyloseq.RData")
abundantTaxa = filter_taxa(mothur, function(x) sum(x == 0) <= 4, TRUE)
basedOnGenus <- as.data.frame(tax_table(abundantTaxa)) %>% 
  filter(!str_detect(Genus, 'uncultured'), !str_detect(Genus, 'unclassified'))

knownTaxa = subset_taxa(abundantTaxa, Genus %in% basedOnGenus$Genus)
basedOnphylums <- as.data.frame(tax_table(knownTaxa)) %>% 
  group_by(Phylum) %>% 
  count() %>% 
  filter( n > 5)

workingTaxa = subset_taxa(knownTaxa, Phylum %in% basedOnphylums$Phylum)
workingTaxa <- prune_taxa(taxa_sums(knownTaxa) > 5, knownTaxa)

#######################################################

dev.off()
graphics.off()
#par("mar")
#par(mar=c(1,1,1,1))
par(mar = rep(0,4))

theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

######################################################

wt <- prune_taxa(names(sort(taxa_sums(workingTaxa), decreasing = TRUE)),
                 workingTaxa)
# sum-normalize
wt.norm <- transform_sample_counts( wt, function(x) x/sum(x) )
# generate stacked bar plot

p = plot_bar(wt.norm, "Depth_m", fill="Phylum")
p + geom_bar(aes(color = Phylum, fill = Phylum), stat = 'identity', position = 'stack')


#####################################

## Perform correlations
D <-  phyloseq_to_metagenomeSeq(workingTaxa)
# agg <- aggregateByTaxonomy(D, lvl="Phylum", norm=FALSE, aggfun=colSums)
dim(D)
mat = MRcounts(D)
cors = correlationTest(D, norm=TRUE, log=FALSE)
ind  = correctIndices(nrow(mat))
cormat = as.matrix(dist(x = mat, diag = FALSE))
#cormat[cormat > 0] = 0
cormat[upper.tri(cormat)][ind] = round(cors[, 1], 4)
cormat[lower.tri(cormat)][ind] = round(cors[, 1], 4)
#cormat <- round(1000 * cormat) # (prints more nicely)

## Perform hclust
bthlcust <- pvclust(cormat, method.dist="correlation",
                    nboot=1000, parallel=TRUE)
plot(bthlcust)
pvrect(bthlcust)

dend <- as.dendrogram(bthlcust)
# Color the branches based on the clusters:
dend <- color_branches(dend)
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend, hang_height=0.1)
# cut the dendrogram based on height:
OTUs = cutree(dend = dend, k = 100, FUN = labels)
for(otu in OTUs)
{
  print(as.data.frame(tax_table(prune_taxa(otu, workingTaxa)))["Genus"])
  n <- readline(prompt="Enter an integer: ")
}

pt = prune_taxa(OTUs[[1]], workingTaxa)
plot_richness(pt, x="O2_uM")

heatmap(cormat)

# And plot:
# circlize_dendrogram(dend)

ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
lm.fit (group ~ weight)

###############

rarefy_even_depth(subGenTaxa, sample.size = 500)

estimate_richness(workingTaxa, split = TRUE, measures = NULL)

p = plot_richness(workingTaxa, x="Depth_m", color="Depth_m")
p + geom_point(size=5, alpha=0.7)

plot_richness(workingTaxa, x="O2_uM")
