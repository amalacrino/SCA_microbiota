# Characterization  of the Sugarcane Aphid (*Melanaphis sacchari* Zehntner) Microbiota in the Continental US

### *Holt JR, Styer A, White J, Armstrong JS, Nibouche S, Costet L, Malacrinò A, Antwi JB, Wulff J, Peterson G,  McLaren N, Medina RF*

#### Journal, 20xx. DOI: 

This is code to replicate the analyses reported in our manuscript. Code developed by Antonino Malacrinò (antonino.malacrino 'at' gmail 'dot' com) | 2019



### Abstract 





### Library Description

DNA samples were sent to the Molecular Research DNA Lab (MR. DNA, Shallowater, TX, USA) for metabarcoding analysis targeting the V3-V4 16S rRNA bacterial region. Samples were sequenced on an Illumina MiSeq platform (Illumina, San Diego, CA, USA) using the MiSeq Reagent Kit v3 300PE chemistry. 

Data handling was carried out using QIIME 1.9 (Caporaso et al. 2010, Caporaso et al. 2012), quality-filtering reads (Phred >= 25), binning OTUs using open-reference OTU-picking through UCLUST algorithm, and discarding chimeric sequences discovered with USEARCH 6.1 (Edgar 2010). Taxonomy was assigned to each OTU through the BLAST method using Greengenes database for 16S rRNA (Caporaso et al. 2012). 



### Data analysis

#### Data

Raw data is available at NCBI SRA under Bioproject PRJNA578411

#### Metadata

| Sample_ID       | Sample_type | Location     | Host_plant    |
| --------------- | ----------- | ------------ | :------------ |
| AL.GS.2014      | Sample      | Alabama      | Grain_sorghum |
| Col.GS.2015     | Sample      | Colony       | Grain_sorghum |
| Col.GS.NSS.2015 | Sample      | Colony       | Grain_sorghum |
| FL.GS.2014      | Sample      | Florida      | Grain_sorghum |
| FL.GS.2015      | Sample      | Florida      | Grain_sorghum |
| LA.GS.2014      | Sample      | Louisiana    | Grain_sorghum |
| TX.GS.2013      | Sample      | Texas        | Grain_sorghum |
| SA.GS.2013      | Sample      | South_Africa | Grain_sorghum |
| AL.SC.2014      | Sample      | Alabama      | Sugarcane     |
| FL.SC.2013      | Sample      | Florida      | Sugarcane     |
| FL.SC.2014      | Sample      | Florida      | Sugarcane     |
| LA.SC.2007.2009 | Sample      | Louisiana    | Sugarcane     |
| LA.SC.2015      | Sample      | Louisiana    | Sugarcane     |
| TX.SC.2015      | Sample      | Texas        | Sugarcane     |
| Control.1*      | Control     | Control      | Control       |
| Control.2**     | Control     | Control      | Control       |

*DNA extraction without aphids

**DNA extraction from solution after washing all aphids

#### Load packages

```R
library("phyloseq")
library("ggplot2")
library("ape")
library("picante")
library("reshape2")
library("ggpubr")
library("decontam")
library("DESeq2")
library("limma")
library("dplyr")
library("broom")
library("data.table")
library("scales")
library("ggrepel")
library("car")
library("lme4")
```

#### Import data

```R
data <- import_biom("otutable.biom", parseFunction = parse_taxonomy_greengenes)
metadata <- import_qiime_sample_data("metadata.txt")
tree <- read.tree("tree.tre")
tree <- root(tree, 1, resolve.root = T)
GM <- merge_phyloseq(data, metadata, tree)
```

#### Remove taxa associated as contaminants

```R
GM <- subset_taxa(GM, Class !="Chloroplast")
GM <- subset_taxa(GM, Genus !="Staphylococcus")
GM <- subset_taxa(GM, Genus !="Haemophilus")
GM <- subset_taxa(GM, Genus !="Oenothera")
GM <- subset_taxa(GM, Genus !="Gardnerella")
GM <- subset_taxa(GM, Genus !="Granulicatella")
GM <- subset_taxa(GM, Genus !="Leptotrichia")
GM <- subset_taxa(GM, Genus !="Prevotella")
GM <- subset_taxa(GM, Genus !="[Prevotella]")
GM <- subset_taxa(GM, Genus !="Ruminococcus")
GM <- subset_taxa(GM, Genus !="Streptococcus")
GM <- filter_taxa(GM, function (x) {sum(x > 0) > 1}, prune=TRUE) #remove singletons

sample_data(GM)$is.neg <- sample_data(GM)$Sample_type == "Control"
contamdf.prev <- isContaminant(GM, method="prevalence", neg="is.neg", threshold = 0.05)
table(contamdf.prev$contaminant)
cont.remove <- subset(contamdf.prev, contaminant == "TRUE")
cont.remove <- row.names(cont.remove)
allTaxa = taxa_names(GM)
allTaxa <- allTaxa[!(allTaxa %in% cont.remove)]
GM <-  prune_taxa(allTaxa, GM)
```

#### Normalize data and subtract from non-template controls

```R
diagdds = phyloseq_to_deseq2(GM, ~ Host_plant)
diagdds = estimateSizeFactors(diagdds)
diagdds = estimateDispersions(diagdds)
diagvst = getVarianceStabilizedData(diagdds)
diagdds.c <- removeBatchEffect(diagvst)

diagdds.c[,1] <- diagdds.c[,1]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,2] <- diagdds.c[,2]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,4] <- diagdds.c[,4]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,5] <- diagdds.c[,5]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,6] <- diagdds.c[,6]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,8] <- diagdds.c[,8]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,9] <- diagdds.c[,9]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,10] <- diagdds.c[,10]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,11] <- diagdds.c[,11]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,12] <- diagdds.c[,12]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,13] <- diagdds.c[,13]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,14] <- diagdds.c[,14]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[,15] <- diagdds.c[,15]-mean(diagdds.c[,3]-diagdds.c[,7])
diagdds.c[diagdds.c<0] <- 0

otu_table(GM) <- otu_table(diagdds.c, taxa_are_rows = TRUE)
```

#### Remove non-template control group

```R
GM <- subset_samples(GM, Sample_type!= "Control")
```

#### Remove singletons

```R
GM <- filter_taxa(GM, function (x) {sum(x > 0) > 1}, prune=TRUE)
```

#### PERMANOVA

Test for differences in the microbial community of aphids collected on different host plants, controlling for the sampling location.

```R
sampledf <- data.frame(sample_data(GM))
dist.mat <- phyloseq::distance(GM, method = "wunifrac")
perm <- how(nperm = 999)
set.seed(100)
pmv <- adonis2(dist.mat ~ Host_plant, strata = Location, data = sampledf, permutations = perm)
pmv
```

#### NMDS

Visualize how samples distribute according to the composition of their microbiota.

```R
nmds_ord <- ordinate(physeq = GM2, method = "NMDS", distance = "wunifrac", formula = ~ Host_plant)

#Build a graph using phyloseq::plot_ordination
```

#### Test for differences in bacterial taxa

This runs a mixed-effect linear model to test for differences in relative abundance within each bacterial genus between  different host plants.

```R
phy <- transform_sample_counts(GM, function(x) x/sum(x))
glom <- tax_glom(phy, taxrank = 'Genus')
dat <- data.table(psmelt(glom))
dat$Genus <- as.character(dat$Genus)
dat[, mean := mean(Abundance, na.rm = TRUE), by = "Genus"]
dat[(mean <= 0.01), Genus := "Others"]
dat <- dat[which(dat$Genus !="Others")]

#Build a graph using "dat" dataframe                               
                               
dat2 <- dat[, c(1:3,10,11,22)]

group.host <- group_by(dat2, Genus)
list.bact <-unique(c(as.character(group.host$Genus)))

model_calculator <- sapply(list.bact,  
                           function(x){
                             data.s <- group.host[which(group.host$Genus==x),]
                             model <- lmer(Abundance ~ Type + (1|Location), data = data.s)
                             aaa <-  Anova(model)
                             return(aaa)},
                           simplify = FALSE,USE.NAMES = TRUE)

res <- do.call(rbind, model_calculator)
res <- setDT(res, keep.rownames = TRUE)[]
colnames(res)[1] <- "Genus"

bbb <- data.table(dat2)
bbb <- bbb[,list(mean=mean(Abundance),sd=sd(Abundance)),by=c("Genus","Type")]
bbb <- dcast(data = bbb,formula = Genus~Type,fun.aggregate = sum, value.var = c("mean","sd"))

ccc <- merge(bbb, res, by="Genus")
```

#### References

Caporaso, J. G., C. L. Lauber, W. A. Walters, D. Berg-Lyons, J. Huntley, N. Fierer, S. M. Owens, J. Betley, L. Fraser, M. Bauer, N. Gormley, J. A. Gilbert, G. Smith, and R. Knight. 2012. Ultra-high-throughput microbial community analysis on the Illumina HiSeq and MiSeq platforms. The ISME Journal 6: 1621-1624.

Caporaso, J. G., J. Kuczynski, J. Stombaugh, K. Bittinger, F. D. Bushman, E. K. Costello, N. Fierer, A. G. Peña, J. K. Goodrich, J. I. Gordon, G. A. Huttley, S. T. Kelley, D. Knights, J. E. Koenig, R. E. Ley, C. A. Lozupone, D. McDonald, B. D. Muegge, M. Pirrung, J. Reeder, J. R. Sevinsky, P. J. Turnbaugh, W. A. Walters, J. Widmann, T. Yatsunenko, J. Zaneveld, and R. Knight. 2010. QIIME Allows Analysis of High-Throughput Community Sequencing Data. Nature Methods 7: 335-336.

Edgar, R. C. 2010. Search and clustering orders of magnitude faster than BLAST. Bioinformatics 26: 2460-2461.
