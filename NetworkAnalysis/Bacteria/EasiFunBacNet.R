# R script for Cross-domain network analyis
## Example: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0393-0#Sec10
## Package: https://github.com/zdk123/SpiecEasi

## Loading library
library(phyloseq)
library(SpiecEasi)
library(igraph)

## Import Bacterial data
## Bacterial ASV(OTU) table
Bac.otus <- read.table("MC2018.network.BAC.otu_table.csv",
                   header=T,sep=",",row.names = 1)
Bac.otumat <- as(as.matrix(Bac.otus), "matrix")
Bac.OTU = otu_table(Bac.otumat, taxa_are_rows = TRUE)

## Bacterial Taxanomy table
Bac.taxmat <- read.table("BAC.network.taxonomy.csv", header=T,sep=",",row.names=1)
Bac.taxmat <- as(as.matrix(Bac.taxmat),"matrix")
Bac.TAX = tax_table(Bac.taxmat)

## Bacterial Metadata
meta = read.table("Mojave_mappingfile_8-Aug-2018.txt",
                  header=TRUE,row.names=1,
                  sep="\t",stringsAsFactors=FALSE)
sampleData <- sample_data(meta)

## Construct Bacterail Phyloseq Object
Bac.physeq = phyloseq(Bac.OTU,Bac.TAX,sampleData)

## Prune ASVs (discard any ASV present less than 5 times)
Bac.physeq.prune = prune_taxa(taxa_sums(Bac.physeq) > 5, Bac.physeq)

## Rarefy data
set.seed(711)
Bac.physeq.prune.rarefy = rarefy_even_depth(Bac.physeq.prune, sample.size = 37435, replace = FALSE, trimOTUs = TRUE)
Bac.physeq.prune.rarefy

## Remove other kingdoms (keep bacteria)
Bac.physeq.prune.rarefy = subset_taxa(Bac.physeq.prune.rarefy, Kingdom != "Rhizaria")
Bac.physeq.prune.rarefy = subset_taxa(Bac.physeq.prune.rarefy, Kingdom != "Chromista")
Bac.physeq.prune.rarefy = subset_taxa(Bac.physeq.prune.rarefy, Kingdom != "Unassigned")
Bac.physeq.prune.rarefy = subset_taxa(Bac.physeq.prune.rarefy, Kingdom != "D_0__Eukaryota")
Bac.physeq.prune.rarefy

## Bacterial phyloseq object for surface subset
Bac.physeq.prune.rarefy.SF = subset_samples(Bac.physeq.prune.rarefy, Layer=="Surface")

## Bacterial phyloseq object for subsurface subset
Bac.physeq.prune.rarefy.SUB = subset_samples(Bac.physeq.prune.rarefy, Layer=="Subsurface")

## Import Fungal data
## Fungal OTU table
Fun.otus <- read.table("MC2018.network.FG.otu_table.csv",header=T,sep=",",row.names=1)
Fun.otumat <- as(as.matrix(Fun.otus), "matrix")
Fun.OTU = otu_table(Fun.otumat, taxa_are_rows = TRUE)

## Fungal Taxonomy table
Fun.taxmat <- read.table("FG.network.taxonomy.csv", header=T,sep=",",row.names=1)
Fun.taxmat <- as(as.matrix(Fun.taxmat),"matrix")
Fun.TAX = tax_table(Fun.taxmat)

## Construct fungal phyloseq object 
Fun.physeq = phyloseq(Fun.OTU,Fun.TAX,sampleData)

## Remove singletons
Fun.physeq.prune = prune_taxa(taxa_sums(Fun.physeq) > 1, Fun.physeq)

## Rarefy fungal data
set.seed(1)
Fun.physeq.prune.rarefy = rarefy_even_depth(Fun.physeq.prune, sample.size = 6842, replace = FALSE, trimOTUs = FALSE)
Fun.physeq.prune.rarefy

## Remove other kingdoms (keep fungi)
Fun.physeq.prune.rarefy = subset_taxa(Fun.physeq.prune.rarefy, Kingdom != "Rhizaria")
Fun.physeq.prune.rarefy = subset_taxa(Fun.physeq.prune.rarefy, Kingdom != "Chromista")
Fun.physeq.prune.rarefy = subset_taxa(Fun.physeq.prune.rarefy, Kingdom != "Unassigned")
Fun.physeq.prune.rarefy = subset_taxa(Fun.physeq.prune.rarefy, Kingdom != "D_0__Eukaryota")
Fun.physeq.prune.rarefy = subset_taxa(Fun.physeq.prune.rarefy, Kingdom != "NA")
Fun.physeq.prune.rarefy

## Fungal phyloseq object for surface subset
Fun.physeq.prune.rarefy.SF = subset_samples(Fun.physeq.prune.rarefy, Layer=="Surface")

## Fungal phyloseq object for subsurface subset
Fun.physeq.prune.rarefy.SUB = subset_samples(Fun.physeq.prune.rarefy, Layer=="Subsurface")

## Bacterial phyloseq object for surface subset (top 100 ASV)
Bac.physeq.prune.rarefy.SF.top100 = prune_taxa(names(sort(taxa_sums(Bac.physeq.prune.rarefy.SF), TRUE))[1:38], Bac.physeq.prune.rarefy.SF)

## Fungal phyloseq object for surface subset (top 100 ASV)
Fun.physeq.prune.rarefy.SF.top100 = prune_taxa(names(sort(taxa_sums(Fun.physeq.prune.rarefy.SF), TRUE))[1:38], Fun.physeq.prune.rarefy.SF)

## Network analysis using SpiecEasi
pargs <- list(rep.num = 50, seed = 10010, ncores = 6)
FG.Bac.network <- spiec.easi(list(Bac.physeq.prune.rarefy.SF.top100, Fun.physeq.prune.rarefy.SF.top100), method='glasso', nlambda=40,
              lambda.min.ratio=1e-2, pulsar.params = pargs)

## Merged Bacterial and Fungal Phyloseq objects for combined taxnomic data
FG.BAC.top100.merged = merge_phyloseq(Bac.physeq.prune.rarefy.SF.top100, Fun.physeq.prune.rarefy.SF.top100)
FG.BAC.top100.merged

## Create igraph object from SpiecEasi network
easi.fun.bac.net <- adj2igraph(getRefit(FG.Bac.network), rmEmptyNodes=TRUE, vertex.attr=list(name=taxa_names(FG.BAC.top100.merged)))

## Extract dataframe from igraph
net.df = as_long_data_frame(easi.fun.bac.net)
cdf = data.frame(OTU1 = net.df$`ver[el[, 1], ]`, OTU2 = net.df$`ver2[el[, 2], ]`, WEIGHT = as.vector(net.df$weight), stringsAsFactors = TRUE )
cdf

## Save dataframe in csv for future use in circlize graph
write.csv(cdf, file="EasiFunBacNet.csv")

