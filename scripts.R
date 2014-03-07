# CODE CHUNK 1: Making the required manifest files to read files etc;

gen.manifest <- function(){
manifest1 <- read.table("batch1_manifest.txt",sep="\t",header=T)
manifest2 <- read.table("batch2_manifest.txt",sep="\t",header=T)
manifest3 <- read.table("batch3_manifest.txt",sep="\t",header=T)
all.manifest <- rbind(manifest1,manifest2,manifest3)
return(all.manifest)
} # function to generate the manifest file first
full.manifest <- gen.manifest() # makes the manifest file

rm(gen.manifest) # removes the function to make the manifest file

full.manifest <- full.manifest[full.manifest$Platform.Type!="METADATA"||full.manifest$Platform.Type!="RNASeq",] # removes all the metadata associated in the manifest files, as well as RNAseq data
all.in.present <- list.files() # this is a list of all the files in the data directory, in which other redundant data type has been deleted
full.manifest <- full.manifest[full.manifest$File.Name%in%all.in.present,] # reduces
full.manifest <- full.manifest[,-c(2:4)] # removes redundant information for the samples
full.manifest <- full.manifest[-c(grep("nocnv_hg19.seg.txt|rsem.genes.results",full.manifest$File.Name)),] # gets rid of all the nocnv_hg19 files and rsem results prior to normalization
samples.id <- as.vector(unique(full.manifest$Sample)) # gets a list of all the samples that are represented in the manifest file
# to ensure that the samples have a paired CNV data and RNAseq data- there should be 2 files for each, we can lapply over the list of all samples
c.file <- function(x){
    return(nrow(full.manifest[full.manifest$Sample==x,]))
    }
c.in <- unlist(lapply(samples.id,c.file))
c.outs <- samples.id[which(c.in!=2)]
full.manifest <- full.manifest[!full.manifest$Sample%in%c.outs,] # gets rid of all samples which are not paired, or have weird annotations

rm(c.file,c.in,c.outs,all.in.present,samples.id) # clear up not needed objects on the workspace

# CODE CHUNK 2: RNA seq data analysis

##########################################################################
# Pseudo codes here
# 1. combine all the RNAseq data into 1 single frame
# 2. identify the genes to perform clustering along (single the common genes will dominate an all vs all clustering)
# 3. perform clustering to identify signatures of tfg beta aberration
# 4. identify samples which we think reasonably represent those with tgfb abberation
# 5. use those to construct a molecular signature profile
##########################################################################

rna.seq.manifest <- full.manifest[full.manifest$Platform.Type=="RNASeqV2",] # gets only RNAseqV2 objects
rna.seq.fp <- paste0(getwd(),"/",rna.seq.manifest$File.Name) # generates the file paths
read.rnaseq <- function(file.paths){
    f <- read.table(file=file.paths,header=T)
    all.dat <- f$normalized_count
    return(all.dat)
    } # function to read in a RNASeq file from a path as a list
rna.seq <- lapply(rna.seq.fp,read.rnaseq) # lapply over the paths
names(rna.seq) <- rna.seq.manifest$Sample
rna.seq <- as.data.frame(rna.seq) # converts to a dataframe
get.gn <- function(){
    f <- read.table(rna.seq.fp[1],header = T)
    gn <- f$gene_id
    return(gn)
    } # function to get the gene names
g.n <- get.gn() # there is a duplicate in SLC35E2; remove the second entry
g.n <- g.n[-16273]
rna.seq <- rna.seq[-16273,]
rownames(rna.seq) <- g.n # assign gene name as rownames of the dataframe
write.csv(rna.seq,"combinedRNAseq.csv") # exports as a csv file

rm(g.n,get.gn,rna.seq.fp,rna.seq.manifest,read.rnaseq) # removes the unwanted workspace objects

# CODE CHUNK 3: Using gene set information for clustering

# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' #
# Comparative method to get clsutering
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' #
# there is no need for us to have an absolutely healthy control; actually since we are interested in identifying the most responsive systems, we can just simply find the samples which show the largest change in the TGFB system
# rescale samples using z transform so that the mean of each sample is 0 and SD is 1
# calculate the mean for each and every gene, then after that, calculate a logFC

####### Code chunk below starts from a csv file that was exported earlier for the entire gene expression array
rna.seq <- read.csv("combinedRNAseq.csv")# calls in the saved RNAseq data from earlier
rownames(rna.seq) <- rna.seq[,1]
rna.seq <- rna.seq[,-1]
tgf.up <- read.table("~/Downloads/tgf_up.txt",sep="\t")
tgf.up <- tgf.up[-(1:2),]
tgf.up <- as.vector(tgf.up)
tgf.down <- read.table("~/Downloads/tgf_down.txt",sep="\t")
tgf.down <- tgf.down[-(1:2),]
tgf.down <- as.vector(tgf.down)
tgf.all <- unique(c(tgf.up,tgf.down)) # set of all genes identified in tgf beta resp

rna.seq.ztransformed <- apply(log2(rna.seq+1),2,scale) # scales along the columns by sample
rownames(rna.seq.ztransformed) <- rownames(rna.seq)
colnames(rna.seq.ztransformed) <- colnames(rna.seq)
rna.seq.z.tgf <- rna.seq.ztransformed[rownames(rna.seq.ztransformed)%in%tgf.all,] # gets all the TGFb response genes
#rna.seq.z.tgf.in <- rna.seq.z.tgf[which(apply(rna.seq.tgf.fc,1,sd)> 3.253159),]
rna.seq.tgf.means <- apply(rna.seq.z.tgf,1,mean) #finds the mean for each gene
rna.seq.tgf.fc <- apply(rna.seq.z.tgf,2,function(x) x/rna.seq.tgf.means) # computes the log FC of the data
# we only want to use the first 50 most variable response genes because they are most likely to be indicative of the subtypes; maximize the SD
rna.seq.tgf.variable <- rna.seq.z.tgf[which(apply(rna.seq.tgf.fc,1,sd)> 3.253159),] #20 most variable genes
tgf.dist.variable <- dist(t(rna.seq.tgf.variable),method = "euclidean") # calculates a distance matrix sing euclidean distance
tgf.clust <- hclust(tgf.dist.variable,method="ward") # does clustering with the wards agglomerative method
clust.assignments <- cutree(tgf.clust,k=2) # this will break them up into 2 clusters, when the maximum distance is used with ward clustering algorithm
clusterassignments <- gsub(1,"blue",clust.assignments)
clusterassignments <- gsub(2,"red",clusterassignments)

require(gplots)
heatmap.2(rna.seq.tgf.variable,trace="none",hclustfun = function(x) hclust(x,method = "ward"),distfun = function(x) dist(x,method = "euclidean"),ColSideColors = clusterassignments) # plots a heatmap with the color bars on the top indicating the membership
# NOT useful #
# not run #
#require(FactoMineR)
#pcaresults <- PCA(t(rna.seq.tgf.variable),graph=TRUE)
#pcaplotting <- pcaresults$ind$coord
#write.table(pcaplotting,"pca.csv",sep="\t")# performs PCA and stores for PC 1 to PC 5 only to be exported; then after some mild parsing, use Adelene's Plot_PCA codes to plot the PCA contour map (saved as PCAplot.png)
# NOT useful #


# CODE CHUNK 4: constructing a representative profile

# Using medians of z-transformed counts
rna.seq.blue <- rna.seq.ztransformed[,which(clusterassignments=="blue")] #gets all the genes for samples in the blue cluster
rna.seq.z.profile <- apply(rna.seq.blue,1,median) # this is the profile based on the median


