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
rna.seq.tgf <- rna.seq[rownames(rna.seq)%in%tgf.all,] # expression data for all gene
tgf.dist <- dist(t(rna.seq.tgf),method = "euclidean") # calculates a distance matrix sing euclidean distance
tgf.clust <- hclust(tgf.dist,method="ward") # does clustering with the wards agglomerative method
clust.assignments <- cutree(tgf.clust,k=4) # this will break them up into 4 clusters, when the maximum distance is used with ward clustering algorithm
# simply to generate a column color bar on the top according to their cluster assignments
clust.col <- gsub("1","red",clust.assignments)
clust.col <- gsub("2","blue",clust.col)
clust.col <- gsub("3","orange",clust.col)
clust.col <- gsub("4","green",clust.col)

require(gplots)
heatmap.2(as.matrix(log10(rna.seq.tgf+1)),trace = "none",hclustfun = function(x) hclust(x,method = "ward"),distfun = function(x) dist(x,method = "euclidean"),col=topo.colors(100),ColSideColors = clust.col) # plot heatmap for all the genes of interest in the TGFb response

# Perform bootstrapping to test structure and robustness of clustering
require(ClassDiscovery)
boot <- BootstrapClusterTest(log10(rna.seq.tgf+1),cutHclust,k=4,method="ward",metric="euclid",nTimes = 1000,verbose = FALSE) # uses 1000 iterations to perform bootstrapping of SAMPLE clustering
image(boot,col=blueyellow(64),dendrogram=tgf.clust,ColSideColors=clust.col,RowSideColors=clust.col) # plots the counts


grp1 <- rna.seq.tgf[,which(clust.col=="blue")] # blue cluster
grp2 <- rna.seq.tgf[,which(clust.col!="blue")] # other clusters

# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' #
# the conflating of the other group might be a problem if we cannot show that the blue cluster is the better cluster
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' #

bootstrapscores <- function(x,y,n.iters=1000){
    # function to calculate score of suitability of grp 1 vs grp 2
    # to guarantee reproducibility of results only #
    set.seed(0)
    seeds <- sample(1:100000,n.iters,replace=TRUE) # generates a list of random seeds
    result.list <- rep(NA,n.iters)
    for (i in 1:n.iters){
    set.seed(seeds[i])
     cat("iteration",i," with seed",seeds[i],"\n")
    gl <- sample(1:372,372,replace=TRUE) # gets a list of genes by their rows
    grp1.list <- x[gl,]
    grp2.list <- y[gl,]
    grp1.up <- x[rownames(grp1.list)%in%tgf.up,]
    grp2.up <- y[rownames(grp2.list)%in%tgf.up,]
    grp1.down <- x[rownames(grp1.list)%in%tgf.down,]
    grp2.down <- y[rownames(grp2.list)%in%tgf.down,]
    up.r <- rep(NA,nrow(grp1.up))
    down.r <- rep(NA,nrow(grp1.down))
    for (j in 1:nrow(grp1.up)){
        z <- wilcox.test(as.numeric(grp1.up[j,]),as.numeric(grp2.up[j,]),alternative = "greater")
        up.r[j] <- z$p.value
    }
    for (k in 1:nrow(grp1.down)){
        z <- wilcox.test(as.numeric(grp1.down[k,]),as.numeric(grp2.down[k,]),alternative = "less")
        down.r[k] <- z$p.value
    }
    b <- sum(up.r,down.r)
    result.list[i] <- b
    } # closes the for loop over the iterations
   return(result.list/372)
}
bluescore <- bootstrapscores(grp1,grp2,1000) # bootstrap for the blue cluster
# both runs will have the same set of seeds, which means both will have the same set of genes compared throughout the simulations - fair comparison but not sure how robust though #
miscscore <- bootstrapscores(grp2,grp1,1000) # bootstrap for other clusters

# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' #
# high p-value = high likelihood of similarity between grp 1 and 2, mean less significant difference between the expression of the 2 groups -
#we want to maximize the significance of difference in the a priori expected direction, hence want to minimize the sum of all the p-values
# for instuitive interpretation take the complement of the p-value (1-p)
# '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' #

bluecom <- 1-bluescore # to calculate the reverse for interpretability
miscom <- 1-miscscore
boxplot(data.frame(Blue=bluecom,Others=miscom)) # boxplot to visualize scores







