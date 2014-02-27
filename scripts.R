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
boot <- BootstrapClusterTest(log10(rna.seq.tgf+1),cutHclust,k=4,method="ward",metric="euclid",nTimes = 1000) # uses 1000 iterations to perform bootstrapping of SAMPLE clustering
image(boot,col=blueyellow(64)) # plots the counts

# to perform a boot strap test for the score::
# score defined as-> if in grp 1 is higher/lower relative to their expected expression values then score as +1 else 0
# Question: Are the similarity scores really significantly different
# logic -> choose  n number of points, then calculate the statistical significance with a paired t-test; robustness of the scores

# code for calculating similarity index
rna.seq.tgf1 <- rna.seq.tgf[,which(clust.col=="blue")]
rna.seq.tgf2 <- rna.seq.tgf[,-c(which(clust.col=="blue"))]

bootstrapscores <- function(x,y,n.iter=10000){
    # initialize the counters
    #x is grp 1, y is grp 2
    j <- 1
    grp1.list <- grp2.list <- rep(NA,n.iter)
    while (j<=n.iter){
        # randomly select some genes
        g.l <- sample(x=1:372,size=372,replace=TRUE)
        cat("iteration",j,"of", n.iter,"\n")
        sm.a <- x[g.l,] # grp 1 gene expression matrix with draw
        sm.b <- y[g.l,] # grp 2 gene expression matrix with draw
        m.a <- apply(log10(sm.a+1),1,median) # median group 1
        m.b <- apply(log10(sm.b+1),1,median) # median group 2
        ups <- which(rownames(rna.seq.tgf)%in%tgf.up)
        downs <- which(rownames(rna.seq.tgf)%in%tgf.down)
        ups.1 <- m.a[ups] # upregulated, grp 1
        ups.2 <- m.b[ups] # upregulated, grp 2
        downs.1 <- m.a[downs] # downregulated, grp 1
        downs.2 <- m.b[downs] # downregulated, grp 2
        score.grp1 <- length(which(ups.1>ups.2))+length(which(downs.1<downs.2))
        score.grp2 <- length(which(ups.2>ups.1))+length(which(downs.2<downs.1))
        grp1.list[j] <- score.grp1
        grp2.list[j] <- score.grp2
        #cat(score.grp1,score.grp2,"\n")
        j <- j+1 # increases the counter j
        }
    res <- data.frame(Grp1=grp1.list,Grp2=grp2.list) # makes a data frame
    return(res)
    }
bootstraprun <- bootstrapscores(rna.seq.tgf1,rna.seq.tgf2,10000) # does 10000 iterations
