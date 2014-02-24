# READ in csv file with RNA seq data from L93 onwards! can skip the first few parts
# CODE CHUNK 1: read in files as objects

gen.manifest <- function(){
manifest1 <- read.table("batch1_manifest.txt",sep="\t",header=T)
manifest2 <- read.table("batch2_manifest.txt",sep="\t",header=T)
manifest3 <- read.table("batch3_manifest.txt",sep="\t",header=T)
all.manifest <- rbind(manifest1,manifest2,manifest3)
return(all.manifest)
}
full.manifest <- gen.manifest()

rm(gen.manifest)

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

rm(c.file,c.in,c.outs,all.in.present,samples.id)
# CODE CHUNK 2: RNA seq data analysis

##########################################################################
# Pseudo codes here
# 1. combine all the RNAseq data into 1 single frame
# 2. identify the genes to perform clustering along (single the common genes will dominate an all vs all clustering)
# 3. perform clustering to identify signatures of tfg beta aberration
# 4. identify samples which we think reasonably represent those with tgfb abberation
# 5. use those to construct a molecular signature profile
##########################################################################
rna.seq.manifest <- full.manifest[full.manifest$Platform.Type=="RNASeqV2",]
rna.seq.fp <- paste0(getwd(),"/",rna.seq.manifest$File.Name)
# all the data will be stored in a single list
read.rnaseq <- function(file.paths){
    f <- read.table(file=file.paths,header=T)
    all.dat <- f$normalized_count
    return(all.dat)
    }
rna.seq <- lapply(rna.seq.fp,read.rnaseq)
# lapply(rna.seq,function(x){return(length(unlist(x)))})
names(rna.seq) <- rna.seq.manifest$Sample
rna.seq <- as.data.frame(rna.seq)
get.gn <- function(){
    f <- read.table(rna.seq.fp[1],header = T)
    gn <- f$gene_id
    return(gn)
    }
g.n <- get.gn() # there is a duplicate in SLC35E2; remove the second entry
g.n <- g.n[-16273]
rna.seq <- rna.seq[-16273,]
rownames(rna.seq) <- g.n
write.csv(rna.seq,"combinedRNAseq.csv")

rm(g.n,get.gn,rna.seq.fp,rna.seq.manifest,read.rnaseq)

# read in genesets
tgf.up <- read.table("~/Downloads/tgf_up.txt",sep="\t")
tgf.up <- tgf.up[-(1:2),]
tgf.up <- as.vector(tgf.up)
tgf.down <- read.table("~/Downloads/tgf_down.txt",sep="\t")
tgf.down <- tgf.down[-(1:2),]
tgf.down <- as.vector(tgf.down)
tgf.all <- unique(c(tgf.up,tgf.down)) # set of all genes identified in tgf beta response in experimental studies

rna.seq.tgf <- rna.seq[rownames(rna.seq)%in%tgf.all,] # expression data for all genes involved in tgfbeta

tgf.scale <- scale(rna.seq.tgf) # scales data using Z transform
#tgf.cov <- apply(tgf.scale,1,function(x){y <-sd(x)/mean(x)})
#tgf.mean <- apply(tgf.scale,1,mean)

#rna.seq.tgf.up <- rna.seq[rownames(rna.seq)%in%tgf.up,]
#rna.seq.tgf.down <- rna.seq[rownames(rna.seq)%in%tgf.down,]

#rm(tgf.up,tgf.down,tgf.all)

# code chunk for hclustering

# Perform bootstrapping o test structure
#require(ClassDiscovery)
#boot <- BootstrapClusterTest(log10(rna.seq.tgf+1),cutHclust,k=4,method="ward",metric="euclid",nTimes = 1000) # uses 1000 iterations to perform bootstrapping of SAMPLE clustering
#image(boot,col=blueyellow(64)) # plots the counts
#require(gplots)

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
tgf.dist <- dist(t(rna.seq.tgf),method = "euclidean")
tgf.clust <- hclust(tgf.dist,method="ward")
clust.assignments <- cutree(tgf.clust,k=4) # this will break them up into 2 clusters, when the maximum distance is used with ward clustering algorithm
clust.col <- gsub("1","red",clust.assignments)
clust.col <- gsub("2","blue",clust.col)
clust.col <- gsub("3","orange",clust.col)
clust.col <- gsub("4","green",clust.col)

# plot the heatmap
#heatmap.2(as.matrix(log10(rna.seq.tgf+1)),trace = "none",hclustfun = function(x) hclust(x,method = "ward"),distfun = function(x) dist(x,method = "euclidean"),col=topo.colors(100),ColSideColors = clust.col) # plot heatmap for all the genes of interest in the TGFb response

# code for calculating similarity index
rna.seq.tgf1 <- rna.seq.tgf[,which(clust.col=="blue")]
rna.seq.tgf2 <- rna.seq.tgf[,-c(which(clust.col=="blue"))]
rna.seq.med1 <- apply(log10(rna.seq.tgf1+1),1,median) # summarize median for group 1
rna.seq.med2 <- apply(log10(rna.seq.tgf2+1),1,median) # summarize median for other groups
# ID the genes that are up regulated by TGFb
ups <- which(rownames(rna.seq.tgf)%in%tgf.up)
# ...and those down regulated by TGFb
downs <- which(rownames(rna.seq.tgf)%in%tgf.down)
# break them up, we will now have 4 lists
#ups.1 <- rna.seq.med1[ups] # upregulated, grp 1
#ups.2 <- rna.seq.med2[ups] # upregulated, grp 2
#downs.1 <- rna.seq.med1[downs] # downregulated, grp 1
#downs.2 <- rna.seq.med2[downs] # downregulated, grp 2
# using a simple method to determine
#score.grp1 <- length(which(ups.1>ups.2))+length(which(downs.1<downs.2))
#score.grp2 <- length(which(ups.2>ups.1))+length(which(downs.2<downs.1))

# to perform a boot strap test for the score ->
# Question: Are the similarity scores really significantly different
# logic -> choose n samples for each grp, then after that, select n number of points, then calculate the statistical significance with a paired t-test

bootstrapscores <- function(x,y,n.iter=10000){
    # initialize the counters
    #x is grp 1, y is grp 2
    j <- 1
    grp1.list <- grp2.list <- rep(NA,n.iter)
    while (j<=n.iter){
        # randomly select some genes
        g.l <- sample(x=1:1000,size=372,replace=TRUE)
        cat("iteration",j,"of", n.iter,"\n")
        # to select the samples ##
        # to build up the scores
        sm.a <- x[g.l,]
        #print(dim(sm.a))
        sm.b <- y[g.l,]
        #print(dim(sm.b))
        # to calculate the medians
        m.a <- apply(log10(sm.a+1),1,median)
        m.b <- apply(log10(sm.b+1),1,median)
        ups <- which(rownames(rna.seq.tgf)%in%tgf.up)
# ...and those down regulated by TGFb
        downs <- which(rownames(rna.seq.tgf)%in%tgf.down)
        ups.1 <- m.a[ups] # upregulated, grp 1
        ups.2 <- m.b[ups] # upregulated, grp 2
        downs.1 <- m.a[downs] # downregulated, grp 1
        downs.2 <- m.b[downs] # downregulated, grp 2
        score.grp1 <- length(which(ups.1>ups.2))+length(which(downs.1<downs.2))
        score.grp2 <- length(which(ups.2>ups.1))+length(which(downs.2<downs.1))
        grp1.list[j] <- score.grp1
        grp2.list[j] <- score.grp2
        n.genes[j] <- n.
        #cat(score.grp1,score.grp2,"\n")
        j <- j+1 # increases the counter
        }
    res <- data.frame(Grp1=grp1.list,Grp2=grp2.list)
    return(res)
    }
bootstraprun <- bootstrapscores(rna.seq.tgf1,rna.seq.tgf2,10000) # does 10000 iterations
