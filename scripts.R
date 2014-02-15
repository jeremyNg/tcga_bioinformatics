
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

rna.seq.tgf.up <- rna.seq[rownames(rna.seq)%in%tgf.up,]
rna.seq.tgf.down <- rna.seq[rownames(rna.seq)%in%tgf.down,]

rm(tgf.up,tgf.down,tgf.all)

tgf.dist <- dist(t(log2(rna.seq.tgf+1)),method = "maximum")
tgf.clust <- hclust(tgf.dist,method="ward")
clust.assignments <- cutree(tgf.clust,k=4) # this will break them up into 2 clusters, when the maximum distance is used with ward clustering algorithm
clust.col <- gsub("1","red",clust.assignments)
clust.col <- gsub("2","blue",clust.col)
clust.col <- gsub("3","orange",clust.col)
clust.col <- gsub("4","green",clust.col)


rna.seq.tgf.dif <- rna.seq.tgf[,which(clust.col%in%c("blue"))]
rna.seq.tgf.other <- rna.seq.tgf[,-which(clust.assignments==2)]

rna.seq.dif.med <- apply(rna.seq.tgf.dif,1,median)
rna.seq.others.med <- apply(rna.seq.tgf.other,1,median)
dif <- abs(rna.seq.dif.med-rna.seq.others.med)

rna.seq.1 <- rna.seq.tgf[,-which(clust.assignments==1)]
rna.seq.3 <- rna.seq.tgf[,-which(clust.assignments==3)]
rna.seq.4 <- rna.seq.tgf[,-which(clust.assignments==4)]
dif2 <- abs(rna.seq.dif.med-apply(rna.seq.1,1,median))
dif3 <- abs(rna.seq.dif.med-apply(rna.seq.3,1,median))
dif4 <- abs(rna.seq.dif.med-apply(rna.seq.4,1,median))
plot(dif,type="l",ylab="Difference (medians)",xlab="Gene",main="Difference in median counts",ylim=c(0,15000))
par(new=TRUE)
plot(dif2,type="l",ylab="Difference (medians)",xlab="Gene",main="Difference in median counts",ylim=c(0,15000),col="red")
par(new=TRUE)
plot(dif3,type="l",ylab="Difference (medians)",xlab="Gene",main="Difference in median counts",ylim=c(0,15000),col="blue")
par(new=TRUE)
plot(dif4,type="l",ylab="Difference (medians)",xlab="Gene",main="Difference in median counts",ylim=c(0,15000),col="green")
require(gplots)

heatmap.2(as.matrix(log10(rna.seq.tgf+1)),ColSideColors=clust.col,trace = "none",hclustfun = function(x) hclust(x,method = "ward"),distfun = function(x) dist(x,method = "maximum")) # plot heatmap for all the genes of interest in the TGFb response
coef.variation <- apply(rna.seq.tgf,1,function(x) sd(x)/mean(x)) #calculate a coefficient of variation for each gene
rna.seq.tgf.high.coef <- rna.seq.tgf[rownames(rna.seq.tgf)%in%rownames(rna.seq.tgf)[which(coef.variation>mean(coef.variation))],]
heatmap.2(as.matrix(log10(rna.seq.tgf.high.coef+1)),ColSideColors = clust.col,trace = "none",hclustfun = function(x) hclust(x,method = "ward"),distfun = function(x) dist(x,method = "maximum"),col=topo.colors(100)) # plot heatmap for all the genes of interest in the TGFb response

heatmap.2(as.matrix(log10(rna.seq.tgf.dif+1)),trace = "none",hclustfun = function(x) hclust(x,method = "ward"),distfun = function(x) dist(x,method = "maximum"),col=topo.colors(100),ColSideColors = clust.cols) # plot heatmap for all the genes of interest in the TGFb response
rm(rna.seq.z.std)
rm(rna.seq.tgf.high.coef,coef.variation)
rm(clust.assignments,tgf.clust,tgf.dist,genes)


# to calculate NXN matrix
compute.cor <- function(x){
    results <- array(dim=c(ncol(x),ncol(x)))
    colnames(results) <- colnames(x) # names by samples
    rownames(results) <- colnames(results)
    for (i in 1:ncol(x)){
        for (j in 1:ncol(x)){
            results[i,j] <- results[j,i] <- cor(x[,i],x[,j],method="pearson")
            }
        }
    return(results)
    }
cor.matrix <- compute.cor(rna.seq.tgf)
cor.matrix.z <- apply(cor.matrix,1,function(x){z <- (x-mean(x))/sd(x)}) # to perform z-normalization of the correlation coefficients

# are the correlation coefficients the same across all samples?
require(reshape)
require(ggplot2)

cor.melt <- melt(cor.matrix)
ggplot(cor.melt,aes(as.factor(X1),X2,group=X1))+geom_tile(aes(fill=value))+scale_fill_gradient(high="red",low="yellow")+ggtitle("Correlation of expression between samples")+xlab("Sample")+ylab("Sample")

rm(cor.melt,cor.matrix,compute.cor)

