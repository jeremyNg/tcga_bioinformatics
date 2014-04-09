# more important workspace objects
1. Manifest files -> full.manifest
2. List of samples in the cluster->All clusters in red and blue in clusterassignments
3. Correlations of the breast and non breast cell lines
4. Correlations of CNV

require(multicore) # for parallelizations
require(gplots)
require(reshape2)
require(ggplot2)
require(CNTools)

gen.manifest <- function(){
manifest1 <- read.table("batch1_manifest.txt",sep="\t",header=T)
manifest2 <- read.table("batch2_manifest.txt",sep="\t",header=T)
manifest3 <- read.table("batch3_manifest.txt",sep="\t",header=T)
all.manifest <- rbind(manifest1,manifest2,manifest3)
return(all.manifest)
} # function to generate the manifest file first
full.manifest <- gen.manifest() # makes the manifest file

rm(gen.manifest) # removes the function to make the manifest file

# removing all the metadata associated in the manifest files, as well as RNAseq data
full.manifest <- full.manifest[full.manifest$Platform.Type!="METADATA"||full.manifest$Platform.Type!="RNASeq",]
# this is a list of all the files in the data directory, in which other redundant data type has been deleted
all.in.present <- list.files()
full.manifest <- full.manifest[full.manifest$File.Name%in%all.in.present,]
full.manifest <- full.manifest[,-c(2:4)] # removes redundant information for the samples
full.manifest <- full.manifest[-c(grep("nocnv_hg19.seg.txt|rsem.genes.results",full.manifest$File.Name)),] # gets rid of all the nocnv_hg19 files and rsem results prior to normalization
samples.id <- as.vector(unique(full.manifest$Sample)) # gets a list of all the samples that are represented in the manifest file
# to ensure that the samples have a paired CNV data and RNAseq data- there should be 2 files for each, we can lapply over the list of all samples
c.file <- function(x){
    return(nrow(full.manifest[full.manifest$Sample==x,]))
    }
c.in <- unlist(lapply(samples.id,c.file))
c.outs <- samples.id[which(c.in!=2)]
# gets rid of all samples which are not paired, or have weird annotations
full.manifest <- full.manifest[!full.manifest$Sample%in%c.outs,]
# clear up unwanted objects in the workspace
rm(c.file,c.in,c.outs,all.in.present,samples.id)


# CODE CHUNK 2: RNA seq data analysis


# gets only RNAseqV2 objects
rna.seq.manifest <- full.manifest[full.manifest$Platform.Type=="RNASeqV2",]
# generates the file paths
rna.seq.fp <- paste0(getwd(),"/",rna.seq.manifest$File.Name)
# function to read in a RNASeq file from a path as a list
read.rnaseq <- function(file.paths){
    f <- read.table(file=file.paths,header=T)
    all.dat <- f$normalized_count
    return(all.dat)
    }
# lapply over the paths with 4 CPU
rna.seq <- mclapply(rna.seq.fp,read.rnaseq,mc.cores = 4)
names(rna.seq) <- rna.seq.manifest$Sample
rna.seq <- as.data.frame(rna.seq)

# function to get the gene names
get.gn <- function(){
    f <- read.table(rna.seq.fp[1],header = T)
    gn <- f$gene_id
    return(gn)
    }
g.n <- get.gn() # there is a duplicate in SLC35E2; remove the second entry
g.n <- g.n[-16273]
rna.seq <- rna.seq[-16273,]
# assign gene name as rownames of the dataframe
rownames(rna.seq) <- g.n
#exporting RNASeq as a CSV file
write.csv(rna.seq,"combinedRNAseq.csv")
# removes the unwanted workspace objects
rm(g.n,get.gn,rna.seq.fp,rna.seq.manifest,read.rnaseq)

# calls in the saved RNAseq data from earlier in a new session
rna.seq <- read.csv("combinedRNAseq.csv")
rownames(rna.seq) <- rna.seq[,1]
rna.seq <- rna.seq[,-1]
tgf.up <- read.table("~/Downloads/tgf_up.txt",sep="\t")
tgf.up <- tgf.up[-(1:2),]
tgf.up <- as.vector(tgf.up)
tgf.down <- read.table("~/Downloads/tgf_down.txt",sep="\t")
tgf.down <- tgf.down[-(1:2),]
tgf.down <- as.vector(tgf.down)
# set of all genes identified in tgf beta response
tgf.all <- unique(c(tgf.up,tgf.down))
# clears workspace of unwanted objects
rm(tgf.up,tgf.down)

# scales along the columns by sample
rna.seq.ztransformed <- apply(log2(rna.seq+1),2,scale)
rownames(rna.seq.ztransformed) <- rownames(rna.seq)
colnames(rna.seq.ztransformed) <- colnames(rna.seq)

# getting all the TGFb response genes
rna.seq.z.tgf <- rna.seq.ztransformed[rownames(rna.seq.ztransformed)%in%tgf.all,]
rna.seq.tgf.means <- apply(rna.seq.z.tgf,1,mean) #finds the mean for each gene
rna.seq.tgf.fc <- apply(rna.seq.z.tgf,2,function(x) x/rna.seq.tgf.means) # computes the log FC of the data
# we only want to use the first 20 most variable response genes because they are most likely to be indicative of the subtypes

rna.seq.tgf.variable <- rna.seq.z.tgf[which(apply(rna.seq.tgf.fc,1,sd)> 3.253159),] #20 most variable genes

rm(rna.seq.tgf.fc,rna.seq.tgf.means)
tgf.dist.variable <- dist(t(rna.seq.tgf.variable),method = "euclidean") # calculates a distance matrix using euclidean distance
tgf.clust <- hclust(tgf.dist.variable,method="ward") # does clustering with the wards agglomerative method
clust.assignments <- cutree(tgf.clust,k=2) # this will break them up into 2 clusters

clusterassignments <- gsub(1,"blue",clust.assignments)
clusterassignments <- gsub(2,"red",clusterassignments)


heatmap.2(rna.seq.tgf.variable,trace="none",hclustfun = function(x) hclust(x,method = "ward"),distfun = function(x) dist(x,method = "euclidean"),ColSideColors = clusterassignments) # plots a heatmap with the color bars on the top indicating the membership


# Using medians of z-transformed counts
rna.seq.blue <- rna.seq.ztransformed[,which(clusterassignments=="blue")] #gets all the genes for samples in the blue cluster
rna.seq.z.profile <- apply(rna.seq.blue,1,median) # this is the profile based on the median
rna.seq.z.profile <- data.frame(Name=rownames(rna.seq),Reads=rna.seq.z.profile)
rownames(rna.seq.z.profile) <- rna.seq.z.profile$Name #This is the reference profile

#clearing of workspace objects
rm(clust.assignments,rna.seq.tgf.fc,rna.seq.tgf.means,rna.seq.tgf.variable)
rm(rna.seq.blue,rna.seq.z.tgf,rna.seq.ztransformed,tgf.all,tgf.clust,tgf.dist.variable,tgf.up,tgf.down)
# WE only keep the rna.seq.z.profile and rna.seq


# Reading in CCLE data which we had done some manipulation with AWK
ccle.expression <- read.table("CCLE_gene_expression_awk.txt",sep="\t",header = TRUE)
# ID the cell lines that are listed as being BREAST
breast.cell.lines <- ccle.expression[,grep("BREAST",colnames(ccle.expression))]
breast.cell.lines$Description <- ccle.expression$Description #adds the gene names
others.cell.lines <- ccle.expression[,-grep("BREAST",colnames(ccle.expression))]
others.cell.lines$Description <- ccle.expression$Description
rm(ccle.expression)


rna.seq.present <- rownames(rna.seq.z.profile)[rownames(rna.seq.z.profile)%in%breast.cell.lines$Description]
# ID the genes that are common to both platforms first
rna.seq.z.profile.ccle <- rna.seq.z.profile[rownames(rna.seq.z.profile)%in%breast.cell.lines$Description,]
breast.cell.lines.expression <- breast.cell.lines[breast.cell.lines$Description%in%rna.seq.present,]
others.cell.lines.expression <- others.cell.lines[others.cell.lines$Description%in%rna.seq.present,]
others.cell.lines <- others.cell.lines[-14027] # removes the second instance of TTL
breast.cell.lines.expression <-breast.cell.lines.expression[-14027,] # removes the second instance of TTL
rm(rna.seq.present,breast.cell.lines,others.cell.lines)
others.correlation.matrix <- merge(rna.seq.z.profile.ccle,others.cell.lines.expression,by.x="Name",by.y="Description") # merges them together in a single frame for other samples
others.correlation.matrix <- others.correlation.matrix[,-3]
others.correlations <- apply(log2(others.correlation.matrix[,3:ncol(others.correlation.matrix)]),2,function(x) cor(x,y=others.correlation.matrix$Reads,method = "spearman"))
breast.correlation.matrix <- merge(rna.seq.z.profile.ccle,breast.cell.lines.expression,by.x="Name",by.y="Description") # merges them together in a single frame
breast.correlations <- apply(log2(breast.correlation.matrix[,3:ncol(breast.correlation.matrix)]),2,function(x) cor(x,y=breast.correlation.matrix$Reads,method="spearman")) # returns a list of correlations of breast cell lines
all.corrs <- c(breast.correlations,others.correlations)
rm(others.correlation.matrix,breast.correlation.matrix)

full.correlation<- data.frame(Type=c(rep("Breast",59),rep("Non-breast",length(others.correlations))),Corr=all.corrs,CellLine=c(names(breast.correlations),names(others.correlations)))# create a dataframe for plotting
rm(all.corrs)
ggplot(data=full.correlation,aes(x=factor(Type),y=Corr))+geom_boxplot()+xlab("Cell line type")+ylab("Correlation coefficient")+ggtitle("Correlation coefficients of groups")  # boxplot with ggplot
ggsave(file="~/Desktop/tcga_bioinformatics/presentations/images/correlations.png")
# get the SCGB1D2 expression
# CCLE
breast.scgb <- breast.cell.lines.expression[breast.cell.lines.expression$Description=="SCGB1D2",]
breast.ccle.scgb <- breast.scgb[-60]
others.scgb <- others.cell.lines.expression[others.cell.lines.expression$Description=="SCGB1D2",]
others.scgb <- others.scgb[-c(1:2)]
x$Expression <- as.numeric(c(brca.ccle.expression.scgb,others.ccle.expression.scgb))
ggplot(data=x,aes(x=factor(Type),y=Expression))+geom_boxplot()+xlab("Cell line type")+ylab("Expression of SCGB1D2")+ggtitle("SCGB1D2 expression of groups (CCLE)")  # boxplot with ggplot
ggsave(file="~/Desktop/tcga_bioinformatics/presentations/images/scgb1d2expression.png")
# TCGA
rna.seq.red <- rna.seq[,which(clusterassignments=="red")]
rna.seq.blue <- rna.seq[,which(clusterassignments=="blue")]
rna.seq.blue.scgb <- log2(rna.seq.blue[rownames(rna.seq.blue)=="SCGB1D2",]+1)
rna.seq.red.scgb <- log2(rna.seq.red[rownames(rna.seq.blue)=="SCGB1D2",]+1)
cluster.status <- c(rep("Blue",ncol(rna.seq.blue)),rep("Red",ncol(rna.seq.red)))
u <- data.frame(Cluster=cluster.status,RNAseq=c(unlist(as.list(rna.seq.blue.scgb)),unlist(as.list(rna.seq.red.scgb))))
v <- melt(u)
ggplot(v,aes(x=Cluster,y=value))+geom_boxplot()+ggtitle("log2 RNASeq reads for each cluster")+ylab("log2 RNASeq reads")
ggsave(file="~/Desktop/tcga_bioinformatics/presentations/images/rnaseqscgb1d2.png")
rm(u,v)
rm(rna.seq.blue,rna.seq.red,others.scgb,breast.scgb)
#CODE CHUNK 8: Getting the sample IDs in the group, and their associated CNV files
samples.in.group <- colnames(rna.seq.blue)
samples.in.group <- gsub("[.]","-",samples.in.group)  # converts it into a gsub call
manifest.files.group <- full.manifest[full.manifest$Sample%in%samples.in.group,] # gets the full list of files for samples in the group
manifest.files.CNV <- manifest.files.group[manifest.files.group$Platform.Type=="CNV (SNP Array)",]
manifest.files.CNV <- manifest.files.CNV[grep("[0-9].hg19.seg.txt",manifest.files.CNV$File.Name),] #keep only the HG19 files
library(CNTools)
data(geneInfo)
ccle.cnv <- read.table("/Volumes/HDD1/tcga/Datasets/CCLE_data/CNV_by_gene.txt",header = TRUE) # read in CCLE dataset
ccle.genes <- ccle.cnv$geneName # get the gene names
ccle.genes.info <- geneInfo[geneInfo$genename%in%ccle.genes,] # create a new gene info list for performing the reduce set
tcga.files <- as.vector(manifest.files.CNV$File.Name)
rm(manifest.files.CNV)
# set a seed first, then we incrementally add the tables to generate a huge table for performing a reduce set object with CN tools
tcga.cnv <- read.table(tcga.files[1],header=TRUE)
for (j in 2:length(tcga.files)){
    cat("Adding",tcga.files[j],"to list\n")
    f <- read.table(tcga.files[j],header=TRUE)
    tcga.cnv <- rbind(tcga.cnv,f)
    }
rm(tcga.files)
tcga.seg <- CNSeg(tcga.cnv,id="Sample",chromosome = "Chromosome",start = "Start",end="End",segMean="Segment_Mean")
tcga.reduced <- getRS(tcga.seg, by="gene",imput = FALSE, XY=FALSE,geneMap = geneInfo,what = "median")
rm(tcga.seg)
tcga.reduced.matrix <- rs(tcga.reduced) # gets the matrix of values
tcga.cnv.ccle <- tcga.reduced.matrix[tcga.reduced.matrix$genename%in%ccle.genes,] # gets only the genes also in CCLE
ccle.cnv.tcga <- ccle.cnv[ccle.cnv$geneName%in%tcga.cnv.ccle$genename,] # overlaps both
tcga.cnv.ccle <- tcga.cnv.ccle[-6979,] # removes a PRG2 which is in a difference chromosome (chr 19)
rownames(tcga.cnv.ccle) <- tcga.cnv.ccle$gene
rownames(ccle.cnv.tcga) <- ccle.cnv.tcga$geneName
tcga.cnv.ccle <- tcga.cnv.ccle[,-c(1:4)]
# to construct a CNV profile for tcga data
tcga.cnv.profile <- apply(tcga.cnv.ccle[,2:ncol(tcga.cnv.ccle)],1,mean)
tcga.cnv.profile <- data.frame(genename=tcga.cnv.ccle$genename,median=tcga.cnv.profile)
ccle.cnv.tcga <- ccle.cnv.tcga[-c(2:4)]
tcga.ccle.merged.cnv <- merge(tcga.cnv.profile,ccle.cnv.tcga,by.x="genename",by.y="geneName") # merges everything so that we can simply do an apply to compute correlation
rownames(tcga.ccle.merged.cnv) <- tcga.ccle.merged.cnv$genename
tcga.ccle.merged.cnv <- tcga.ccle.merged.cnv[,-1]
rm(tcga.cnv.profile,ccle.cnv.tcga,tcga.cnv.ccle,tcga.reduced,tcga.reduced.matrix,tcga,cnv,geneInfo)
correlations.cnv <- apply(tcga.ccle.merged.cnv[,2:ncol(tcga.ccle.merged.cnv)],2,function(x){cor(tcga.ccle.merged.cnv$median,x,method="pearson")}) #computes the correlation using the pearson correlation method
correlations.cnv <- data.frame(CellLine=names(correlations.cnv),Correlation=correlations.cnv)
#correlations.breast <- correlations.cnv[grep("BREAST",colnames(tcga.ccle.merged.cnv)[2:ncol(tcga.ccle.merged.cnv)])]
#correlations.others <- correlations.cnv[-grep("BREAST",colnames(tcga.ccle.merged.cnv)[2:ncol(tcga.ccle.merged.cnv)])]
rm(ccle.cnv,ccle.genes.info,ccle.genes,manifest.files.group)

# to clear away RNA seq data stuff to save memory
rm(breast.cell.lines.expression,rna.seq,rna.seq.z.profile,rna.seq.z.profile.ccle,others.cell.lines.expression,breast.correlation,others.correlation)
correlations.expression <- full.correlation
rm(full.correlation,full.manifest,clusterassignments,tcga.ccle.merged.cnv)


#cnvcorrelations.plot <- data.frame(Type=c(rep("Breast",length(correlationsbreast)),rep("Non-Breast",length(correlationsothers))),Cor=c(correlationsbreast,correlationsothers))
#ggplot(cnvcorrelations.plot,aes(x=factor(Type),y=Cor))+geom_boxplot()+xlab("Origin")+ylab("Correlation coefficient")+ggtitle("Correlation of CNV profile")
#ggsave(file="~/Desktop/tcga_bioinformatics/presentations/images/CNVcorrelations.png")

# Reading in the MAF files
getMAF <- function(){
    x <- read.table("batch1.maf",sep="\t",header=TRUE)
    y <- read.table("batch2.maf",sep="\t",header=TRUE)
    z <- read.table("batch3.maf",sep="\t",header=TRUE)
    mafs <- rbind(x,y,z)
    mafs
    }
allMAFs <- getMAF()
rm(getMAF)
cclemaf <- read.delim("CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",sep="\t",header = T) # reads in the latest MAF file available for use
cclemaf$detail <- paste0(cclemaf$Hugo_Symbol,"-",cclemaf$Variant_Classification)
cell.lines <- unique(cclemaf$Tumor_Sample_Barcode)  # get a list of all the cell lines
genes.present <- unique(cclemaf$Hugo_Symbol)
allMAFs <- allMAFs[allMAFs$Hugo_Symbol%in%genes.present,]
allMAFs$detailed <- paste0(allMAFs$Hugo_Symbol,"-",allMAFs$Variant_Classification)
allMAFs$Tumor_Sample_Barcode <- gsub("01[A-B]-.*$","01",allMAFs$Tumor_Sample_Barcode)
    # perform a match
x <- allMAFs[allMAFs$Tumor_Sample_Barcode%in%samples.in.group,]# this is a DF of only those in the cluster
n <- length(unique(x$Tumor_Sample_Barcode)) # this is the n
get.frequency <- function(gene){
        w <- x[x$Hugo_Symbol==gene,]
        s <- length(unique(gsub("01A-$","01A",w$Tumor_Sample_Barcode))) # extracts the no of unique samples with the mutation
        s
        } # defines subroutine
genes <- as.vector(unique(x$Hugo_Symbol)) # gets the gene list
count <- mclapply(genes,getM,mc.cores = 4)
pyx <- unlist(count)/n # this is p(y|x)
mutation.frequency <- data.frame(Gene=genes,P.Mutate=pyx,P.WT=1-pyx) # returns a dataframe
rm(pyx,genes,count,get.frequency,x,n) # clears the workspace
mutation.frequency2 <- mutation.frequency[order(mutation.frequency$P.Mutate,decreasing = TRUE),]


cell.lines.all <- as.vector(unique(cclemaf$Tumor_Sample_Barcode))
cell.lines.breast <- cell.lines.all[grep("BREAST",cell.lines.all)]

# assume that the total number of mutation is a product already present mutation and an accelerated rate of mutation in the cancer cell line
# assume that the increased rate is a function of the total mutation; that is, mu = log2N where N is the no of mutation
# we also assume that there are new random mutations that arise just by chance along time that follows a uniform distribution
# each trial is a set of bernoulli trial, where N=mu(sum(rbinom(1,1,Fi))) for i in 1:i genes
calculatebytrials <- function(x){
    cl <- cclemaf[cclemaf$Tumor_Sample_Barcode==x,] # gets the cell lines only
    mutatedgenes <- unique(cl$Hugo_Symbol)
    vecsofp <- mutation.frequency[mutation.frequency$Gene%in%mutatedgenes,]
    vecsofp <- vecsofp$P.Mutate # this is the vector of p values
    n.mutation <- length(mutatedgenes)
    mu <- log2(n.mutation) # this is the multiplier
    # we simulate over 10,000 trials each
    v=c()
    set.seed(n.mutation) #sets a seed for reproducibility
    for (i in 1:10000){
        n.r <- runif(1,0,n.mutation) # 136 is the average per cell line
        g <- rbinom(n.mutation,1,vecsofp) # single simulation
        v[i] <- (mu*sum(g))+n.r # this is the mutation estimate after the consideration of mutation acceleration
        }
    w <- ceiling(v)
    est <- ceiling(median(w))
    est
    }


getMC <- function(x){
    # x is the cell line name
    # we draw the frequency from mutation frequency
    # This is the function to get the number of mutations per cell line
    cl <- cclemaf[cclemaf$Tumor_Sample_Barcode==x,] # gets the cell lines only
    mutatedgenes <- unique(cl$Hugo_Symbol)
    vecsofp <- mutation.frequency[mutation.frequency$Gene%in%mutatedgenes,]
    vecsofp <- vecsofp$P.Mutate # this is the vector of p values
    n.mutation <- length(mutatedgenes)
    n.mutation
    }
estimate.mutations <- mclapply(cell.lines.breast,calculatebytrials,mc.cores = 4)
mutations.count <- mclapply(cell.lines.breast,getMC,mc.cores=4)
mutation <- data.frame(CellLine=cell.lines.breast,Esimated.Mutations=unlist(estimate.mutations),No.mutations=unlist(mutations.count))
P.Score <- dpois(m$No.mutations,m$Esimated.Mutations) # we model the number of mutations as a poisson distribution
mutation$P.Score <- P.Score
rm(calculatebytrials,getMC)

# to combine the scores for the breast cases
cnv <- correlations.cnv[grep("BREAST",correlations.cnv$CellLine),]
gene.expression <- correlations.expression[grep("BREAST",correlations.expression$CellLine),2:3]
scoring <- merge(cnv,gene.expression,by="CellLine")
scoring <- merge(scoring,mutation,by="CellLine")
scoring$Total <- apply(scoring[,c(2,3,6)],1,sum)
scoring$Total <- scoring$Total+2
scoring <- scoring[order(scoring$Total,decreasing = TRUE),]

