# CODE CHUNK 1: Making the required manifest files to read files etc;
# require(multicore) # for parallelizations

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
rna.seq.z.profile <- data.frame(Name=rownames(rna.seq),Reads=rna.seq.z.profile)
rownames(rna.seq.z.profile) <- rna.seq.z.profile$Name
# CODE CHUNK 5: Read in CCLE gene expression data
ccle.expression <- read.table("CCLE_gene_expression_awk.txt",sep="\t",header = TRUE)
# ID the cell lines that are listed as being BREAST
brca.cell.lines <- ccle.expression[,grep("BREAST",colnames(ccle.expression))]
brca.cell.lines$Description <- ccle.expression$Description
others.cell.lines <- ccle.expression[,-grep("BREAST",colnames(ccle.expression))]
others.cell.lines$Description <- ccle.expression$Description

# CODE CHUNK 6: Comparing the CCLE and TCGA data using the spearman correlation coefficient
rna.seq.present <- rownames(rna.seq.z.profile)[rownames(rna.seq.z.profile)%in%ccle.expression$Description]
rna.seq.z.profile.ccle <- rna.seq.z.profile[rownames(rna.seq.z.profile)%in%ccle.expression$Description,]# ID the genes that are common to both platforms first
ccle.breast.expression <- brca.cell.lines[brca.cell.lines$Description%in%rna.seq.present,]
others.cell.lines.expression <- others.cell.lines[others.cell.lines$Description%in%rna.seq.present,]
others.cell.lines <- others.cell.lines[-14027]
ccle.breast.expression <- ccle.breast.expression[-14027,] # removes the second instance of TTL
others.correlation.matrix <- merge(rna.seq.z.profile.ccle,others.cell.lines.expression,by.x="Name",by.y="Description") # merges them together in a single frame for other samples
others.correlation.matrix <- others.correlation.matrix[,-3]
others.correlations <- apply(log2(others.correlation.matrix[,3:ncol(others.correlation.matrix)]),2,function(x) cor(x,y=others.correlation.matrix$Reads,method = "spearman"))
correlation.matrix <- merge(rna.seq.z.profile.ccle,ccle.breast.expression,by.x="Name",by.y="Description") # merges them together in a single frame
correlations <- apply(log2(correlation.matrix[,3:ncol(correlation.matrix)]),2,function(x) cor(x,y=correlation.matrix$Reads,method="spearman")) # returns a list of correlations
all.corrs <- c(correlations,others.correlations)

# CODE CHUNK 7: plotting boxplot of correlation coefficients with ggplot
require(ggplot2)
x <- data.frame(Type=c(rep("Breast",59),rep("Non-breast",length(others.correlations))),Corr=all.corrs)# create a dataframe for plotting
ggplot(data=x,aes(x=factor(Type),y=Corr))+geom_boxplot()+xlab("Cell line type")+ylab("Correlation coefficient")+ggtitle("Correlation coefficients of groups")  # boxplot with ggplot
ggsave(file="~/Desktop/tcga_bioinformatics/presentations/images/correlations.png")
# get the SCGB1D2 expression
# CCLE
brca.ccle.expression.scgb <- brca.cell.lines[brca.cell.lines$Description=="SCGB1D2",]
brca.ccle.expression.scgb <- brca.ccle.expression.scgb[-60]
others.ccle.expression.scgb <- others.cell.lines.expression[others.cell.lines$Description=="SCGB1D2",]
others.ccle.expression.scgb <- others.ccle.expression.scgb[-c(1:2)]
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

#CODE CHUNK 8: Getting the sample IDs in the group, and their associated CNV files
samples.in.group <- colnames(rna.seq.blue)
samples.in.group <- gsub("[.]","-",samples.in.group)  # converts it into a gsub call
manifest.files.group <- full.manifest[full.manifest$Sample%in%samples.in.group,] # gets the full list of files for samples in the group
manifest.files.CNV <- manifest.files.group[manifest.files.group$Platform.Type=="CNV (SNP Array)",]
manifest.files.CNV <- manifest.files.CNV[grep("[0-9].hg19.seg.txt",manifest.files.CNV$File.Name),] #keep only the HG19 files

data(geneInfo)
ccle.cnv <- read.table("/Volumes/HDD1/tcga/Datasets/CCLE_data/CNV_by_gene.txt",header = TRUE) # read in CCLE dataset
ccle.genes <- ccle.cnv$geneName # get the gene names
ccle.genes.info <- geneInfo[geneInfo$genename%in%ccle.genes,] # create a new gene info list for performing the reduce set
require(CNTools)
require(GenomicRanges)
tcga.files <- as.vector(manifest.files.CNV$File.Name)
# set a seed first, then we incrementally add the tables to generate a huge table for performing a reduce set object with CN tools

for (j in 2:length(tcga.files)){
    cat("Adding",tcga.files[j],"to list\n")
    f <- read.table(tcga.files[j],header=TRUE)
    tcga.cnv <- rbind(tcga.cnv,f)
    }
tcga.seg <- CNSeg(tcga.cnv,id="Sample",chromosome = "Chromosome",start = "Start",end="End",segMean="Segment_Mean")
tcga.reduced <- getRS(tcga.seg, by="gene",imput = FALSE, XY=FALSE,geneMap = geneInfo,what = "median")
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
correlations.cnv <- apply(tcga.ccle.merged.cnv[,2:ncol(tcga.ccle.merged.cnv)],2,function(x){cor(tcga.ccle.merged.cnv$median,x,method="pearson")}) #computes the correlation using the pearson correlation method
correlationsbreast <- correlations.cnv[grep("BREAST",colnames(tcga.ccle.merged.cnv)[2:ncol(tcga.ccle.merged.cnv)])]
correlationsothers <- correlations.cnv[-grep("BREAST",colnames(tcga.ccle.merged.cnv)[2:ncol(tcga.ccle.merged.cnv)])]
cnvcorrelations.plot <- data.frame(Type=c(rep("Breast",length(correlationsbreast)),rep("Non-Breast",length(correlationsothers))),Cor=c(correlationsbreast,correlationsothers))
ggplot(cnvcorrelations.plot,aes(x=factor(Type),y=Cor))+geom_boxplot()+xlab("Origin")+ylab("Correlation coefficient")+ggtitle("Correlation of CNV profile")
ggsave(file="~/Desktop/tcga_bioinformatics/presentations/images/CNVcorrelations.png")

# note: Code chunk for running CNV landscape was not included, but can be built from the workspace object that is loaded via load("workspace.Rdata")

# CODE CHUNK 9: Working with somatic mutations
# using level 2 MAF (uncurated data)

require(multicore) # for mclapply call in the calculation of priors to speed up computation
getMAF <- function(){
    x <- read.table("batch1.maf",sep="\t",header=TRUE)
    y <- read.table("batch2.maf",sep="\t",header=TRUE)
    z <- read.table("batch3.maf",sep="\t",header=TRUE)
    mafs <- rbind(x,y,z)
    mafs
    }
allMAFs <- getMAF()
rm(getMAF)
allMAFs <- allMAFs[with(allMAFs,order(Chrom)),] # orders the chromosome number
allMAFs$detailed <- paste0(allMAFs$Hugo_Symbol,"-",allMAFs$Variant_Classification)
calPy <- function(x){
    # get the genes
    genes <- as.vector(unique(allMAFs$detailed))
    cat("Priors for",length(genes),"genes will be estimated\n")
    # perform the summarizations
    n <- 935 # this is the number of samples
    # define subroutine to get mutations by the genes
    getM <- function(gene){
        #cat("Estimating prior for ",gene,"\n")
        w <- x[x$detailed==gene,]
        s <- length(unique(gsub("01A-$","01A",w$Tumor_Sample_Barcode))) # extracts the no of unique samples with the mutation
        s
        }
    counts <- mclapply(genes,getM,mc.cores = 4)
    py <- unlist(counts)/n
    res <- data.frame(Gene=genes,PriorFreq=py)
    return(res)
    }
# to calculate p(y|x)
# getMAF for only those in the cluster
getMAF <- function(){
    #llMAFs$Tumor_Sample_Barcode <- gsub("01[A-B]-.*$","01",allMAFs$Tumor_Sample_Barcode)
    # perform a match
    x <- allMAFs[allMAFs$Tumor_Sample_Barcode%in%samples.in.group,]# this is a DF of only those in the cluster
    n <- length(unique(x$Tumor_Sample_Barcode)) # this is the n
    getM <- function(gene){
        w <- x[x$detailed==gene,]
        s <- length(unique(gsub("01A-$","01A",w$Tumor_Sample_Barcode))) # extracts the no of unique samples with the mutation
        s
        } # defines subroutine
    genes <- as.vector(unique(x$detailed)) # gets the gene list
    cat("Calculating prior for ",length(genes),"genes \n")
    count <- mclapply(genes,getM,mc.cores = 4)
    pyx <- unlist(count)/n # this is p(y|x)
    res <- data.frame(Gene=genes,Prior=pyx) # returns a dataframe
    res
    }

# THESE ARE THE PARAMETERS FOR COMPUTATION OF CONDITIONAL PROBABILITY #
pYX <- getMAF() # gets the P(Y|X) for the genes //
pY <- 568/936 # this is the probability of being TGFb sensitive
priors <- calPy(allMAFs) # this is a table of priors for all gene-mutation type combination # p(x)
#######################################################################
cclemaf <- read.table("CCLE_Oncomap3.maf",sep="\t",header=TRUE)
calculatePOST <- function(cellline){
    cat("Computing for cell line",cellline,"\n")
    cclemaf2 <- cclemaf[cclemaf$Tumor_Sample_Barcode==cellline,]
    mtype <- unique(cclemaf2$detailed)  # this specifies the mutationtype
    #print(mtype)
    x <- priors[priors$Gene%in%mtype,] # gets all the probabilities px -> probability of mutation occuring
    #print(x)
    if(nrow(x)==0){
        px <- 0
        } else
            {px <- prod(x$PriorFreq)}
    yx <- pYX[pYX$Gene%in%mtype,] # gets pY|X -> prob of mutation occuring given its a mutant

    if(nrow(yx)==0){
        pyx <- 0
        } else
            {  pyx <- prod(yx$Prior)
  }

    #print(pyx)
    py <- pY
    #print(py)
    post <- pyx*px/py
    cat("P(Y|X)=",pyx,",PX=",px,"Score=",post,"\n")
    return(post) #returns the posterior distribution
    }
cclemaf$detailed <- paste0(cclemaf$Hugo_Symbol,"-",cclemaf$Variant_Classification)
cls <- as.vector(unique(cclemaf$Tumor_Sample_Barcode))
mutationscores
mutationscores <- mclapply(cls,calculatePOST,mc.cores=4)
# These are the scoring components
correlations.cnv <- as.data.frame(correlations.csv)# object that stores CNV correlations for all the cell lines # already in workspace
correlations.cnv$CL <- rownames(correlations.cnv)
all.corrs <- as.data.frame(all.corrs) # this is the
all.corrs$CL <- rownames(all.corrs)
mutationscores <- data.frame(CL=cls,MScores=unlist(mutationscores))
scoringmatrix <- merge(correlations.cnv,all.corrs,"CL","CL")
scoringmatrix <- merge(scoringmatrix,mutationscores,"CL","CL")
scoringmatrix$TotalScore <- apply(scoringmatrix[,2:4],1,function(x)sum(x)+2)
# Code Chunk 10: Selection of cell lines and representing them
# probably more visualizations
