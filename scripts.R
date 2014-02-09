# CODE CHUNK 1: Name parsing- to check which samples are going to be used

batch1 <- read.table("batch1_manifest.txt",header = T,sep="\t")
batch2 <- read.table("batch2_manifest.txt",header = T,sep="\t")
batch3 <- read.table("batch3_manifest.txt",header = T,sep="\t")
manifest_all <- rbind(batch1,batch2,batch3)
manifest_all <- manifest_all[manifest_all$Platform.Type!="METADATA",]

# to perform subsetting and splitting, and to determine which are the needed files
# all platform types
avail.types <- as.vector(unique(manifest_all$Platform))
# all IDs
sele.sample <- function(){
    lapply(avail.types[1:3],function(x){g <- manifest_all[manifest_all$Platform==x,];return(unique(g$Sample))})
    }
sample.by.type <- sele.sample()
names(sample.by.type) <- avail.types[1:3]

# to verify which IDs have data for both CNV and RNA-Seq
cnv.data <- sample.by.type$Genome_Wide_SNP_6
rnaseqv1.data <- sample.by.type$IlluminaHiSeq_RNASeq
rnaseqv2.data <- sample.by.type$IlluminaHiSeq_RNASeqV2
cnv.data.seq1 <- as.vector(cnv.data[cnv.data%in%rnaseqv1.data]) # cnv data with RNA seq 1
cnv.data.seq2 <- as.vector(cnv.data[cnv.data%in%rnaseqv2.data]) # cnv data with RNA seq 2
cnv.with.seq <- unique(c(cnv.data.seq1,cnv.data.seq2)) # returns all the IDs with both CNV data and accompanying RNA seq (either 1 or 2) data; for id-ing which are present in somatic mutations
# all the ones with seq and CNV data are with RNA seq 2 as well, so we don't need the RNASeqV1 data
# to copy thus only all the files for RNAseq2 and CNV and mutation into one single folder (for the ease of reading in the files later)

# CODE CHUNK 2: TO ISOLATE ALL THE FILES TO BE READ IN- remove all the RNASeq1 data first
manifest.in <- manifest_all[manifest_all$Platform!="IlluminaHiSeq_RNASeq",]
manifest.in <- manifest.in[manifest.in$Sample%in%cnv.with.seq,]
files.in <- manifest.in$File.Name
files.in <- as.vector(files.in[-grep("exon|isoform|junction|genes.normalize|hg18",files.in)])
manifest.final <- manifest.in[manifest.in$File.Name%in%files.in,]

# CODE CHUNK 3 : To read in the files
object.names <- paste0(manifest.final$Platform.Type,"-",manifest.final$Sample)
read.files <- function(fn,ol){
    f.in <- read.table(file=fn)
    assign(ol,f.in,envir = globalenv())
    return(NULL)
    }
mapply(read.files,fn=files.in,ol=object.names)




















