library(primatR)
library(GenomicRanges)
library(karyoploteR)

# This script is an example of what was used to identify hotspots. Paths have
# been change for demonstration purposes. Parameters remain the same.

seqlens <- read.csv('/Users/callummacphillamy/Projects/REFERENCES/UOA_WAGYU/UOA_Wagyu_1.withY.fa.fai',
                    header=FALSE, sep='\t') 
# Insert column of 1s for the start position
seqlens$start <- 1

# Drip rows that start with 'mat'
seqlens <- seqlens[!grepl('mat', seqlens$V1),]
seqlens <- seqlens[!grepl('MT', seqlens$V1),]

# Create a Seqinfo object

seqinfo <- Seqinfo(seqnames=seqlens$V1, seqlengths=seqlens$V2)


SV_BED <- "/path/to/SV.bed"

x <- read.table(SV_BED,
col.names=c('chrom','start','end','name','score','strand'), header=FALSE, sep='\t')

x_gr <- makeGRangesFromDataFrame(x, seqinfo = seqinfo)
x_gr <- trim(x_gr, use.names=TRUE)
res <- hotspotter(x_gr, bw=200000, num.trial = 2000)

res_df <- data.frame(res)
write.table(res_df, file="/path/to/SV_hotspots.tsv")

res_bed <- res_df[, c('seqnames','start','end','num.events','pvalue','strand')]
write.table(res_bed, file="/path/to/SV_hotspots.bed",
            sep='\t',
            col.names=FALSE,
            row.names=FALSE,
            quote=FALSE)
