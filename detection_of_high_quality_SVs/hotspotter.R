library(primatR)
library(GenomicRanges)
library(karyoploteR)

# This script is an example of what was used to identify hotspots. Paths have
# been change for demonstration purposes. Parameters remain the same.

# Update paths and uncomment the below lines to run on your machine
# CHROM_SIZES <- "/path/to/UOA_Wagyu_1.withY.fa.fai"
# SV_BED <- "/path/to/SV.bed" # Bed file of SVs from your caller(s)
# OUT_HOTSPOT_TSV <- "/path/to/desired/save/location/hotspot.tsv"
# OUT_HOTSPOT_BED <- "/path/to/desired/save/location/hotspot.bed"

seqlens <- read.csv(CHROM_SIZES,
                    header=FALSE, sep='\t') 

# Insert column of 1s for the start position
seqlens$start <- 1

# Drip rows that start with 'mat'
seqlens <- seqlens[!grepl('mat', seqlens$V1),]
seqlens <- seqlens[!grepl('MT', seqlens$V1),]

# Create a Seqinfo object

seqinfo <- Seqinfo(seqnames=seqlens$V1, seqlengths=seqlens$V2)

# Load the SV bed file
x <- read.table(SV_BED,
col.names=c('chrom','start','end','name','score','strand'), header=FALSE, sep='\t')

# Make a GRanges object from the SV bed file dataframe
x_gr <- makeGRangesFromDataFrame(x, seqinfo = seqinfo)
x_gr <- trim(x_gr, use.names=TRUE)

# Run hotspotter
res <- hotspotter(x_gr, bw=200000, num.trial = 2000)

res_df <- data.frame(res)
write.table(res_df, file=OUT_HOTSPOT_TSV)

res_bed <- res_df[, c('seqnames','start','end','num.events','pvalue','strand')]
write.table(res_bed, file=OUT_HOTSPOT_BED,
            sep='\t',
            col.names=FALSE,
            row.names=FALSE,
            quote=FALSE)
