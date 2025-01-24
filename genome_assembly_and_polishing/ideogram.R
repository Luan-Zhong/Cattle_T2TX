## 2024.07.03
## For final review seminar and manuscript
## This Rscript does the following:
## 1. create an ideogram for each assemblies containing telomere, centromere and gaps

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)

dirin <- "/Users/polenpineda/Documents/TulixWagyu/all_Rscript_in/ideogram/"
dirout <- "/Users/polenpineda/Documents/TulixWagyu/all_Rscript_out/ideogram_2/"

## input data 
### length (input is from the plot_length.R output)
all_len <- read_delim("/Users/polenpineda/Documents/TulixWagyu/all_Rscript_out/plot_length/all_chr_length.txt")
all_len <- all_len %>%
  mutate(seqnames = paste0(asm,"_",chr))
### gaps
gap_wagyu <- read_delim(paste0(dirin,"UOA_Wagyu_1.coor"), col_names = FALSE)
gap_tuli <- read_delim(paste0(dirin,"UOA_Tuli_1.coor"), col_names = FALSE)
gap_ars <- read_delim(paste0(dirin,"ARS_UCD_v2.0.coor"), col_names = FALSE)
gap_wagyu$seqnames <- paste0("Wagyu_", gap_wagyu$X1)
gap_tuli$seqnames <- paste0("Tuli_", gap_tuli$X1)
gap_ars$seqnames <- paste0("ARS_UCD2.0_", gap_ars$X1)
all_gap <- rbind(gap_wagyu, gap_tuli, gap_ars) 
names(all_gap) <- c("chr","start","end","seqnames")
all_gap <- all_gap %>%
  filter(chr %in% c(paste(1:29), "X", "Y"))

### centromeres are from the centromere_boundary.R output

ctr_wagyu <- read_delim("/Users/polenpineda/Documents/TulixWagyu/all_Rscript_out/final_repeats_2/wagyu_centromere_location.txt")
ctr_ars <- read_delim("/Users/polenpineda/Documents/TulixWagyu/all_Rscript_out/centromere_boundary/ars_ctr_location_2.txt")
ctr_wagyu$seqnames <- paste0("Wagyu_", ctr_wagyu$chromosome)
# ctr_tuli$seqnames <- paste0("Tuli_", ctr_tuli$chromosome)
ctr_ars$seqnames <- paste0("ARS_UCD2.0_", ctr_ars$chromosome)
# all_ctr <- rbind(ctr_wagyu, ctr_tuli, ctr_ars)

stop()
## Plot ideogram for the assemblies with telomeres and gaps
## idea here https://www.earlham.ac.uk/articles/marking-telomeres-simple-ideogram-r

### create ideogram for per assemblies with telomere, centromere and gaps
## create factors for every chr
all_len$seqnames <- factor(
  all_len$seqnames, levels = unique(all_len$seqnames), ordered = TRUE
)

levels_seqnames <- unique(c(levels(all_len$seqnames)))

for (i in c("ARS_UCD2.0_", "Wagyu_", "Tuli_")) {
  asm_level <- levels_seqnames[grep(i, levels_seqnames)]
  asm_len <- all_len %>% filter(grepl(i, seqnames)) %>% mutate(seqnames = factor(seqnames, levels = asm_level))
  asm_tlm <- all_tlm %>% filter(grepl(i, seqnames)) %>% mutate(seqnames = factor(seqnames, levels = asm_level))
  asm_ctr <- all_ctr %>% filter(grepl(i, seqnames)) %>% mutate(seqnames = factor(seqnames, levels = asm_level))
  asm_gap <- all_gap %>% filter(grepl(i, seqnames)) %>% mutate(seqnames = factor(seqnames, levels = asm_level))
  tlm <- asm_tlm %>%
    gather(telomere, count, p_count, q_count) %>%
    mutate(start = ifelse(telomere == "p_count" & !is.na(count), 1, length-(count*6)),
           end = ifelse(telomere == "p_count" & !is.na(count), 1+(count*6), length))
  
  asm_gff <- tlm %>% select(seqnames, start, end) %>% mutate(type = "telomere") %>% na.omit()
  asm_gff <- asm_gff %>% mutate(seqnames=factor(seqnames, levels=levels(asm_gff$seqnames)))
  
  plot <- ggplot(asm_len, aes(x=factor(seqnames), y=length)) +
    geom_rect(aes(ymax=length+1),
              ymin=1,
              xmin=as.numeric(asm_len$seqnames)-0.2,
              xmax=as.numeric(asm_len$seqnames)+0.2,
              fill="gray",
              colour="gray") +
    scale_y_continuous(
      limits=c(0, ceiling(max(asm_len$length)/2e8)*2e8),
      labels=label_number(
        accuracy=1,
        scale=1e-6,
        suffix="Mbp")
    ) +
    coord_flip(clip="off") +
    theme(axis.text.y=element_text(
      colour="black",
      size=8,
      margin=margin(r=10)
    ),
    axis.text.x=element_text(size=7),
    axis.ticks.y=element_blank(),
    axis.title=element_blank(),
    axis.line.x=element_line(),
    panel.grid.major=element_blank(), 
    panel.grid.minor=element_blank(), 
    panel.background=element_blank())
  
  plot_fin <- plot +
    geom_rect(data=asm_ctr,
              aes(ymin=start, ymax=end, xmin=as.numeric(seqnames)-0.2, xmax=as.numeric(seqnames)+0.2),
              fill="steelblue",
              colour="steelblue",
              inherit.aes=FALSE) +
    geom_rect(data=asm_gff,
              aes(ymin=start, ymax=end, xmin=as.numeric(seqnames)-0.3, xmax=as.numeric(seqnames)+0.3),
              fill="red",
              colour="red",
              inherit.aes=FALSE) +
    geom_rect(data=asm_gap,
              aes(ymin=start, ymax=end, xmin=as.numeric(seqnames)-0.2, xmax=as.numeric(seqnames)+0.2),
              fill="black",
              colour="black",
              inherit.aes=FALSE)
  filename <- paste0(dirout,"ideogram_asm_", i, "plot.png")
  ggsave(filename, plot_fin, width = 7, height = 4, units = "in", dpi = 300)
}
stop()

### get the gene density
require(RIdeogram)

wagyu_karyotype <- read.table("/Users/polenpineda/Documents/TulixWagyu/RIdeogram/input_files/wagyu_karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
wagyu_gene_density <- GFFex(input = "/Users/polenpineda/Documents/TulixWagyu/RIdeogram/input_files/Bos_taurus_gca040286185v1.ASM4028618v1.110.chr.gff3.gz", karyotype = "/Users/polenpineda/Documents/TulixWagyu/RIdeogram/input_files/wagyu_karyotype.txt", feature = "gene", window = 1000000)
wagyu_gene_marker <- read.table("/Users/polenpineda/Documents/TulixWagyu/RIdeogram/input_files/gene_marker.txt", sep = "\t", header = T, stringsAsFactors = F)

wagyu_karyotype$Chr <- factor(
  wagyu_karyotype$Chr, levels = unique(wagyu_karyotype$Chr), ordered = TRUE
)
levels_seqnames <- unique(c(levels(wagyu_karyotype$Chr)))

wagyu_len <- wagyu_karyotype  %>% mutate(Chr = factor(Chr, levels = levels_seqnames))
wagyu_tlm <- wagyu_gene_marker %>% filter(Type == "telomere") %>% mutate(Chr = factor(Chr, levels = levels_seqnames))
wagyu_ctr <- wagyu_gene_marker %>% filter(Type == "centromere") %>% mutate(Chr = factor(Chr, levels = levels_seqnames))
wagyu_gap <- wagyu_gene_marker %>% filter(Type == "gap") %>% mutate(Chr = factor(Chr, levels = levels_seqnames))
wagyu_gene <- wagyu_gene_density %>% mutate(Chr = factor(Chr, levels = levels_seqnames))

### Create the base
plot_wagyu <- ggplot(wagyu_len, aes(x=factor(Chr), y=End)) +
  geom_rect(aes(ymax=End+1),
            ymin=1,
            xmin=as.numeric(wagyu_len$Chr)-0.2,
            xmax=as.numeric(wagyu_len$Chr)+0.2,
            fill="lightgray",
            colour="lightgray") +
  scale_y_continuous(
    limits=c(0, ceiling(max(wagyu_len$End)/2e8)*2e8),
    labels=label_number(
      accuracy=1,
      scale=1e-6,
      suffix=" Mb")
  ) +
  coord_flip(clip="off") +
  theme(axis.text.y=element_text(
    colour="black",
    size=8,
    margin=margin(r=10)
  ),
  axis.text.x=element_text(size=7),
  axis.ticks.y=element_blank(),
  axis.title=element_blank(),
  axis.line.x=element_line(),
  panel.grid.major=element_blank(), 
  panel.grid.minor=element_blank(), 
  panel.background=element_blank())
plot_wagyu

wagyu_gene_density$Chr <- ifelse(wagyu_gene_density$Chr == "X", 30, as.numeric(wagyu_gene_density$Chr))

plot_wagyu_gene <- plot_wagyu +
  geom_rect(data=wagyu_gene_density,
            aes(ymin=Start, ymax=End, xmin=as.numeric(Chr)-0.22, xmax=as.numeric(Chr)+(Value/400)),
            fill="darkgray",
            inherit.aes=FALSE)
plot_wagyu_gene

plot_wagyu_fin <- plot_wagyu +
  geom_rect(data=wagyu_ctr,
            aes(ymin=Start, ymax=End, xmin=as.numeric(Chr)-0.2, xmax=as.numeric(Chr)+0.2),
            fill="royalblue",
            alpha=0.5,
            inherit.aes=FALSE)+
  geom_rect(data=wagyu_gene_density,
            aes(ymin=Start, ymax=End, xmin=as.numeric(Chr)-0.22, xmax=as.numeric(Chr)+(Value/400)),
            fill="forestgreen",
            alpha=0.5,
            inherit.aes=FALSE)+
  geom_rect(data=wagyu_tlm,
            aes(ymin=Start, ymax=End, xmin=as.numeric(Chr)-0.3, xmax=as.numeric(Chr)+0.3),
            fill="red",
            colour="red",
            inherit.aes=FALSE) +
  geom_rect(data=wagyu_gap,
            aes(ymin=Start, ymax=End, xmin=as.numeric(Chr)-0.2, xmax=as.numeric(Chr)+0.2),
            fill="black",
            colour="black",
            inherit.aes=FALSE)
plot_wagyu_fin

filename <- paste0(dirout,"ideogram_wagyu_final_plot.png")
ggsave(filename, plot_fin, width = 7, height = 4, units = "in", dpi = 300)

plot_wagyu

ars_fai <- "/Users/polenpineda/Documents/TulixWagyu/reference_genome/ARS-UCD2.0_frozen.fna.fai"
ars_fai <- read_delim(ars_fai, col_names = c("Chr","End"))
ars_karyotype <- ars_fai %>% filter(Chr %in% c(1:29,"X","Y")) %>% mutate(Start = 0) %>% select(Chr, Start, End)

## additional bases

additional_bases <- merge(wagyu_karyotype, ars_karyotype, by = "Chr")
additional_bases <- additional_bases %>% mutate(additional_bases = End.x - End.y)

plot_additional_bases <- additional_bases %>%
  ggplot(aes(x=additional_bases/1e6, y=Chr)) +
  geom_bar(stat = "identity", fill = "gray") +
  theme_classic() +
  labs(y = NULL, x = "Size (Mb)") +
  theme(legend.position = "none")
plot_additional_bases

wagyu_karyotype$asm <- "UOA_Wagyu_1"
ars_karyotype$asm <- "ARS-UCD2.0"
additional_bases_bind <- rbind(wagyu_karyotype, ars_karyotype)
Y_row <- data.frame(Chr = "Y", Start = 0, End = 0, asm = "UOA_Wagyu_1")
additional_bases_bind <- rbind(additional_bases_bind, Y_row)


plot_additional_bases_bind <- additional_bases_bind %>%
  ggplot(aes(x = End / 1e6, y = Chr, fill = asm)) +  # Fill by 'asm' variable
  geom_bar(stat = "identity", position = "dodge") +  # Use position = "dodge" for side-by-side bars
  theme_classic() +
  labs(y = NULL, x = "Size (Mb)", fill = "Assembly") +
  scale_fill_manual(values = c(
    "ARS-UCD2.0" = "#F88379",
    "UOA_Wagyu_1" = "#0892d0"
  ))
plot_additional_bases_bind
