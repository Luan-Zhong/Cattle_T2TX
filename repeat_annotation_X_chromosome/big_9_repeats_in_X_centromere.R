library(readr)
library(dplyr)
library(ggplot2)
library(scales)

dir <- "/Users/polenpineda/Documents/TulixWagyu/centromere/x_chromosome/big_9_repeats/Rscript_out/"
#blast <- "/Users/polenpineda/Documents/TulixWagyu/centromere/x_chromosome/big_9_repeats/UOA_Wagyu_1_vs_X_39560000-39950000.fa.85.blast.tsv"
blast <- "/Users/polenpineda/Documents/TulixWagyu/centromere/x_chromosome/big_9_repeats/UOA_Wagyu_1_vs_X_318Kbp.fa.85.blast.tsv"

blast_df <- read_delim(blast, comment = "#", col_names = c("qry", "ref", "perc_id", "align_len",
                                                                       "mismatches", "gap_opens", "qry_start", "qry_end",
                                                                       "ref_start", "ref_end", "eval", "bitscore"))

# filter 95% percent ID
blast_df_filter_95 <- blast_df %>%
  filter(perc_id > 95, ref == "X")

blast_df_filter_100 <- blast_df %>%
  filter(perc_id == 100, ref == "X")

xctr <- data.frame(start = 38e6, end = 50e6)

# alignment coverage within the centromeric region
blast_df_filter_95_ctr <- blast_df %>%
  mutate(new_align_len = ref_end - ref_start,
         tmp_start = ifelse(ref_start > ref_end, ref_end, ref_start),
         tmp_end = ifelse(ref_start > ref_end, ref_start, ref_end),
         oldref_start = ref_start,
         oldref_end = ref_end,
         ref_start = tmp_start,
         ref_end = tmp_end) %>%
  filter(perc_id > 95, ref == "X",
         ref_start > xctr$start, ref_end < xctr$end) %>%
  arrange(ref_start)

x_axis <- 1:390e3
counts <- sapply(x_axis, function(x) sum(blast_df_filter_95_ctr$qry_start <= x & blast_df_filter_95_ctr$qry_end >= x))
cov_blast_df_filter_95_ctr <- data.frame(x = x_axis, count = counts)
coverage_390kb_95 <- ggplot(cov_blast_df_filter_95_ctr, aes(x = x, y = count)) +
  geom_line() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = comma) +
  labs(x = "position (bp)", y = "coverage (bp)",
       title = "Coverage of the aligned 390Kb to the x centromere region",
       caption = "Blast result with sequence identity filtered at 95%")
coverage_390kb_95
ggsave(filename = "coverage_aligned_390Kb_x_ctr.jpg", plot = coverage_390kb_95, path = dir, width = 8, height = 4, units = "in", dpi = 300)


## create a histogram/density of the hit sizes

histo_blast_95_density_result <- density(blast_df_filter_95_ctr$align_len, bw = .01)
print(histo_blast_95_density_result)
histo_blast_95_density_peak <- histo_blast_95_density_result$x[which.max(histo_blast_95_density_result$y)]
print(histo_blast_95_density_peak)

histo_blast_95 <- ggplot(blast_df_filter_95_ctr, aes(x = align_len)) +
  geom_density(bw = 0.01) +
  scale_x_continuous(labels = comma) +
  labs(title = paste0("Density plot peak at ",histo_blast_95_density_peak),
       x = "Alignment length",
       y = "Density")
histo_blast_95



stop()

## Focusing on the bigger ones for now
## filtering the alignment length to be greater than 10k
big_blast_df_filter_95_ctr <- blast_df %>%
  mutate(new_align_len = ref_end - ref_start) %>%
  filter(perc_id > 95, ref == "X", align_len > 10e3,
         ref_start > xctr$start, ref_end < xctr$end) %>%
  arrange(ref_start)

mean(big_blast_df_filter_95_ctr$perc_id)
  

x_axis <- 1:318440
counts <- sapply(x_axis, function(x) sum(big_blast_df_filter_95_ctr$qry_start <= x & big_blast_df_filter_95_ctr$qry_end >= x))
cov_big_blast_df_filter_95_ctr <- data.frame(x = x_axis, count = counts)
coverage_390kb_95_10kalign <- ggplot(cov_big_blast_df_filter_95_ctr, aes(x = x, y = count)) +
  geom_line() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(limits = c(1, 15), breaks = c(1:15)) +
  labs(x = "Position (bp)", y = "Coverage (bp)",
       caption = "Blast result with sequence identity filtered at 95% and >10% alignment cover")
coverage_390kb_95_10kalign
ggsave(filename = "coverage_aligned_390Kb_x_ctr_10k_align.jpg", plot = coverage_390kb_95_10kalign, path = dir, width = 8, height = 4, units = "in", dpi = 300)


location_big_9 <- blast_df

## i should look at the interval of the hits

big_blast_df_filter_95_ctr <- big_blast_df_filter_95_ctr %>%
  arrange(ref_start) %>%
  mutate(tmp_start = ifelse(ref_start > ref_end, ref_end, ref_start),
         tmp_end = ifelse(ref_start > ref_end, ref_start, ref_end)) %>% 
  mutate(diff = tmp_start - lag(tmp_end), diff_qry = qry_start - lag(qry_end))


### How about the repeatmodeler repeats, what are within this big repeat?

blast_consensei <- "/Users/polenpineda/Documents/TulixWagyu/centromere/x_chromosome/big_9_repeats/X_39560000-39950000_vs_consensi.classified.autosome.sats.fa.85.blast.renamed.tsv"
blast_consensei_df <- read_delim(blast_consensei, comment = "#", col_names = c("qry", "ref", "perc_id", "align_len",
                                                           "mismatches", "gap_opens", "qry_start", "qry_end",
                                                           "ref_start", "ref_end", "eval", "bitscore"))
blast_consensei_df_filter <- blast_consensei_df %>%
  filter(perc_id >= 95)
## look at the alignment coverage?
x_axis <- 1:390e3
counts <- sapply(x_axis, function(x) sum(blast_consensei_df_filter$qry_start <= x & blast_consensei_df_filter$qry_end >= x))
cov_blast_consensei_df_filter <- data.frame(x = x_axis, count = counts)
coverage_390kb_95_consensei <- ggplot(cov_blast_consensei_df_filter, aes(x = x, y = count)) +
  geom_line() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = comma) +
  labs(x = "position (bp)", y = "coverage (bp)",
       title = "Coverage of the aligned consensus repeats from repmod to the 390Kbp",
       caption = "Blast result with sequence identity filtered at 95%")
coverage_390kb_95_consensei
ggsave(filename = "coverage_aligned_390Kb_consensei.jpg", plot = coverage_390kb_95_consensei, path = dir, width = 8, height = 4, units = "in", dpi = 300)

## seems like repeatmodeler missed this big repeat

## self alignment if there are repeats inside

blast_318_self <- "/Users/polenpineda/Documents/TulixWagyu/centromere/x_chromosome/big_9_repeats/X_318Kbp_vs_X_318Kbp.fa.85.blast.tsv"
blast_318_self_df <- read_delim(blast_318_self, comment = "#", col_names = c("qry", "ref", "perc_id", "align_len",
                                                                               "mismatches", "gap_opens", "qry_start", "qry_end",
                                                                               "ref_start", "ref_end", "eval", "bitscore"))
blast_318_self_df_filter <- blast_318_self_df %>%
  filter(perc_id >= 95)

x_axis <- 1:318440
counts <- sapply(x_axis, function(x) sum(blast_318_self_df_filter$qry_start <= x & blast_318_self_df_filter$qry_end >= x))
cov_blast_318_self_df_filter <- data.frame(x = x_axis, count = counts)
coverage_blast_318_self_df_filter <- ggplot(cov_blast_318_self_df_filter, aes(x = x, y = count)) +
  geom_line() +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(limits = c(1, 6), breaks = c(1:6)) +
  labs(x = "position (bp)", y = "coverage (bp)",
       title = "Coverage of the aligned consensus repeats from repmod to the 390Kbp",
       caption = "Blast result with sequence identity filtered at 95%")
coverage_blast_318_self_df_filter
ggsave(filename = "coverage_self_318Kb_consensei.jpg", plot = coverage_blast_318_self_df_filter, path = dir, width = 8, height = 4, units = "in", dpi = 300)

## what are the sizes of these repeats?


summary_cov_blast_318_self_df_filter <- cov_blast_318_self_df_filter %>%
  group_by(count) %>%
  summarise(freq = n())

### inverted repeat within the 318Kbp sequence


### create a gff file from the blast results

gff_blast <- blast_df_filter_95_ctr %>%
  mutate(chr = ref,
         source = "blast",
         start = ref_start,
         end = ref_end,
         len = ref_end - ref_start,
         score = perc_id,
         strand = ifelse(len < 0, "-", "+"),
         frame = "0",
         attribute = ".") %>%
  select(chr, source, qry, start, end, score, strand, frame, attribute)
write_tsv(gff_blast, file = paste0(dir, "x_ctr_380kbp_blast.gff"), col_names = FALSE)


### is it found in another chromosomes?
# the big chunk is not found in other chromosomes
blastf_df_auto_x <- blast_df %>%
  filter(perc_id > 95, align_len > 10e3, qry_start > 56526, qry_end < 374965)

blastf_df_auto_x <- blast_df %>%
  filter(perc_id > 99, qry_start > 56526, qry_end < 374965)

blastf_df_auto_x_summary <- blastf_df_auto_x %>%
  group_by(ref) %>%
  summarise(count = n(),
            perc_id_mean = mean(perc_id),
            perc_id_sd = sd(perc_id),
            align_len_mean = mean(align_len),
            align_len_min = mean(align_len),
            align_len_max = max(align_len))

