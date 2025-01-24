## PSPineda
## Count the number of rDNA copies from blast result

library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(tidyr)

get_first_digit <- function(x) {
  x <- abs(x)
  floor(x / 10^(floor(log10(x))))
}

dirout <- "/Users/polenpineda/Documents/TulixWagyu/rdna/rdna_35Kb/Rscript_output_2"
## INPUT
## blast
blast_rdna <- "/Users/polenpineda/Documents/TulixWagyu/rdna/rdna_35Kb/UOA_Wagyu_1_vs_rdna_rep_chr11_split.fa.85.blast.tsv"
blast_rdna_df <- read_delim(blast_rdna, comment = "#", col_names = c("qry", "ref", "perc_id", "align_len",
                                                                     "mismatches", "gap_opens", "qry_start", "qry_end",
                                                                     "ref_start", "ref_end", "eval", "bitscore"))

## fai
rdna_fai <- "/Users/polenpineda/Documents/TulixWagyu/rdna/rdna_35Kb/rdna_rep_chr11_split.fa.fai"
rdna_fai <- read_delim(rdna_fai, col_names = c("qry","rdna_size"))
rdna_fai <- rdna_fai %>% select(qry, rdna_size)

## wagyu fai
wagyu_fai <- "/Users/polenpineda/Documents/TulixWagyu/all_Rscript_in/plot_length/UOA_Wagyu_1.fa.fai"
wagyu_fai <- read_delim(wagyu_fai, col_names = c("ref","wagyu_size"))
wagyu_fai <- wagyu_fai %>% select(ref, wagyu_size)

##
unplaced_fai <- "/Users/polenpineda/Documents/TulixWagyu/reference_genome/Wagyu/each_chr/unplaced.fa.fai"
unplaced_fai_df <- read_delim(unplaced_fai, col_names = c("ref","size"))

## PROCESSS and CLEAN FILE
blast_rdna_df <- blast_rdna_df %>%
  mutate(location = ifelse(grepl("^mat", ref), "unplaced", ref))

blast_rdna_df <- merge(blast_rdna_df, rdna_fai, by = "qry")

blast_rdna_df_filter <- blast_rdna_df %>%
  mutate(align_cover = qry_end - qry_start + 1,
         align_len_per = align_len/rdna_size) %>%
  filter(perc_id > 95, align_len_per > 0.80)

## ANALYSIS

blast_18S <- blast_rdna_df_filter %>% filter(qry == "11_18S")
blast_18S_axis <- 1:1872
counts_18S <- lapply(blast_18S_axis, function(x) {
  subset <- blast_18S[blast_18S$qry_start <= x & blast_18S$qry_end >= x, ]
  count <- nrow(subset)
  perc_id <- mean(subset$perc_id)  # Calculate perc_id based on your actual column
  return(data.frame(x = x, count = count, perc_id = perc_id))
})
blast_18S_axis_per1bp <- do.call(rbind, counts_18S)

blast_28S <- blast_rdna_df_filter %>% filter(qry == "11_28S")
blast_28S_axis <- 1:1751
counts_28S <- lapply(blast_28S_axis, function(x) {
  subset <- blast_28S[blast_28S$qry_start <= x & blast_28S$qry_end >= x, ]
  count <- nrow(subset)
  perc_id <- mean(subset$perc_id)  # Calculate perc_id based on your actual column
  return(data.frame(x = x, count = count, perc_id = perc_id))
})
blast_28S_axis_per1bp <- do.call(rbind, counts_28S)

blast_IGS <- blast_rdna_df_filter %>% filter(qry == "11_IGS")
blast_IGS_axis <- 1:24378
counts_IGS <- lapply(blast_IGS_axis, function(x) {
  subset <- blast_IGS[blast_IGS$qry_start <= x & blast_IGS$qry_end >= x, ]
  count <- nrow(subset)
  perc_id <- mean(subset$perc_id)  # Calculate perc_id based on your actual column
  return(data.frame(x = x, count = count, perc_id = perc_id))
})
blast_IGS_axis_per1bp <- do.call(rbind, counts_IGS)

blast_repeat <- blast_rdna_df_filter %>% filter(qry == "11_repeat_59")
blast_repeat_axis <- 1:51
counts_repeat <- lapply(blast_repeat_axis, function(x) {
  subset <- blast_repeat[blast_repeat$qry_start <= x & blast_repeat$qry_end >= x, ]
  count <- nrow(subset)
  perc_id <- mean(subset$perc_id)  # Calculate perc_id based on your actual column
  return(data.frame(x = x, count = count, perc_id = perc_id))
})
bblast_repeat_axis_per1bp <- do.call(rbind, counts_repeat)

## get coordinates of the rDNA



## unplaced

blast_rdna_df_filter_unplaced <- blast_rdna_df_filter %>%
  filter(location == "unplaced" )

# blast_rdna_df_filter_unplaced_grouped <- blast_rdna_df_filter_unplaced %>%
#   group_by(ref) %>%
#   slice(c(which.min(ref_start), which.max(ref_end))) %>%
#   select(qry, ref, ref_start, ref_end, align_len_per)
# write_tsv(blast_rdna_df_filter_unplaced_grouped, file = paste0(dirout,"/blast_rdna_df_filter_unplaced_grouped.tsv"))

## just get all the 18S

blast_rdna_df_filter_unplaced_18S <- blast_rdna_df_filter_unplaced %>%
  filter(qry == "11_18S") %>%
  arrange(ref, ref_start) %>%
  group_by(ref) %>%
  mutate(rdna_end = lead(ref_start -1),
         rdna_complete_size = abs(ref_start - rdna_end)+1)

rdna_location <- blast_rdna_df_filter_unplaced_18S %>%
  select(ref, start = ref_start, end = rdna_end, rdna_complete_size)
rdna_location <- merge(rdna_location, wagyu_fai, by = "ref")
write_tsv(rdna_location, file = paste0(dirout,"/rdna_location.tsv"))

## count how many complete (~7Kb) rdna are in the unplaced?

blast_rdna_df_filter_unplaced_grouped <- blast_rdna_df_filter_unplaced %>%
  filter(!(qry %in% c("11_repeat_59", "11_ITS1", "11_ITS2", "11_IGS"))) %>%
  group_by(ref, qry) %>%
  summarise(count = n())
blast_rdna_df_filter_unplaced_grouped_filter <- blast_rdna_df_filter_unplaced_grouped %>%
  group_by(ref) %>%
  filter(n() >= 3) %>%
  mutate(final_count = if_else(qry == 28, count / 2, count))

summary_blast_split_df_unplaced_grouped <- blast_rdna_df_filter_unplaced_grouped_filter %>%
  group_by(ref) %>%
  summarise(count = sum(final_count)/3) %>%
  mutate(count = get_first_digit(count))

summary_blast_split_df_unplaced_grouped <- merge(summary_blast_split_df_unplaced_grouped, unplaced_fai_df, by = "ref")
sum(summary_blast_split_df_unplaced_grouped$count)
sum(summary_blast_split_df_unplaced_grouped$size)
# 190 copies
# 7,649,577 total size
# in 100 scaffolds

## count copies of the repeat
blast_repeat_59 <- blast_rdna_df_filter_unplaced %>%
  filter(qry == "11_repeat_59") %>%
  group_by(ref, qry) %>%
  summarise(count = n())

summary_blast_split_df_unplaced_grouped <- merge(summary_blast_split_df_unplaced_grouped, blast_repeat_59, by = "ref", all = TRUE)
summary_blast_split_df_unplaced_grouped <- summary_blast_split_df_unplaced_grouped %>%
  mutate(average_repeat_copy = count.y/count.x)

mean_perc_id <- blast_rdna_df %>%
  group_by(qry) %>%
  summarise(mean_perc_id = mean(perc_id, na.rm = TRUE),
            mean_ave_len = mean(align_len),
            median_len = median(align_len),
            max = max(align_len))
mean_perc_id

## in autosomes
blast_rdna_df_filter_autosomes <- blast_rdna_df_filter %>%
  filter(location != "unplaced" )
blast_rdna_df_filter_filter_grouped <- blast_rdna_df_filter_autosomes %>%
  filter(!(qry %in% c("11_repeat_59", "11_ITS1", "11_ITS2", "11_IGS"))) %>%
  group_by(ref, qry) %>%
  summarise(count = n())
blast_rdna_df_filter_filter_grouped_filter <- blast_rdna_df_filter_filter_grouped %>%
  group_by(ref) %>%
  filter(n() >= 3) %>%
  mutate(final_count = if_else(qry == 28, count / 2, count))

summary_blast_split_df_filter_grouped <- blast_rdna_df_filter_filter_grouped_filter %>%
  group_by(ref) %>%
  summarise(count = sum(final_count)/3) %>%
  mutate(count = get_first_digit(count))

summary_blast_split_df_filter_grouped <- merge(summary_blast_split_df_filter_grouped, wagyu_fai, by = "ref")
sum(summary_blast_split_df_filter_grouped$count)
sum(summary_blast_split_df_filter_grouped$size)
# 9 copies in 4 chromosomes
# chr 2,4,11,25


## PLOT
plot_unplaced_length_id <- ggplot(blast_rdna_df_filter, aes(x = align_len, y = perc_id, color = location)) +
  geom_point() + 
  theme_classic() +
  scale_color_manual(values = c("unplaced" = "orange")) +
  labs(title = "rDNA blast result with alignment length and sequence identity",
       x = "alignment length",
       y = "sequence identity")
plot_unplaced_length_id

plot_blast_18S_axis_per1bp <- ggplot(blast_18S_axis_per1bp, aes(x = x, y = count, color = perc_id)) +
  geom_line(size = 0.5) +
  scale_color_viridis(option = "H", limits = c(85, 100)) +
  theme_classic() +
  labs(title = "rDNA blast result base copy in the repeat (coverage)",
       x = "alignment position",
       y = "base copy number",
       color = "sequence identity")
plot_blast_18S_axis_per1bp

plot_blast_28S_axis_per1bp <- ggplot(blast_28S_axis_per1bp, aes(x = x, y = count, color = perc_id)) +
  geom_line(size = 0.5) +
  scale_color_viridis(option = "H", limits = c(85, 100)) +
  theme_classic() +
  labs(title = "rDNA blast result base copy in the repeat (coverage)",
       x = "alignment position",
       y = "base copy number",
       color = "sequence identity")
plot_blast_28S_axis_per1bp

plot_blast_IGS_axis_per1bp <- ggplot(blast_IGS_axis_per1bp, aes(x = x, y = count, color = perc_id)) +
  geom_line(size = 0.5) +
  scale_color_viridis(option = "H", limits = c(85, 100)) +
  theme_classic() +
  labs(title = "rDNA blast result base copy in the repeat (coverage)",
       x = "alignment position",
       y = "base copy number",
       color = "sequence identity")
plot_blast_IGS_axis_per1bp

plot_bblast_repeat_axis_per1bp <- ggplot(bblast_repeat_axis_per1bp, aes(x = x, y = count, color = perc_id)) +
  geom_line(size = 0.5) +
  scale_color_viridis(option = "H", limits = c(85, 100)) +
  theme_classic() +
  labs(title = "rDNA blast result base copy in the repeat (coverage)",
       x = "alignment position",
       y = "base copy number",
       color = "sequence identity")
plot_bblast_repeat_axis_per1bp

## SRF 51bp repeat
# CACTTGGCCTCCGGAGGGCGACCGAGCCCGGTCGACCAGCCGTCCCGCCGG
## chr 11 repeat
# GGTCGACCAGCCGTCCCGCCGGCACTTGGCCTCCGGAGGGCGGGCCGGCCC


#### get the percent identity per sections

summary_percID <- blast_rdna_df_filter %>%
  group_by(qry) %>%
  summarise(perc_id = mean(perc_id),
            count = n(),
            total_align_len = sum(align_len),
            mean_align_len = mean(align_len),
            sd_align_len = sd(align_len))

## get ETS, IGS coordinates (28S to 18S)
## this contains the rDNA sections information: blast_rdna_df_filter

blast_rdna_df_filter_reorient <- blast_rdna_df_filter %>%
  mutate(tmp_size = ref_end - ref_start,
         strand = ifelse(tmp_size < 0, "-", "+"),
         tmp_start = ifelse(tmp_size < 0, ref_end, ref_start),
         tmp_end = ifelse(tmp_size < 0, ref_start, ref_end),
         size = tmp_end - tmp_start) %>%
  select(-ref_start, -ref_end) %>%
  mutate(ref_start = tmp_start, ref_end = tmp_end)


rdna_ETS_IGS <- blast_rdna_df_filter_reorient %>% filter(qry %in% c("11_18S", "11_28S")) %>%
  group_by(ref) %>%
  arrange(ref, ref_start) %>%
  mutate(group_num = cumsum(qry == "11_28S")) %>%
  ungroup() %>%
  mutate(group_name = paste0(ref,"_",group_num))

rdna_ETS_IGS_filter <- rdna_ETS_IGS %>%
  group_by(group_name) %>%
  summarise(section = paste0(qry, collapse = ","),
            chr = unique(ref),
            start = first(ref_end),
            end = last(ref_start),
            group_name = unique(group_name),
            strand = unique(strand),
            len = end - start) %>%
  mutate(region = "3ETS_IGS_5ETS")

write_tsv(rdna_ETS_IGS_filter, file = paste0(dirout,"/rdna_ETS_IGS_filter.tsv"))



### HERE ONWARD ARE RESULTS FROM BLASTING THE CONSENSUS 35S RDNA (final rDNA consensus)
con_blast_rdna <- "/Users/polenpineda/Documents/TulixWagyu/blast/out/UOA_Wagyu_1.withY_vs_rDNA_35S_split.85.blast.tsv"
con_blast_rdna_df <- read_delim(con_blast_rdna, comment = "#", col_names = c("qry", "ref", "perc_id", "align_len",
                                                                     "mismatches", "gap_opens", "qry_start", "qry_end",
                                                                     "ref_start", "ref_end", "eval", "bitscore"))

## fai
con_rdna_fai <- "/Users/polenpineda/Documents/TulixWagyu/rdna/rdna_35Kb/cattle_consensus_final/rDNA_35S_split.fasta.fai"
con_rdna_fai <- read_delim(con_rdna_fai, col_names = c("qry","rdna_size"))
con_rdna_fai <- con_rdna_fai %>% select(qry, rdna_size)

## PROCESSS and CLEAN FILE
con_blast_rdna_df <- con_blast_rdna_df %>%
  mutate(location = ifelse(grepl("^mat", ref), "unplaced", ref))

con_blast_rdna_df <- merge(con_blast_rdna_df, con_rdna_fai, by = "qry")

con_blast_rdna_df_filter <- con_blast_rdna_df %>%
  mutate(align_cover = qry_end - qry_start + 1,
         align_len_per = align_len/rdna_size) %>%
  filter(perc_id > 95, align_len_per > 0.80)

## ANALYSIS
## 190 in unplaced + 9 in autosomes = 199 copies in total
con_summary_rDNA_region <- con_blast_rdna_df_filter %>%
  group_by(qry) %>%
  summarise(count = n(),
            contig_count = n_distinct(ref),
            perc_id = mean(perc_id),
            total_align = sum(align_len),
            mean_align = mean(align_len),
            mean_cover = mean(align_len_per),
            contigs = n_distinct(ref))

con_ETS_rDNA_region <- con_blast_rdna_df_filter %>%
  filter(qry %in% c("3ETS","5ETS")) %>%
  group_by(ref) %>%
  arrange(ref, ref_start) %>%
  select(qry, ref, perc_id, align_len, ref_start, ref_end) %>%
  mutate(diff = ref_start - lag(ref_end),
         strand = ifelse((ref_end - ref_start)<0, "-", "+"))

con_IGS_coordinates_forward <- con_ETS_rDNA_region %>%
  filter(strand == "+") %>%
  group_by(ref) %>%
  mutate(group_num = cumsum(qry == "3ETS")) %>%
  ungroup()

con_IGS_coordinates_forward_summary <- con_IGS_coordinates_forward %>%
  group_by(ref, group_num) %>%
  summarise(section = paste0(qry, collapse = ","),
            chr = unique(ref),
            start = first(ref_end),
            end = last(ref_start),
            group_num = unique(group_num),
            strand = unique(strand),
            len = end - start) %>%
  filter(len >0, len < 35e3)

con_IGS_coordinates_reverse <- con_ETS_rDNA_region %>%
  filter(strand == "-") %>%
  group_by(ref) %>%
  mutate(group_num = cumsum(qry == "5ETS")) %>%
  ungroup()

con_IGS_coordinates_reverse_summary <- con_IGS_coordinates_reverse %>%
  group_by(ref, group_num) %>%
  summarise(section = paste0(qry, collapse = ","),
            chr = unique(ref),
            start = first(ref_end),
            end = last(ref_start),
            group_num = unique(group_num),
            strand = unique(strand),
            len = end - start) %>%
  filter(len > 20e3, len < 35e3)

con_IGS_coordinates <- rbind(con_IGS_coordinates_forward_summary,con_IGS_coordinates_reverse_summary)
con_IGS_coordinates_bed <- con_IGS_coordinates %>% ungroup %>% select(chr, start, end)

con_IGS_coordinates_forward_summary_bed <- con_IGS_coordinates_forward_summary %>% ungroup %>% select(chr, start, end)
con_IGS_coordinates_reverse_summary_bed <- con_IGS_coordinates_reverse_summary %>% ungroup %>% select(chr, start, end)

write_tsv(con_IGS_coordinates_forward_summary, file = paste0(dirout,"/con_IGS_coordinates_forward_summary.tsv"))
write_tsv(con_IGS_coordinates_reverse_summary, file = paste0(dirout,"/con_IGS_coordinates_reverse_summary.tsv"))
write_tsv(con_IGS_coordinates, file = paste0(dirout,"/con_IGS_coordinates.tsv"))
write_tsv(con_IGS_coordinates_bed, file = paste0(dirout,"/con_IGS_coordinates.bed"), col_names = FALSE)
write_tsv(con_IGS_coordinates_forward_summary_bed, file = paste0(dirout,"/con_IGS_coordinates_forward_summary.bed"), col_names = FALSE)
write_tsv(con_IGS_coordinates_reverse_summary_bed, file = paste0(dirout,"/con_IGS_coordinates_reverse_summary.bed"), col_names = FALSE)

## 51 bp repeat
rdna_order <- c("18S","ITS1","5.8S","ITS2","28S","3ETS","rDNA_51bp","IGS","5ETS")

con_51bp_chr11 <- con_blast_rdna_df_filter %>% filter(ref == "11")
con_51bp_chr11$qry <- factor(con_51bp_chr11$qry, levels = rdna_order)
plot_con_51bp_chr11 <- con_51bp_chr11 %>%
  ggplot(aes(x=ref_start/1e6, y=qry)) +
  geom_point(size = 1, alpha = 0.5) +
  theme_classic() +
  labs(x = "Genomic position (Mb)", y = "rDNA region")
plot_con_51bp_chr11

con_51bp_mat0000056 <- con_blast_rdna_df_filter %>% filter(ref == "mat-0000056")
con_51bp_mat0000056$qry <- factor(con_51bp_mat0000056$qry, levels = rdna_order)
plot_con_51bp_mat0000056 <- con_51bp_mat0000056 %>%
  ggplot(aes(x=ref_start/1e6, y=qry)) +
  geom_point(size = 0.8) +
  theme_classic() +
  labs(x = "Genomic position (Mb)", y = "rDNA region")
plot_con_51bp_mat0000056

## creating heatmap with identity matrix of IGS

library(reshape2)

id_mat_IGS <- "/Users/polenpineda/Documents/TulixWagyu/rdna/rdna_35Kb/extract_IGS_consensus/all/identity_matrix_edited.txt"
id_mat_IGS_df <- read_delim(id_mat_IGS, col_names = TRUE)


data_melt <- melt(id_mat_IGS_df, id.vars = "name")
IGS_heatmap <- data_melt %>%
  arrange(-value)
plot_IGS_heatmap <- IGS_heatmap %>%
  ggplot(aes(y = reorder(name, -value), x = reorder(variable, -value))) +
  geom_tile(aes(fill = value)) +
  scale_fill_viridis("Percentage identity") +
  labs(x = NULL, y = "IGS copies", title = "Sequence identity matrix of the IGS in the rDNA unit", fill = "Sequence identity (%)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
plot_IGS_heatmap
ggsave(filename = paste0(dirout,"/plot_IGS_heatmap.png"), plot = plot_IGS_heatmap, height = 20, width = 20, dpi = 300)

min_val <- min(id_mat_no_col1)
max_val <- max(id_mat_no_col1)
mean_val <- sum(id_mat_no_col1) / (134 * 134)

summary_id_mat_IGS_df <- data.frame(
  Min = min_val,
  Max = max_val,
  Mean = mean_val
)

## 35Kb

id_mat_35 <- "/Users/polenpineda/Documents/TulixWagyu/rdna/rdna_35Kb/extract_30Kb/extract_reorient/identity_matrix_edited.txt"
id_mat_35_df <- read_delim(id_mat_35, col_names = TRUE)


data_melt_35 <- melt(id_mat_35_df, id.vars = "name")
rdna_heatmap <- data_melt_35 %>%
  arrange(-value)
plot_rdna_heatmap <- rdna_heatmap %>%
  ggplot(aes(y = reorder(name, -value), x = reorder(variable, -value))) +
  geom_tile(aes(fill = value)) +
  scale_fill_viridis("Percentage identity") +
  labs(x = NULL, y = "IGS copies", title = "Sequence identity matrix of the 35Kb rDNA unit", fill = "Sequence identity (%)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 40),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
plot_rdna_heatmap
ggsave(filename = paste0(dirout,"/plot_rdna_heatmap.png"), plot = plot_rdna_heatmap, height = 20, width = 20, dpi = 300)
