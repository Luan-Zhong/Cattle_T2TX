## 2024.09.02
## Repeat analysis with repeatmasker output
## The library used for repeatmasker were a custom library:
#### 1. combine the default repeatmasker library with RepBase RepeatMasker Edition - version 20181026
#### 2. extract the fasta file for bos taurus specific reapets using famdb.py "Bos taurus"
#### 3. replace the "Satellite/centromeric" repeats with the custom satellites (SATI-SATVII)
#### 4. run repeatmasker with the custom library in: /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/tools/repeatmasker/RepeatMasker/custom_lib/RepeatMaskerLib.h5.bostaurus-rm.withsats.fa
#### 5. run RM2bed.py in repeatmasker utils as postprocessing to remove overlaps with maximum 40 divergence
#### 6. concatenate the out.bed files

library(dplyr)
library(readr)
library(ggplot2)
library(scales)
library(tidyverse)
library(viridis)
library(hrbrthemes)
library(GenomicRanges)

####### FUNCTIONS #######

fill_gaps <- function(df) {
  df <- df %>%
    arrange(chr, start)

  new_rows <- data.frame(chr=character(), start=integer(), end.x=integer(),
                         family=character(), score=character(), orientation=character(),
                         class=character(), subclass=character(),
                         divergence = integer(), linkage_id=character(), align_len=integer(),
                         class_family=character(), end.y=character(), size=character(),
                         copies_aligned = character(),
                         stringsAsFactors=FALSE)

  for (i in 1:(nrow(df) - 1)) {
    current_end <- df$end.x[i]
    next_start <- df$start[i + 1]

    if (next_start > current_end + 1) {
      new_row <- data.frame(
        chr = df$chr[i],
        start = current_end + 1,
        end.x = next_start - 1,
        family = "Others",
        score = "",
        orientation = "+",
        class = "Others",
        subclass = "Others",
        divergence = -1,
        linkage_id = "",
        align_len = (next_start - 1) - (current_end + 1),
        class_family = "Others",
        end.y = "",
        size = "",
        copies_aligned = "",
        stringsAsFactors = FALSE
      )
      new_rows <- rbind(new_rows, new_row)
    }
  }

  df <- rbind(df, new_rows)
  df <- df %>%
    arrange(chr, start)

  return(df)
}

####### INPUT FILES #######
### INSERT INPUT FILES HERE
dirout <- "/Users/polenpineda/Documents/TulixWagyu/all_Rscript_out/final_repeats_2/"
#repeatmasker
rm <- "/Users/polenpineda/Documents/TulixWagyu/repeat_masker/concat_rm.bed"
rm <- read_delim(rm, col_names = c("chr","start","end","family","score","orientation","class","subclass","divergence","linkage_id"))

rm_tsv <- "/Users/polenpineda/Documents/TulixWagyu/repeat_masker/concat_rm.out.tsv"
rm_tsv <- read_delim(rm_tsv, col_names = c("SW_score","perc_div","perc_del","perc_ins",
                                           "chr","qry_begin","qry_end","qry_left","orientation","query_repeat","family",
                                           "rep_begin","rep_end","rep_left","ID","remarks"))

## divergence is -1.0 because it only works for .aln files
## rm is the bed file (no overlap and is already filtered to <40% divergence)
## rm_tsv is the .out file which has more details than the bed file (contains the perc id)

#cenpa broad peaks
cenpa_broad <- "/Users/polenpineda/Documents/TulixWagyu/CENPA/cenpa_level/from_Callum/LIB203001.stringent.bedgraph.stringent.col4.bedgraph"
cenpa_broad <- read_delim(cenpa_broad, col_names = c("chr","start","end","total_cenpa"))
head(cenpa_broad)

#cenpa max peaks
cenpa_max <- "/Users/polenpineda/Documents/TulixWagyu/CENPA/cenpa_level/from_Callum/LIB203001.stringent.bedgraph.stringent.bed"
cenpa_max <- read_delim(cenpa_max, col_names = c("chr","start","end","total_cenpa","peak_cenpa","peak_region"))
head(cenpa_max)

#cenpa1bp
cenpa_perbase <- "/Users/polenpineda/Documents/TulixWagyu/CENPA_repeatmasker/from_repeatmasker_1base_coverage/concat_CENPA_perbase.tsv"
cenpa_perbase <- read_delim(cenpa_perbase, col_names = c("chr", "start", "end", "query_repeat", "align_len", "strand", "class", "subclass", "divergence", "index", "cenpa_level"))

#satellites
sat_fai <- "/Users/polenpineda/Documents/TulixWagyu/all_Rscript_in/centromere_boundary/final_sats.fa.fai"
sat_len_df <- read_delim(sat_fai, col_names = c("family", "size", "start", "base_line", "bytes_line"))
sat_len_df <- sat_len_df %>% filter(family != "telomere") %>% select(family,size) %>% mutate(query_repeat = family)

#wagyu genome
wagyu_fai <- "/Users/polenpineda/Documents/TulixWagyu/all_Rscript_in/plot_length/UOA_Wagyu_1.fa.fai"
wagyu_fai_df <- read_delim(wagyu_fai, col_names = c("chr", "size", "start", "base_line", "bytes_line"))
wagyu_fai_df <- wagyu_fai_df %>% select(chr,size)
wagyu_genome_size_chr <- wagyu_fai_df %>%
  filter(chr %in% c(1:29, "X")) %>%
  summarise(size = sum(size))
# 3031482482
unplaced_genome_size <- wagyu_fai_df %>%
  filter(!chr %in% c(1:29, "X")) %>%
  summarise(size = sum(size))
# 110927292

#satellite color
sat_color <- data.frame(family = c("SATI", "SATII", "SATIII", "SATVI", "SATV", "SATVI", "SATVII"),
                        color = c("#dd314b", "#00008b", "darkgreen", "#f4895f", "limegreen", "#7D0DC3", "#CCCC00"))

####### file cleaning #######
# check if there's overlap filtered? ##none
check_ovlp <- rm %>%
  arrange(chr, start) %>%
  mutate(diff = start - lag(end)) %>%
  filter(diff < 0)

## set up and clean the input files

### cleaning the columns in the repeatmasker output file
rm_tsv$new_rep_begin <- as.integer(gsub("[()]", "", rm_tsv$rep_begin))
rm_tsv$new_rep_end <- as.integer(gsub("[()]", "", rm_tsv$rep_end))
rm_tsv$new_rep_left <- as.integer(gsub("[()]", "", rm_tsv$rep_left))
rm_tsv$new_qryleft <- as.integer(gsub("[()]", "", rm_tsv$qry_left))

rm_tsv_raw <- rm_tsv

## set up the file for the repeat bed
rm_df <- rm %>%
  mutate(align_len = end - start +1,
         class_family = paste0(class, "/", subclass),
         repeat_type = ifelse(class == "Satellite", family, class))
rm_df <- merge(rm_df, sat_len_df, by = "family", all = TRUE)
rm_df <- rm_df %>% mutate(cover_align = ifelse(is.na(align_len / size), 1, align_len / size)) %>% filter(cover_align >= 0.80) %>%
  select(-query_repeat)
rm_df_auto <- rm_df %>% filter(chr %in% 1:29)
rm_df_auto$chr <- as.numeric(rm_df_auto$chr)
rm_df_unplaced <- rm_df %>% filter(!chr %in% c(1:29, "X"))
rm_df_auto_sats <- rm_df_auto %>% filter(class == "Satellite")

## get copy numbers and sequence identity (I will use the .out file from repeatmasker)
rm_tsv <- merge(rm_tsv, sat_len_df, by = "query_repeat", all = TRUE)
rm_tsv <- rm_tsv %>% mutate(
  align_len = qry_end-qry_begin,
  family = family.x,
  cover_align = ifelse(is.na(align_len / size), 1, align_len / size)) %>% filter(cover_align >= 0.80) %>%
  select(-family.y)
rm_tsv <- rm_tsv %>%
  mutate(
    align_len = qry_end - qry_begin + 1,
    perc_id = 1 - (perc_div / 100),
    tmp = family
  ) %>%
  filter(perc_id > 0.6) %>%
  separate(tmp, into = c("class", "subclass"), sep = "/", remove = FALSE) %>%
  mutate(subclass = coalesce(subclass, class),
         repeat_type = ifelse(class == "Satellite", query_repeat, class),
         category = ifelse(
           class == "Satellite", query_repeat,
           ifelse(class %in% c("LINE", "SINE", "LTR"), class, "Others")))
unique(rm_tsv$repeat_type)


rm_tsv_auto <- rm_tsv %>% filter(chr %in% 1:29) %>% mutate(chr = as.numeric(chr))
rm_tsv_unplaced <- rm_tsv %>% filter(!chr %in% c(1:29, "X"))
rm_tsv_auto_sats <- rm_tsv_auto %>% filter(family == "Satellite/centr") %>% mutate(chr = as.numeric(chr))

### ANALYSIS ###

####### CTR BOUNDARY #######
#### get the centromere boundary
### use l-apply to all the autosomes to get the centromere boundaries
chrom_list <- c(1:29)

group_func <- function(chrom) {
  df <- rm_df_auto_sats %>%
    filter(chr == chrom) %>%
    arrange(start) %>%
    mutate(diff = start - lag(end),
           cover = end-start+1)

  df$group <- cumsum(!is.na(df$diff) & df$diff > 100000)

  df_summary <- df %>%
    group_by(group) %>%
    summarise(chromosome = chr,
              start = min(start),
              end = max(end),
              rep_size = sum(cover),
              loc_size = max(end) - min(start) + 1,
              asm = "wagyu")
  assign(paste0("zctr_raw_loc_", chrom), df, envir = .GlobalEnv)
  assign(paste0("zctr_loc_", chrom), df_summary, envir = .GlobalEnv)
}

invisible(lapply(chrom_list, group_func))

## summarise all the repeat-rich regions and only get the most repeat blocks
ctr_all_repclusters <- do.call(rbind, lapply(chrom_list, function(chrom) get(paste("zctr_loc_", chrom, sep = ""))))
gg <- ggplot(ctr_all_repclusters, aes(x = start/1e6, y = loc_size / 1e6, color = chromosome)) +
  geom_line(size = 0.3) +
  labs(title = "Repeat cover size per chromosomes", x = "Position (Mb)", y = "Repeat cover (Mb)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)
gg

ctr_all_repclusters$diff <- c(NA, abs(diff(ctr_all_repclusters$rep_size)))
ctr_all_repclusters_df_filter <- ctr_all_repclusters %>%
  group_by(chromosome) %>%
  filter(diff > 20000 | is.na(diff)) %>%
  filter(row_number() < n())
mean(ctr_all_repclusters_df_filter$rep_size)

ctr_all_repclusters_df_filter_final <- ctr_all_repclusters_df_filter %>%
  group_by(chromosome) %>%
  summarise(start = min(start), end = max(end)) %>%
  mutate(chromosome = as.numeric(chromosome)) %>%
  arrange(chromosome)
sum(ctr_all_repclusters_df_filter_final$end)

####### ctr loc #######
ctr_loc <- ctr_all_repclusters_df_filter_final %>% select(chr = chromosome, end)

## checking the centromere locations if they are indeed within the repeatclusters
## check chr14, 25, 27

rm_df_25 <- rm_df %>% filter(chr == 25)
rm_df_27 <- rm_df %>% filter(chr == 27)
rm_df_23 <- rm_df %>% filter(chr == 23, end < 14104924)

plot_ctrcheck_chr <- zctr_raw_loc_27 %>%
  ggplot(aes(x = start/1e6, y = repeat_type, color = repeat_type)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c(
    "SATVII" = "#CCCC00",
    "SATII" = "#00008b",
    "SATIV" = "#f4895f",
    "SATVI" = "#7D0DC3",
    "SATI" = "#dd314b",
    "SATV" = "limegreen",
    "SATIII" = "darkgreen")) +
  geom_hline(data = ctr_all_repclusters_df_filter_final %>%
               filter(chromosome == unique(zctr_raw_loc_27$chr)),
             aes(yintercept = end / 1e6),
             linetype = "dashed", color = "red")
plot_ctrcheck_chr

## end of centromere boundary analysis

####### Prepare files with ctr #######
rm_tsv_auto <- left_join(rm_tsv_auto,ctr_loc, by = "chr")
rm_tsv_auto <- rm_tsv_auto %>% mutate(copies_aligned = if_else(is.na(size), 1, round(align_len / size)))

rm_tsv_auto_sats <- rm_tsv_auto %>% filter(family == "Satellite/centr") %>% mutate(chr = as.numeric(chr)) %>%
  mutate(SAT_category = case_when(
    repeat_type == "SATI" & perc_id >= 0.95 ~ "SAT1a",
    repeat_type == "SATI" & perc_id < 0.95 & perc_id > 0.90 ~ "SAT1b",
    repeat_type == "SATI" & perc_id <= 0.90 ~ "SAT1c",
    repeat_type == "SATII" & perc_id >= 0.95 ~ "SATIIa",
    repeat_type == "SATII" & perc_id < 0.95 ~ "SATIIb",
    repeat_type == "SATVI" & perc_id >= 0.98 ~ "SATVIa",
    repeat_type == "SATVI" & perc_id < 0.98 ~ "SATVIb",
    TRUE ~ repeat_type
  ))
rm_tsv_auto_sats_ctr <- rm_tsv_auto_sats %>% filter(qry_end <= end)

rm_df_auto <- left_join(rm_df_auto, ctr_loc, by = "chr")
rm_df_auto <- rm_df_auto %>% mutate(copies_aligned = if_else(is.na(size), 1, round(align_len / size)))

## the rm_df_auto do not have overlaps while rm_tsv_auto have overlaps but contains divergence
## get the divergence of the repeat satellites from rm_tsv_auto to rm_df_auto
### 4556793-4569165 = 12372

## input files pre-processing done

## START OF ANALYSIS

### FINAL FILES TO USE IN SUBSEQUENT ANALYSIS
### rm_df       for all of the chromosomes and unplaced (however, this does not contain divergence identity)
### columns:  class         15 types (general repeat family names)
###           sublass       45 types (what types are the general repeat family names)
###           repeat_type   21 types (contains the class name but the satellite are satellite repeats (eg. SATI, SATII))
###           family        12852 types (the specific repeat types)
###           class_family  55 types (has the class and subclass)

## save the final bedfile of the repeatmasker satellite repeats

####### SAVING RM CLEAN FILES #######
rm_df_auto_sats <- rm_tsv_auto %>% filter(family == "Satellite/centr") %>%
  mutate(copies_group = cut(copies_aligned,
                            breaks = c(0, 100, 200, 300, 400, 500, 1000, Inf),
                            labels = c("1-100", "101-200", "201-300", "301-400", "401-500", "501-1000", "1000+")))

rm_df_sats <- rm_df %>% filter(class_family == "Satellite/centr")
rm_df_sats_ctr <- merge(rm_df_sats, ctr_loc, by = "chr")
rm_df_sats_ctr <- rm_df_sats_ctr %>% filter(end.x <= end.y ) %>% select(-end.y, end = end.x, -subclass, -linkage_id, -divergence, -family)


write_tsv(rm_df_sats_ctr, file = paste0(dirout, "repeatmasker_clean_bovsats_only_in_ctr.tsv"), col_names = TRUE)
write_tsv(rm_df, file = paste0(dirout, "repeatmasker_clean.tsv"))
write_tsv(rm_tsv_auto, file = paste0(dirout, "repeatmasker_clean_with_divergence_but_has_overlaps.tsv"))

## T2T chromosomes

T2T_ctr <- rm_df%>% filter(chr %in% c(9, 10, 21, 23))

####### REPEATS SUMMARY #######
## SUMMARY OF THE REPEATS ENTIRE GENOME WITH UNPLACED
size_genome <- sum(wagyu_fai_df$size) #3142409774 (entire genome)
size_autosomes <- sum(wagyu_fai_df$size[wagyu_fai_df$chr %in% 1:29], na.rm = TRUE) #2857872688 (autosomes)
size_X <- sum(wagyu_fai_df$size[wagyu_fai_df$chr %in% "X"], na.rm = TRUE) #173609794 (X)
size_autosomes_X <- sum(wagyu_fai_df$size[wagyu_fai_df$chr %in% c(1:29,"X")], na.rm = TRUE) #3031482482 (autosomes + X)
size_unplaced <- sum(wagyu_fai_df$size[!(wagyu_fai_df$chr %in% c(1:29, "X"))], na.rm = TRUE) #110927292 (unplaced)
#3142409774-1587976420

## repeats in the entire genome
### 50.53372% repeats in the genome
### 49.58363% repeats in the chromosome
### 76.49811% in the unplaced are repeats
summary_repeats_genome_entire <- rm_df %>%
  summarise(total_repeats = sum(align_len),
            total_repeats_genome = total_repeats/size_genome*100)
print(summary_repeats_genome_entire)
summary_repeats_genome_chr_and_X <- rm_df %>%
  filter(chr %in% c(1:29,"X")) %>%
  summarise(total_repeats = sum(align_len),
            total_repeats_genome = total_repeats/size_autosomes_X*100)
print(summary_repeats_genome_chr_and_X)
summary_repeats_genome_unplaced <- rm_df_unplaced %>%
  summarise(total_repeats = sum(align_len),
            total_repeats_genome = total_repeats/size_genome*100,
            total_repeats_unplaced = total_repeats/size_unplaced*100,
            total_repeats_unplaced_satellites = sum(align_len[class == "Satellite"], na.rm = TRUE),
            perc_repeats_unplaced_satellites = total_repeats_unplaced_satellites/size_unplaced*100,
            perc_repeats_unplaced_others = sum(align_len[class != "Satellite"], na.rm = TRUE)-(332488+7649577))
print(summary_repeats_genome_unplaced)
#84857278/110927292 - 76.49811% are repeats
#76705124/110927292 - 69.14901% are Satellites
#7649577/110927292 - 6.896028% are rDNA
#332488/110927292 - 0.2997351% are telomere region
#170089/110927292 - 0.1533338% are other repeats
## 69.14901+6.896028+0.2997351+0.1533338 = 76.34477
#26291536/110927292 - 23% are telomere contigs (this is the size of the contigs that contains telomere)

## summary of the family and class
summary_repeats_genome_entire_class <- rm_df %>%
  group_by(class) %>%
  summarise(total_repeats = sum(align_len),
            total_repeats_genome = total_repeats/size_genome*100)

summary_repeats_genome_autosomes <- rm_df %>%
  filter(chr %in% c(1:29)) %>%
  group_by(class) %>%
  summarise(total_repeats = sum(align_len),
            total_repeats_genome = total_repeats/(size_autosomes)*100)

summary_repeats_genome_entire_class_2 <- summary_repeats_genome_entire_class %>%
  select(class, total_repeats) %>%
  mutate(class = ifelse(class %in% c('DNA', 'Unknown', 'RC', 'snRNA', 'tRNA', 'srpRNA', 'scRNA', 'ARTEFACT'), 'Others', class)) %>%
  group_by(class) %>%
  summarise(total_repeats = sum(total_repeats, na.rm = TRUE)) %>%
  ungroup() %>%
  bind_rows(tibble(class = "Non-repetitive", total_repeats = summary_repeats_genome_entire$total_repeats)) %>%
  mutate(total_repeats_genome = total_repeats/size_genome*100)


plot_summary_repeats_genome_entire_class <- summary_repeats_genome_entire_class_2 %>%
  ggplot(aes(x = total_repeats_genome, y = reorder(class, total_repeats_genome), fill = class)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "Percentage in the genome", y = "") +
  theme(legend.position = "none") +
  geom_text(aes(label = paste0(round(total_repeats_genome, 2), "%")), hjust = -0.2, color = "black") +
  xlim(0, 60)
plot_summary_repeats_genome_entire_class
ggsave(filename = paste0(dirout, "plot_summary_repeats_genome_entire_class.png"), plot = plot_summary_repeats_genome_entire_class, width = 5, height = 2, dpi = 300)


####### SAT COPIES #######
### count all SAT copies within the genome
summary_rm_tsv_auto_sats_perchr <- rm_tsv_auto_sats %>%
  group_by(chr,SAT_category) %>%
  reframe(copy = as.numeric(sum(copies_aligned)),
            total_len = as.numeric(sum(align_len)),
            mean_percID = as.numeric(mean(perc_id)),
            min_percID = as.numeric(min(perc_id)),
            max_percID = as.numeric(max(perc_id)),
            sd_percID = as.numeric(sd(perc_id)),
            median_perdID = as.numeric(median(perc_id)),
            IQR_percID = as.numeric(IQR(perc_id)))
write_tsv(summary_rm_tsv_auto_sats_perchr, file = paste0(dirout, "summary_bovsats_auto_percID_perchr.txt"))

plot_check_SATI_groups <- rm_tsv_auto_sats %>%
  filter(repeat_type %in% c("SATI")) %>%
  ggplot(aes(x = qry_begin/1e6, y = perc_id, color = SAT_category)) +
  geom_point() + theme_classic() +
  facet_wrap(~ chr, scales = "free_x")
plot_check_SATI_groups

plot_check_SATII_groups <- rm_tsv_auto_sats %>%
  filter(repeat_type %in% c("SATII")) %>%
  ggplot(aes(x = qry_begin/1e6, y = perc_id, color = SAT_category)) +
  geom_point() + theme_classic() +
  facet_wrap(~ chr, scales = "free_x")
plot_check_SATII_groups

plot_check_SATVI_groups <- rm_tsv_auto_sats %>%
  filter(repeat_type %in% c("SATVI")) %>%
  ggplot(aes(x = qry_begin/1e6, y = perc_id, color = SAT_category)) +
  geom_point() + theme_classic() +
  facet_wrap(~ chr, scales = "free_x")
plot_check_SATVI_groups


### count all SAT copies within the centromere

summary_rm_tsv_auto_sats_ctr_perchr <- rm_tsv_auto_sats_ctr %>%
  group_by(chr,SAT_category) %>%
  reframe(copy = as.numeric(sum(copies_aligned)),
          total_len = as.numeric(sum(align_len)),
          mean_percID = as.numeric(mean(perc_id)),
          min_percID = as.numeric(min(perc_id)),
          max_percID = as.numeric(max(perc_id)),
          sd_percID = as.numeric(sd(perc_id)),
          median_perdID = as.numeric(median(perc_id)),
          IQR_percID = as.numeric(IQR(perc_id)))
write_tsv(summary_rm_tsv_auto_sats_ctr_perchr, file = paste0(dirout, "summary_bovsats_auto_percID_perchr_within_ctr.txt"))

summary_rm_tsv_auto_sats_ctr_and_genome_perchr <- merge(summary_rm_tsv_auto_sats_ctr_perchr, summary_rm_tsv_auto_sats_perchr, by = c("chr","SAT_category"), all = TRUE)
write_tsv(summary_rm_tsv_auto_sats_ctr_and_genome_perchr, file = paste0(dirout, "summary_bovsats_auto_percID_perchr_within_ctr_and_genome.txt"))

####### CENP-A analysis #######

####### CENP-A file prep #######
rm_tsv_filter_cenpa <- rm_df_auto_sats %>% filter(family == "Satellite/centr",
                                                copies_aligned > 0)
cenpa_sats <- rm_tsv_filter_cenpa %>%
  select(chr, query_repeat, start = qry_begin, end = qry_end, align_len, orientation, perc_id, ctr_end = end, size, copies_aligned)

cenpa_broad <- cenpa_broad %>%
  mutate(cenpa_align = end - start +1,
         cenpa_per_bases = total_cenpa/cenpa_align)

cenpa_max <- cenpa_max %>%
  separate(peak_region, into = c("chr","maxpeak_start","maxpeak_end"), sep = "[:-]") %>%
  mutate(start = as.numeric(maxpeak_start), end = as.numeric(maxpeak_end)) %>%
  select(chr, start, end, peak_cenpa) %>%
  mutate(cenpa_align = end - start +1)

cenpa_max_chr <- cenpa_max %>%
  filter(chr %in% c(1:29,"X"))

####### CENP-A intersect #######
## get the cenpa broad peaks per satellite repeats
## warnings are from the unplaced
cenpa_sats_gr1 <- GRanges(seqnames = cenpa_sats$chr,
                          ranges = IRanges(start = cenpa_sats$start, end = cenpa_sats$end))
cenpa_broad_gr2 <- GRanges(seqnames = cenpa_broad$chr,
                           ranges = IRanges(start = cenpa_broad$start, end = cenpa_broad$end))
cenpa_broad_overlaps <- findOverlaps(cenpa_sats_gr1, cenpa_broad_gr2)
cenpa_broad_intersected_ranges <- pintersect(cenpa_sats_gr1[queryHits(cenpa_broad_overlaps)],
                                             cenpa_broad_gr2[subjectHits(cenpa_broad_overlaps)])
cenpa_broad_intersected_df <- as.data.frame(cenpa_broad_intersected_ranges)

cenpa_sats_gr1_indices <- data.frame(cenpa_sats_gr1_index = seq_along(cenpa_sats_gr1))
cenpa_broad_gr2_indices <- data.frame(cenpa_broad_gr2_index = seq_along(cenpa_broad_gr2))
cenpa_sats_with_index <- cbind(cenpa_sats, cenpa_sats_gr1_index = seq_len(nrow(cenpa_sats)))
cenpa_broad_with_index <- cbind(cenpa_broad, cenpa_broad_gr2_index = seq_len(nrow(cenpa_broad)))

cenpa_broad_intersected_df_merged <- cenpa_broad_intersected_df %>%
  mutate(cenpa_sats_gr1_index = queryHits(cenpa_broad_overlaps), cenpa_broad_gr2_index = subjectHits(cenpa_broad_overlaps)) %>%
  left_join(cenpa_sats_with_index, by = "cenpa_sats_gr1_index") %>%
  left_join(cenpa_broad_with_index, by = "cenpa_broad_gr2_index")

### get the cenpa max peaks per satellite repeats
## warnings are from the unplaced

cenpa_max_gr2 <- GRanges(seqnames = cenpa_max_chr$chr,
                         ranges = IRanges(start = cenpa_max_chr$start, end = cenpa_max_chr$end))
cenpa_max_overlaps <- findOverlaps(cenpa_sats_gr1, cenpa_max_gr2)
cenpa_max_intersected_ranges <- pintersect(cenpa_sats_gr1[queryHits(cenpa_max_overlaps)],
                                           cenpa_max_gr2[subjectHits(cenpa_max_overlaps)])
cenpa_max_intersected_df <- as.data.frame(cenpa_max_intersected_ranges)

cenpa_max_gr2_indices <- data.frame(cenpa_max_gr2_index = seq_along(cenpa_max_gr2))
cenpa_max_with_index <- cbind(cenpa_max_chr, cenpa_max_gr2_index = seq_len(nrow(cenpa_max_chr)))
cenpa_max_intersected_df_merged <- cenpa_max_intersected_df %>%
  mutate(cenpa_sats_gr1_index = queryHits(cenpa_max_overlaps), cenpa_max_gr2_index = subjectHits(cenpa_max_overlaps)) %>%
  left_join(cenpa_sats_with_index, by = "cenpa_sats_gr1_index") %>%
  left_join(cenpa_max_with_index, by = "cenpa_max_gr2_index")

####### Kruskal-Wallis test #######
## analyse the data
## cenpa_broad_intersected_df_merged (this is where the broad peaks are)
## cenpa_max_intersected_df_merged (this is where the max peaks are)
### 1. What is the average CENPA level per satellite repeats?
### 2. Does CENPA have correlation to the sequence identity?
### 3. Kruskal-Wallis rank test to see statistical difference of the CENPA with the satellite repeats

### Process input file: filter when the intersection is not 50% of the satellite repeat size

cenpa_broad_final_raw <- cenpa_broad_intersected_df_merged %>%
  mutate(align_cover = width/size,
         cenpa_by_width = cenpa_per_bases * width) %>%
  filter(align_cover > 0.5, as.numeric(seqnames) < 30)

# normalised_average_CENP_A_level = total cenpa level that intersects with the satellite repeats
cenpa_broad_final <- cenpa_broad_final_raw %>%
  select(seqnames, start.x, end.x, width, query_repeat, start.y, end.y, ctr_end, size, total_cenpa, cenpa_per_bases, cenpa_by_width) %>%
  mutate(copy_row = ceiling(width/size)) %>%
  uncount(copy_row) %>%
  group_by(start.x, end.x, total_cenpa, size) %>%
  mutate(cenpa_per_sat = cenpa_by_width / n()) %>%
  ungroup

cenpa_max_filter <- cenpa_max_intersected_df_merged %>%
  mutate(align_cover = width/size,
         cenpa_per_bases = peak_cenpa/cenpa_align,
         normalised_average_CENP_A_level = cenpa_per_bases * width) %>%
  filter(align_cover > 0.5, as.numeric(seqnames) < 30)

### 1. What is the average CENPA level per satellite repeats?
### 3. Kruskal-Wallis rank test to see statistical difference of the CENPA with the satellite repeats
## Kruskal–Wallis test
# KW_broad_df <- cenpa_broad_final %>%
#   group_by(query_repeat) %>%
#   summarise(
#     count = n(),
#     total = sum(normalised_average_CENP_A_level, na.rm = TRUE),
#     mean = mean(normalised_average_CENP_A_level, na.rm = TRUE),
#     sd = sd(normalised_average_CENP_A_level, na.rm = TRUE),
#     median = median(normalised_average_CENP_A_level, na.rm = TRUE),
#     IQR = IQR(normalised_average_CENP_A_level, na.rm = TRUE)
#   )

cenpa_max_X_ctr <- cenpa_max %>% filter(chr == "X", start >= 38e6, end <= 50e6) %>% mutate(query_repeat = "X_ctr")
cenpa_max_X_53 <- cenpa_max %>% filter(chr == "X", start >= 110.8e6, end <= 111.3e6) %>% mutate(query_repeat = "X_53bp")

cenpa_max_final <- cenpa_max_filter %>%
  select(chr = seqnames, start = start.x, end = end.x, peak_cenpa, cenpa_align, query_repeat)

cenpa_max_all <- rbind(cenpa_max_final, cenpa_max_X_ctr) # , cenpa_max_X_53)

KW_broad_df <- cenpa_broad_final %>%
  group_by(query_repeat) %>%
  summarise(
    count = n(),
    size = sum(width),
    total = sum(cenpa_per_sat, na.rm = TRUE),
    mean = mean(cenpa_per_bases, na.rm = TRUE),
    sd = sd(cenpa_per_bases, na.rm = TRUE),
    median = median(cenpa_per_bases, na.rm = TRUE),
    IQR = IQR(cenpa_per_bases, na.rm = TRUE))
print(KW_broad_df)




kruskal.test(cenpa_per_bases ~ query_repeat, data = cenpa_broad_final)
# data:  cenpa_per_bases by query_repeat
# Kruskal-Wallis chi-squared = 61.976, df = 2, p-value = 3.484e-14

pairwise.wilcox.test(cenpa_broad_final$cenpa_per_bases, cenpa_broad_final$query_repeat,
                     p.adjust.method = "fdr")
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction
#
# data:  cenpa_broad_final$normalised_average_CENP_A_level and cenpa_broad_final$query_repeat
#
# SATIII  SATV
# SATV   < 2e-16 -
#   SATVII 1.1e-10 0.0046
#
# P value adjustment method: fdr

KW_max_df_peak <- cenpa_max_all %>%
  group_by(query_repeat) %>%
  summarise(
    count = n(),
    total = sum(peak_cenpa, na.rm = TRUE),
    mean = mean(peak_cenpa, na.rm = TRUE),
    sd = sd(peak_cenpa, na.rm = TRUE),
    median = median(peak_cenpa, na.rm = TRUE),
    IQR = IQR(peak_cenpa, na.rm = TRUE)
  )

# KW_max_df_level <- cenpa_max_final %>%
#   group_by(query_repeat) %>%
#   summarise(
#     count = n(),
#     total = sum(peak_cenpa, na.rm = TRUE),
#     mean = mean(peak_cenpa, na.rm = TRUE),
#     sd = sd(peak_cenpa, na.rm = TRUE),
#     median = median(peak_cenpa, na.rm = TRUE),
#     IQR = IQR(peak_cenpa, na.rm = TRUE),
#     mean_all = sum(align_len)/sum(peak_cenpa)
#   )

kruskal.test(peak_cenpa ~ query_repeat, data = cenpa_max_all)

pairwise.wilcox.test(cenpa_max_all$peak_cenpa, cenpa_max_all$query_repeat,
                     p.adjust.method = "fdr")

### 2. Does CENPA have correlation to the sequence identity?

plot_box_CENPA_broad_perbase <- cenpa_broad_final %>%
  ggplot(aes(x = query_repeat, y = cenpa_per_bases, color = query_repeat)) +
  geom_boxplot() +
  labs(
    x = "",
    y = "Normalised CENP-A level"
  ) +
  theme_classic() +
  scale_color_manual(values = c(
    "SATVII" = "#CCCC00",
    "SATV" = "limegreen",
    "SATIII" = "darkgreen"), guide = "none")
plot_box_CENPA_broad_perbase

plot_box_CENPA_broad_total <- cenpa_broad_final %>%
  ggplot(aes(x = query_repeat, y = cenpa_per_sat, color = query_repeat)) +
  geom_boxplot() +
  labs(
    x = "",
    y = "CENP-A level"
  ) +
  theme_classic() +
  scale_color_manual(values = c(
    "SATVII" = "#CCCC00",
    "SATV" = "limegreen",
    "SATIII" = "darkgreen"), guide = "none") +
  scale_y_continuous(labels = comma)
plot_box_CENPA_broad_total

ggsave(filename = paste0(dirout, "plot_box_CENPA_broad_total.png"), plot = plot_box_CENPA_broad_total, width = 5, height = 3, dpi = 300)
ggsave(filename = paste0(dirout, "plot_box_CENPA_broad_perbase.png"), plot = plot_box_CENPA_broad_perbase, width = 7, height = 4, dpi = 300)

## SATV is clustering between less than .8 and greater than 0.8

# cenpa_broad_SATV_less80 <- cenpa_broad_final %>%
#   filter(perc_id < .8)
# cenpa_broad_SATV_more80 <- cenpa_broad_final %>%
#   filter(perc_id >= .8)

# summary_cenpa_broad <- cenpa_broad_final %>%
#   group_by(query_repeat) %>%
#   summarise(
#     count = n(),  # Number of rows in each group
#     mean_value = mean(normalised_average_CENP_A_level, na.rm = TRUE),
#     sd_value = sd(normalised_average_CENP_A_level, na.rm = TRUE),
#     min_value = min(normalised_average_CENP_A_level, na.rm = TRUE),
#     max_value = max(normalised_average_CENP_A_level, na.rm = TRUE),
#     median_value = median(normalised_average_CENP_A_level, na.rm = TRUE),
#     IQR_value = IQR(normalised_average_CENP_A_level, na.rm = TRUE)
#   )

####### CENP-A per base #######
cenpa_perbase_sats <- cenpa_perbase %>% filter(class == "Satellite")
cenpa_perbase_sats <- merge(cenpa_perbase_sats, sat_len_df, by = "query_repeat")
cenpa_perbase_sats <- cenpa_perbase_sats %>%
  mutate(copies_align = align_len/size) %>%
  filter(copies_align > 0.8)

#### can't do per base because I now have HOR

####### other CENP-A peaks #######

cenpa_broad_X <- cenpa_broad %>% filter(chr == "X")

summary_cenpa_broad_X <- cenpa_broad_X %>%
  summarise(
    count = n(),
    mean_value = mean(total_cenpa, na.rm = TRUE),
    sd_value = sd(total_cenpa, na.rm = TRUE),
    min_value = min(total_cenpa, na.rm = TRUE),
    max_value = max(total_cenpa, na.rm = TRUE),
    median_value = median(total_cenpa, na.rm = TRUE),
    IQR_value = IQR(total_cenpa, na.rm = TRUE)
  )

cenpa_53bp_mean <- cenpa_broad_X %>% filter(start >= 110.8e6, end <= 111.3e6)
sum(cenpa_53bp_mean$total_cenpa) #1,119,302
sum(cenpa_53bp_mean$cenpa_align) #74321
sum(cenpa_53bp_mean$cenpa_per_bases) #3847.848

cenpa_X_ctr_mean <- cenpa_broad_X %>% filter(start >= 38e6, end <= 50e6)
sum(cenpa_X_ctr_mean$total_cenpa) #378,408
sum(cenpa_X_ctr_mean$cenpa_align) #92515
sum(cenpa_X_ctr_mean$cenpa_per_bases) #1353.733
median(cenpa_X_ctr_mean$cenpa_align) #509

summary_cenpa_broad_X_ctr <- cenpa_broad_X %>% filter(start >= 38e6, end <= 50e6) %>%
  summarise(
    count = n(),
    mean_value = mean(total_cenpa, na.rm = TRUE),
    sd_value = sd(total_cenpa, na.rm = TRUE),
    min_value = min(total_cenpa, na.rm = TRUE),
    max_value = max(total_cenpa, na.rm = TRUE),
    median_value = median(total_cenpa, na.rm = TRUE),
    IQR_value = IQR(total_cenpa, na.rm = TRUE)
  )

summary_cenpa_broad_X_53bp <- cenpa_broad_X %>% filter(start >= 110.8e6, end <= 111.3e6) %>%
  summarise(
    count = n(),
    mean_value = mean(total_cenpa, na.rm = TRUE),
    sd_value = sd(total_cenpa, na.rm = TRUE),
    min_value = min(total_cenpa, na.rm = TRUE),
    max_value = max(total_cenpa, na.rm = TRUE),
    median_value = median(total_cenpa, na.rm = TRUE),
    IQR_value = IQR(total_cenpa, na.rm = TRUE)
  )

plot_CENPA_peaks <- cenpa_broad_X %>%
  ggplot(aes(x = start/1e6, y = total_cenpa)) +
  geom_point(size = 0.01) +
  labs(x = "Position (Mb)", y = "CENP-A broad level") +
  theme_classic() +
  scale_y_continuous(labels = comma) +
  geom_vline(xintercept = 26805429/1e6, linetype = "dashed", color = "blue", size = 0.1) +
  geom_vline(xintercept = 141053523/1e6, linetype = "dashed", color = "blue", size = 0.1) +
  geom_vline(xintercept = 150651618/1e6, linetype = "dashed", color = "blue", size = 0.1)
plot_CENPA_peaks

plot_CENPA_peaks_chr9 <- cenpa_broad %>%
  ggplot(aes(x = start/1e6, y = total_cenpa)) +
  geom_point(size = 0.01) +
  labs(x = "Position (Mb)", y = "CENP-A broad level") +
  theme_classic() +
  scale_y_continuous(labels = comma) +
  geom_vline(xintercept = 26805429/1e6, linetype = "dashed", color = "blue", size = 0.1) +
  geom_vline(xintercept = 141053523/1e6, linetype = "dashed", color = "blue", size = 0.1) +
  geom_vline(xintercept = 150651618/1e6, linetype = "dashed", color = "blue", size = 0.1)
plot_CENPA_peaks_chr9

####### all CENP-A peaks #######
# files:
cenpa_broad_final
cenpa_X_ctr_mean
cenpa_53bp_mean

cenpa_X_ctr_all <- cenpa_X_ctr_mean %>%
  mutate(
    size = 509,
    copy_row = ceiling(cenpa_align/size)) %>%
  uncount(copy_row) %>%
  group_by(start, end, total_cenpa) %>%
  mutate(cenpa_per_sat = total_cenpa / n(),
         query_repeat = "X_ctr") %>%
  ungroup


cenpa_broad_all <- cenpa_broad_final %>%
  mutate(chr = seqnames, start = start.x, end = end.x, cenpa_align = width) %>%
  select(chr, start, end, total_cenpa, cenpa_align, cenpa_per_bases, size, cenpa_per_sat, query_repeat)

cenpa_53bp <- cenpa_53bp_mean %>%
  mutate(
    size = 53,
    copy_row = ceiling(cenpa_align/size)) %>%
  uncount(copy_row) %>%
  group_by(start, end, total_cenpa) %>%
  mutate(cenpa_per_sat = total_cenpa / n(),
         query_repeat = "X_53bp") %>%
  ungroup

cenpa_28bp <- cenpa_broad_X %>%
  filter(start >= 150651618, end <= 150657971) %>%
  mutate(copies = 236) %>%
  uncount(copies) %>%
  mutate(cenpa_per_sat = total_cenpa / 236,
         size = round(cenpa_align / 236),
         query_repeat = "X_28bp")

cenpa_24bp <- cenpa_broad_X %>%
  filter(start >= 141051169, end <= 141072658) %>%
  mutate(copies = 892) %>%
  uncount(copies) %>%
  mutate(cenpa_per_sat = total_cenpa / 892,
         size = round(cenpa_align / 892),
         query_repeat = "X_24bp")

cenpa_18bp <- cenpa_broad_X %>%
  filter(start >= 26805318, end <= 26825765) %>%
  mutate(copies = 1140) %>%
  uncount(copies) %>%
  mutate(cenpa_per_sat = total_cenpa / 1140,
         size = round(cenpa_align / 1140),
         query_repeat = "X_18bp")

# cenpa_all_repeats <- rbind(cenpa_18bp, cenpa_24bp, cenpa_28bp, cenpa_53bp, cenpa_broad_all, cenpa_X_ctr_all)
cenpa_all_repeats <- rbind(cenpa_broad_all, cenpa_X_ctr_all)


plot_box_CENPA_all_repeats_perbase <- cenpa_all_repeats %>%
  ggplot(aes(x = query_repeat, y = cenpa_per_bases, color = query_repeat)) +
  geom_boxplot() +
  labs(
    x = "",
    y = "Normalised CENP-A level"
  ) +
  theme_classic() +
  scale_color_manual(values = c(
    "SATVII" = "#CCCC00",
    "SATII" = "#00008b",
    "SATIV" = "#f4895f",
    "SATVI" = "#7D0DC3",
    "SATI" = "#dd314b",
    "SATV" = "limegreen",
    "SATIII" = "darkgreen",
    "X_ctr" = "#DE3163",
    "X_53bp" = "#0047AB",
    "X_28bp" = "#89CFF0",
    "X_24bp" = "#5D3FD3",
    "X_18bp" = "#9FE2BF"), guide = "none")
plot_box_CENPA_all_repeats_perbase
ggsave(filename = paste0(dirout, "plot_box_CENPA_all_repeats_perbase_3.png"), plot = plot_box_CENPA_all_repeats_perbase, width = 4, height = 2.5, dpi = 900)

write_tsv(cenpa_all_repeats, file = paste0(dirout, "cenpa_all_repeats.tsv"))

KW_broad_df_all <- cenpa_all_repeats %>%
  group_by(query_repeat) %>%
  summarise(
    count = n(),
    total = sum(cenpa_per_sat, na.rm = TRUE),
    mean = mean(cenpa_per_bases, na.rm = TRUE),
    sd = sd(cenpa_per_bases, na.rm = TRUE),
    median = median(cenpa_per_bases, na.rm = TRUE),
    IQR = IQR(cenpa_per_bases, na.rm = TRUE))
print(KW_broad_df_all)


kruskal.test(cenpa_per_bases ~ query_repeat, data = cenpa_all_repeats)
#Kruskal-Wallis chi-squared = 6175.7, df = 7, p-value < 2.2e-16

pairwise.wilcox.test(cenpa_all_repeats$cenpa_per_bases, cenpa_all_repeats$query_repeat,
                     p.adjust.method = "fdr")
# SATIII  SATV    SATVII  X_18bp  X_24bp  X_28bp  X_53bp
# SATV   1.7e-12 -       -       -       -       -       -
#   SATVII 0.00077 1.6e-07 -       -       -       -       -
#   X_18bp < 2e-16 < 2e-16 < 2e-16 -       -       -       -
#   X_24bp < 2e-16 < 2e-16 < 2e-16 < 2e-16 -       -       -
#   X_28bp < 2e-16 < 2e-16 < 2e-16 < 2e-16 < 2e-16 -       -
#   X_53bp 8.8e-07 0.06843 1.7e-07 < 2e-16 < 2e-16 < 2e-16 -
#   X_ctr  < 2e-16 < 2e-16 < 2e-16 < 2e-16 < 2e-16 < 2e-16 < 2e-16

plot_box_CENPA_all_repeats_total_cenpa <- cenpa_all_repeats %>%
  ggplot(aes(x = query_repeat, y = cenpa_per_sat, color = query_repeat)) +
  geom_boxplot() +
  labs(
    x = "",
    y = "Normalised CENP-A level"
  ) +
  theme_classic() +
  scale_color_manual(values = c(
    "SATVII" = "#CCCC00",
    "SATII" = "#00008b",
    "SATIV" = "#f4895f",
    "SATVI" = "#7D0DC3",
    "SATI" = "#dd314b",
    "SATV" = "limegreen",
    "SATIII" = "darkgreen",
    "X_ctr" = "#DE3163",
    "53bp" = "#89CFF0",
    "28bp" = "#0047AB",
    "24bp" = "#5D3FD3",
    "18bp" = "#9FE2BF"), guide = "none")
plot_box_CENPA_all_repeats_total_cenpa

cenpa_broad_chr_raw <- cenpa_broad %>% filter(chr == 10)
cenpa_broad_chr_final <- cenpa_broad_final %>% filter(seqnames == 10)

stop()

####### BUBBLE PLOT PERC #######
## percentage identity, sequence identity of the bovine satellite repeats

plot_bubble_sats_perc_id <- rm_df_auto_sats %>%
  ggplot(aes(x = query_repeat, y = perc_id*100, size = copies_group, color = query_repeat)) +
  geom_point(alpha = 0.6) +
  scale_size_manual(name = "Copies",
                    values = c("1-100" = 1, "101-200" = 2, "201-300" = 3, "301-400" = 4, "401-500" = 5, "501-1000" = 6, "1000+" = 7)) +
  theme_classic() +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines")) +
  labs(y = "Sequence Identity (%)", color = "Repeat Satellites", x = "") +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines")) +
  scale_color_manual(values = c(
    "SATVII" = "#CCCC00",
    "SATII" = "#00008b",
    "SATIV" = "#f4895f",
    "SATVI" = "#7D0DC3",
    "SATI" = "#dd314b",
    "SATV" = "limegreen",
    "SATIII" = "darkgreen"),
    guide = "none")
plot_bubble_sats_perc_id
ggsave(filename = paste0(dirout, "plot_bubble_sats_perc_id.png"), plot = plot_bubble_sats_perc_id, width = 5, height = 2, dpi = 300)


plot_perc_SATI_SATVI <- rm_tsv_auto_sats %>%
  filter(repeat_type %in% c("SATI", "SATVI"),
         chr %in% c(9,10,21,23)) %>%
  ggplot() +
  geom_rect(aes(xmin = (qry_begin/1e6)+0.3, xmax = (qry_end/1e6)-0.3, ymin = (perc_id*100)+0.3, ymax = (perc_id*100)-0.3, fill = repeat_type)) +
  facet_wrap(~ chr, scales = "free_x") +
  theme_classic() +
  scale_fill_manual(values = c(
    "SATVI" = "#7D0DC3",
    "SATI" = "#dd314b",
    guide = "none")) +
  labs(x = "Position (Mb)", y = "Percent identity", fill = "Satellite repeat") +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray", size = 0.2) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "gray", size = 0.2)
plot_perc_SATI_SATVI
ggsave(filename = paste0(dirout, "plot_perc_SATI_SATVI.png"), plot = plot_perc_SATI_SATVI, width = 5, height = 3, dpi = 300)

plot_perc_SATII <- rm_tsv_auto_sats_ctr %>%
  filter(repeat_type %in% c("SATII"),
         chr %in% c(9,10,21,23)) %>%
  ggplot() +
  geom_rect(aes(xmin = (qry_begin/1e6), xmax = (qry_end/1e6), ymin = (perc_id*100)+0.3, ymax = (perc_id*100)-0.3, fill = repeat_type)) +
  facet_wrap(~ chr, scales = "free_x") +
  theme_classic() +
  scale_fill_manual(values = c(
    "SATII" = "#00008b",
    guide = "none")) +
  labs(x = "Position (Mb)", y = "Percent identity", fill = "Satellite repeat") +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray", size = 0.2)
plot_perc_SATII
ggsave(filename = paste0(dirout, "plot_perc_SATII.png"), plot = plot_perc_SATII, width = 5, height = 3, dpi = 300)


# ## 2023.10.14 look at the repeat blocks after the centromere:
#
# rm_afterctr_7 <- rm_df %>% filter(chr == 7, start >= 15252090, end.x <= 15378858)
# rm_afterctr_7_filter <- rm_afterctr_7 %>% filter(start >= 15253487, end.x <= 15265499)
# # repeatblock size = 126,768
# # end of simple repeat: 15253487-15265499 (12 012kb of AC rich) (9752 is the total align of the simple repeat)
# # > sum(rm_afterctr_7$align_len)
# # [1] 39,128 kb
# # conclusion: It's a simple_repeat
#
# rm_afterctr_14 <- rm_df %>% filter(chr == 14, start >= 3727365, end.x <= 7271880)
# # align_len of repeat is: 1,258,788
# # size of block is: 3,544,515
# # unannotated: 2,285,727
#
# rm_afterctr_22 <- rm_df %>% filter(chr == 22, start >= 11566283, end.x <= 12143750)
# # align_len of repeat is: 181,893
# # size of block is: 577,467
# # contains 5S from 11807240-12025738 = 218,498
# # conclusion: 5S repeatblock
#
# rm_afterctr_25 <- rm_df %>% filter(chr == 25, start >= 6795081, end.x <= 9738224) %>% arrange(start) %>% mutate(diff = start - lag(end.x))
# # align_len of repeat is: 1,160,404
# # size of block is: 2,943,143
# #
#
# rm_afterctr_27 <- rm_df %>% filter(chr == 27, start >= 11965008, end.x <= 17580859) %>% arrange(start) %>% mutate(diff = start - lag(end.x))
# # align_len of repeat is: 1,744,500
# # size of block is: 5,615,851
# ## 2024.10.13 done (see sup for result)

rm_df_9 <- rm_df %>% filter(chr == 9)

####### SUMMARY FOR PAPER #######
## For the paper figure summary of the repeats in the entire genome
summary_genome_entire_class <- summary_repeats_genome_entire_class %>%
  mutate(category = ifelse(total_repeats > 1e6,
                           ifelse(class == "DNA", "Others", class),
                           "Others")) %>%
  # Calculate total_repeats for the "Non-repetitive" category
  bind_rows(data.frame(
    class = "Non-repetitive",
    total_repeats = 3142409774 - sum(summary_repeats_genome_entire_class$total_repeats),
    total_repeats_genome = (3142409774 - sum(summary_repeats_genome_entire_class$total_repeats))/3142409774* 100,
    category = "Non-repetitive"
  )) %>%
  group_by(category) %>%
  summarise(
    perc_genome = round(sum(total_repeats_genome),2),
    repeats_bp = sum(total_repeats),
    class = paste0(class, collapse = ",")
  )

plot_summary_genome_entire_class <- summary_genome_entire_class %>%
  ggplot(aes(x = perc_genome, y = reorder(category, perc_genome), fill = category)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(x = "Percentage in the genome", y = "") +
  theme(legend.position = "none") +
  geom_text(aes(label = paste0(perc_genome, "%")), hjust = -0.2, color = "black") +
  xlim(0, 60)
plot_summary_genome_entire_class

####### GC CONTENT #######
### GC content of the satellite repeats ###

gc_sats <- "/Users/polenpineda/Documents/TulixWagyu/repeat_masker/GC_satellites_rm.bed"
gc_sats_df <- read_delim(gc_sats, comment = "#", col_names = c("chr","start","end","family","align_len",
                                                               "strand","class","subclass","divergence","col","perc_AT","perc_GC","numA","numC","numG","numT",
                                                               "num_N","num_oth","seq_len"))

gc_sats_df <- merge(gc_sats_df, sat_len_df, by = "family")
gc_sats_df_filter <- gc_sats_df %>%
  mutate(copies = round(align_len/size),1) %>%
  filter(copies > 0) %>%
  select(chr,family,copies,align_len,perc_GC,numG,numC) %>%
  mutate(copies_group = cut(copies,
                            breaks = c(0, 100, 200, 300, 400, 500, 1000, Inf),
                            labels = c("1-100", "101-200", "201-300", "301-400", "401-500", "501-1000", "1000+")))

summary_gc_sats_df_filter_sats <- gc_sats_df_filter %>%
  group_by(family) %>%
  summarise(
    count = n(),
    total_align = sum(align_len),
    num_GC = sum(numG, numC),
    total_copies = sum(copies),
    average_perc_GC = sum(perc_GC * copies) / total_copies)


plot_gc_sats_df_filter <- gc_sats_df_filter %>%
  ggplot(aes(x = family, y = perc_GC *100, size = copies_group, color = family)) +
  geom_point(alpha = 0.6) +
  scale_size_manual(name = "Copies",
                    values = c("1-100" = 1, "101-200" = 2,
                               "201-300" = 3, "301-400" = 4,
                               "401-500" = 5, "501-1000" = 6,
                               "1000+" = 7)) +
  theme_classic() +
  theme(legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines")) +
  labs(y = "GC content (%)", color = "Repeat Satellites", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c(
    "SATVII" = "#CCCC00",
    "SATII" = "#00008b",
    "SATIV" = "#f4895f",
    "SATVI" = "#7D0DC3",
    "SATI" = "#dd314b",
    "SATV" = "limegreen",
    "SATIII" = "darkgreen"))
plot_gc_sats_df_filter


## FOR FURTHER ANALYSIS OF THE CENTROMERE REPEATS, USE rm_tsv_auto or rm_df_auto
#### Centromeric satellite repeats

ctr_rm_tsv_auto <- rm_tsv_auto %>%
  filter(end.x < end.y, chr %in% c(1:29)) %>%
  mutate(copies_aligned = round(align_len/size))

####### X_53bp #######
X_53bp <- rm_df %>% filter(chr == "X", start >110.8e6, end <111.3e6)


rm_df_X <- rm %>% filter(chr == "X", class == "Satellite")
rm_tsv_X <- rm_tsv_raw %>% filter(chr == "X", family == "Satellite/centr")

stop()
## files for the bovine satellite repeats

####### saving files #######
## contains the centromeric region boundary
ctr_all_repclusters_df_filter_final
write_tsv(ctr_all_repclusters_df_filter_final, file = paste0(dirout, "wagyu_centromere_location.txt"))

## summary of all the repeats in the genome
summary_repeats_genome_entire
summary_repeats_genome_chr
summary_repeats_genome_unplaced
summary_repeats_genome_entire_class
summary_repeats_genome_entire_repeattype
summary_repeats_genome_unplaced_class
summary_repeats_genome_chr_repeattype
summary_repeats_genome_entire_subclass
summary_genome_entire_class #(this is the data for the figure to be put into the sup file)
plot_summary_genome_entire_class


####### X palindrome vs consensei #######
## blast with the consensus repeat modeler result

blast_consensi <- "/Users/polenpineda/Documents/TulixWagyu/blast/out/consensi.classified.autosome.sats_vs_X_38Mb_to_50Mb_palindrome_95perc_sequence_left.85.blast.edited.tsv"
blast_consensi <- read_delim(blast_consensi, col_names = c("qry", "ref", "perc_id", "align_len",
                                                           "mismatches", "gap_opens", "qry_start", "qry_end",
                                                           "ref_start", "ref_end", "eval", "bitscore"))
consensi_fai <- "/Users/polenpineda/Documents/TulixWagyu/blast/db/consensi.classified.autosome.sats.fa.fai"
consensi_fai <- read_delim(consensi_fai, col_names = c("ref","size"))
consensi_fai <- consensi_fai %>% select(ref,size)

blast_consensi <- merge(blast_consensi, consensi_fai, by = "ref")

blast_consensi_filter <- blast_consensi %>%
  mutate(align_cover = qry_end - qry_start + 1,
         align_len_per = align_len/size) %>%
  filter(perc_id > 95, align_len_per > 0.80) %>%
  separate(ref, into = c("ref","family"), sep = "#")

summary_blast_consensi_filter <- blast_consensi_filter %>%
  group_by(ref) %>%
  summarise(count = n(),
            mean_perc_id = mean(perc_id),
            mean_align_len = mean(align_len),
            sum_align_len = sum(align_len),
            size = mean(size),
            arms = paste0(qry, collapse = ", "),
            family = first(family),
            percentage = (sum_align_len / sum(blast_consensi_filter$align_len)) * 100)

summary_blast_consensi_filter_family <- blast_consensi_filter %>%
  group_by(family) %>%
  summarise(count = n(),
            mean_perc_id = mean(perc_id),
            mean_align_len = mean(align_len),
            sum_align_len = sum(align_len),
            size = mean(size),
            arms = paste0(qry, collapse = ", "),
            family = first(family),
            percentage = (sum_align_len / sum(blast_consensi_filter$align_len)) * 100)

options(scipen = 999)
plot_blast_consensi_filter_family <- summary_blast_consensi_filter_family %>%
  ggplot(aes(x = sum_align_len, y = family, fill = family)) +
  geom_col() +
  theme_classic() +
  geom_text(aes(label = paste0(round(percentage, 2), "%"),
                x = sum_align_len + max(sum_align_len) * 0.02), # Adjust position
            hjust = 0) +
  labs(x = "Total alignment length (bp)", y = "Repeat family") +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 4.5e6), labels = scales::comma)
plot_blast_consensi_filter_family

####### X chromosome repeat #######

## autosome satellite repeats in the X chromosome?

rm_tsv_X_raw <- rm_tsv_raw %>% filter(chr == "X") %>% mutate(align_len = qry_end - qry_begin)
rm_tsv_X_raw_sats <- rm_tsv_X_raw %>% filter(family == "Satellite/centr")

plot_X_raw_sats <- rm_tsv_X_raw_sats %>%
  ggplot(aes(x = qry_begin/1e6, y = align_len, color = query_repeat)) +
  geom_point(size = 0.01) +
  labs(x = "Position (Mb)", y = "Alignment length", title = "X chromosome") +
  geom_vline(xintercept = c(110.8, 111.3), size = 0.1, color = "blue")
plot_X_raw_sats

plot_X_raw_sats_ctr <- rm_tsv_X_raw_sats %>%
  ggplot(aes(x = qry_begin/1e6, y = align_len, color = query_repeat)) +
  geom_point(size = 0.01) +
  labs(x = "Position (Mb)", y = "Alignment length", title = "X chromosome") +
  ylim(0,300) +
  geom_vline(xintercept = 110.8, linetype = "dashed", color = "blue", size = 0.1) +
  geom_vline(xintercept = 111.3, linetype = "dashed", color = "blue", size = 0.1) +
  geom_vline(xintercept = 38, linetype = "dashed", color = "red", size = 0.1) +
  geom_vline(xintercept = 50, linetype = "dashed", color = "red", size = 0.1)
plot_X_raw_sats_ctr

rm_tsv_X_raw_sats_ctr <- rm_tsv_X_raw_sats %>% filter(qry_begin >= 38e6, qry_end <= 50e6)
rm_tsv_X_raw_sats_53bp <- rm_tsv_X_raw_sats %>% filter(qry_begin >= 110.8e6, qry_end <= 111.3e6)

X_ctr_SATIV <- rm_tsv_X_raw_sats_ctr %>% filter(query_repeat == "SATIV")
SATIV_axis <- 1:3808
counts_SATIV_axis <- lapply(SATIV_axis, function(x) {
  subset <- X_ctr_SATIV[X_ctr_SATIV$new_rep_begin <= x & X_ctr_SATIV$new_rep_end >= x, ]
  count <- nrow(subset)
  perc_id <- mean(subset$perc_div)
  return(data.frame(x = x, count = count, perc_id = perc_id))
})
X_ctr_SATIV_axis_per1bp <- do.call(rbind, counts_SATIV_axis)

plot_X_ctr_SATIV_axis_per1bp <- ggplot(X_ctr_SATIV_axis_per1bp, aes(x = x, y = count, color = perc_id)) +
  geom_line(size = 0.5) +
  scale_color_viridis(option = "H", limits = c(85, 100)) +
  theme_classic() +
  labs(title = "X chromosome centromere alignment to SATIV",
       x = "alignment position",
       y = "base copy number",
       color = "sequence identity")
plot_X_ctr_SATIV_axis_per1bp

X_SATIV_noctr <- rm_tsv_X_raw_sats %>% filter(query_repeat == "SATIV",
                                              qry_begin > 50e6 | qry_end < 38e6)
SATIV_axis <- 1:3808
counts_SATIV_axis <- lapply(SATIV_axis, function(x) {
  subset <- X_SATIV_noctr[X_SATIV_noctr$new_rep_begin <= x & X_SATIV_noctr$new_rep_end >= x, ]
  count <- nrow(subset)
  perc_id <- mean(subset$perc_div)
  return(data.frame(x = x, count = count, perc_id = perc_id))
})
X_SATIV_noctr_axis_per1bp <- do.call(rbind, counts_SATIV_axis)

plot_X_SATIV_noctr_axis_per1bp <- ggplot(X_SATIV_noctr_axis_per1bp, aes(x = x, y = count, color = perc_id)) +
  geom_line(size = 0.5) +
  scale_color_viridis(option = "H", limits = c(85, 100)) +
  theme_classic() +
  labs(title = "X chromosome alignment to SATIV",
       x = "alignment position",
       y = "base copy number",
       color = "sequence identity")
plot_X_SATIV_noctr_axis_per1bp


x_rm <- "/Users/polenpineda/Documents/TulixWagyu/repeat_masker/X.fa.out.bed"
x_rm_df <- read_delim(x_rm, col_names = c("chr","start","end","family","score","orientation","class","subclass","divergence","linkage_id"))

x_rm_df <- x_rm_df %>% mutate(align_len = end - start +1)

summary_x_rm_df <- x_rm_df %>%
  group_by(class) %>%
  summarise(total_align = sum(align_len),
            percent_in_X = sum(align_len)/173607011*100)
sum(x_rm_df$align_len)
# 103517007/173607011*100 = 59.6272

summary2_x_rm_df <- x_rm_df %>%
  group_by(subclass) %>%
  summarise(total_align = sum(align_len),
            percent_in_X = sum(align_len)/173607011*100)


stop()

####### ARS-UCD2.0 #######

ars_rm <- "/Users/polenpineda/Documents/TulixWagyu/repeat_masker/concat_ars_rm.bed"
ars_rm_df <- read_delim(ars_rm, col_names = c("chr","start","end","family","score","orientation","class","subclass","divergence","linkage_id"))
ars_rm_df <- ars_rm_df %>%
  mutate(align_len = end - start +1,
       class_family = paste0(class, "/", subclass),
       repeat_type = ifelse(class == "Satellite", family, class))

ars_rm_summary_repeats_genome_entire <- ars_rm_df %>%
  summarise(total_repeats = sum(align_len),
            total_repeats_genome = total_repeats/size_genome*100)

ars_rm_summary_repeats_genome_entire_class <- ars_rm_df %>%
  group_by(class) %>%
  summarise(total_repeats = sum(align_len),
            total_repeats_genome = total_repeats/size_genome*100)

ars_rm_summary_repeats_genome_autosomes <- ars_rm_df %>%
  filter(chr %in% c(1:29)) %>%
  group_by(class) %>%
  summarise(total_repeats = sum(align_len),
            total_repeats_genome = total_repeats/(size_autosomes)*100)
