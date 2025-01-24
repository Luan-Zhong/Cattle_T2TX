library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(RColorBrewer)

## read data
ARS_size <- read_delim("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/tuli_wagyu/telomeres/ARS_UCD_v2.0.fa.fai", col_names = c("id","size","start_loc","A","B"))
ARS_size <- ARS_size %>% select(id, size) %>% mutate(size = size, endloc = cumsum(size), begin = endloc - size)
ARS_tsv <- read_delim("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/tuli_wagyu/telomeres/ARS_UCD2.0_tsv_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)
ARS_bg <- read_delim("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/tuli_wagyu/telomeres/ARS_UCD2.0_bg_tidk-search_telomeric_repeat_windows.bedgraph", col_names = c("id","start","end","value"))

hap1_size <- read_delim("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/tuli_wagyu/telomeres/assembly.haplotype1.fasta_tmp_asm.fasta.fai", col_names = c("id","size","start_loc","A","B"))
hap1_size <- hap1_size %>% select(id, size) %>% mutate(size = size, endloc = cumsum(size), begin = endloc - size)
hap1_tsv <- read_delim("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/tuli_wagyu/telomeres/assembly.haplotype1.fasta_tmp_asm_tsv_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)
hap1_bg <- read_delim("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/tuli_wagyu/telomeres/assembly.haplotype1.fasta_tmp_asm_bedgraph_tidk-search_telomeric_repeat_windows.bedgraph", col_names = c("id","start","end","value"))

hap2_size <- read_delim("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/tuli_wagyu/telomeres/assembly.haplotype2.fasta_tmp_asm.fasta.fai", col_names = c("id","size","start_loc","A","B"))
hap2_size <- hap2_size %>% select(id, size) %>% mutate(size = size, endloc = cumsum(size), begin = endloc - size)
hap2_tsv <- read_delim("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/tuli_wagyu/telomeres/assembly.haplotype2.fasta_tmp_asm_tsv_tidk-search_telomeric_repeat_windows.tsv", col_names = TRUE)
hap2_bg <- read_delim("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/tuli_wagyu/telomeres/assembly.haplotype2.fasta_tmp_asm_bedgraph_tidk-search_telomeric_repeat_windows.bedgraph", col_names = c("id","start","end","value"))


## filter unplaced
ARS_size <- ARS_size %>% filter(id %in% 1:29)
ARS_tsv <- ARS_tsv %>% filter(id %in% 1:29)
ARS_bg <- ARS_bg %>% filter(id %in% 1:29)

hap1_size <- hap1_size %>% filter(id %in% 1:29)
hap1_tsv <- hap1_tsv %>% filter(id %in% 1:29)
hap1_bg <- hap1_bg %>% filter(id %in% 1:29)

hap2_size <- hap2_size %>% filter(id %in% 1:29)
hap2_tsv <- hap2_tsv %>% filter(id %in% 1:29)
hap2_bg <- hap2_bg %>% filter(id %in% 1:29)

# ## filter non-telomeric sites
# ARS_filt <- ARS_tsv %>% group_by(id) %>% slice(c(1:2, (n()-1):n())) %>% 
#   filter(forward_repeat_number >= 50 | reverse_repeat_number >= 50, id!="X") %>% group_by(id, tmp = window == 20000) %>% 
#   summarise(window = min(window), across(contains("repeat_number"), sum)) %>%
#   ungroup() %>% select(-tmp)
# ARS_filt <- merge(ARS_filt, ARS_size, by="id") %>% 
#   mutate(reverse_size = reverse_repeat_number*6, forward_size = forward_repeat_number*6, for_start = size - forward_size, rev_end = 1 + reverse_size,
#          loc_start = ifelse(window==20000, 1, for_start), loc_end = ifelse(window==20000, rev_end, size),
#          telomere_type = ifelse(forward_repeat_number > 50,"TTAGGG","CCCTAA"),
#          telomere_count = ifelse(telomere_type=="TTAGGG",forward_repeat_number,reverse_repeat_number),
#          telomere_size = telomere_count*6, asm="ARS_UCD2.0")
# ARS_filt_final <- ARS_filt %>% select(id,telomere_type, telomere_count, telomere_size, loc_start, loc_end, asm)
# 
# hap1_filt <- hap1_tsv %>% group_by(id) %>% slice(c(1:2, (n()-1):n())) %>% 
#   filter(forward_repeat_number >= 50 | reverse_repeat_number >= 50, id!="X") %>% group_by(id, tmp = window == 20000) %>% 
#   summarise(window = min(window), across(contains("repeat_number"), sum)) %>%
#   ungroup() %>% select(-tmp)
# hap1_filt <- merge(hap1_filt, ARS_size, by="id") %>% 
#   mutate(reverse_size = reverse_repeat_number*6, forward_size = forward_repeat_number*6, for_start = size - forward_size, rev_end = 1 + reverse_size,
#          loc_start = ifelse(window==20000, 1, for_start), loc_end = ifelse(window==20000, rev_end, size),
#          telomere_type = ifelse(forward_repeat_number > 50,"TTAGGG","CCCTAA"),
#          telomere_count = ifelse(telomere_type=="TTAGGG",forward_repeat_number,reverse_repeat_number),
#          telomere_size = telomere_count*6, asm="Wagyu")
# hap1_filt_final <- hap1_filt %>% select(id,telomere_type, telomere_count, telomere_size, loc_start, loc_end, asm)
# 
# hap2_filt <- hap2_tsv %>% group_by(id) %>% slice(c(1:2, (n()-1):n())) %>% 
#   filter(forward_repeat_number >= 50 | reverse_repeat_number >= 50, id!="X") %>% group_by(id, tmp = window == 20000) %>% 
#   summarise(window = min(window), across(contains("repeat_number"), sum)) %>%
#   ungroup() %>% select(-tmp)
# hap2_filt <- merge(hap2_filt, ARS_size, by="id") %>% 
#   mutate(reverse_size = reverse_repeat_number*6, forward_size = forward_repeat_number*6, for_start = size - forward_size, rev_end = 1 + reverse_size,
#          loc_start = ifelse(window==20000, 1, for_start), loc_end = ifelse(window==20000, rev_end, size),
#          telomere_type = ifelse(forward_repeat_number > 50,"TTAGGG","CCCTAA"),
#          telomere_count = ifelse(telomere_type=="TTAGGG",forward_repeat_number,reverse_repeat_number),
#          telomere_size = telomere_count*6, asm="Tuli")
# hap2_filt_final <- hap2_filt %>% select(id,telomere_type, telomere_count, telomere_size, loc_start, loc_end, asm)
# 
# ### CREATING TABLES
# ## create summaries
# all_telo_size <- genome_list <- list(ARS_UCD2.0 = ARS_filt_final$telomere_size,
#                                      Wagyu = hap1_filt_final$telomere_size,
#                                      Tuli = hap2_filt_final$telomere_size)
# all_asm <- rbind(ARS_filt_final, hap1_filt_final, hap2_filt_final)
# 
# chr <- data.frame(c(1:29))
# names(chr) <- "id"
# all_asm$id <- as.numeric(all_asm$id)
# all_asm_chr <- merge(chr, all_asm, by="id", all = TRUE) %>% select(id, telomere_count, asm) %>%
#   group_by(id, asm) %>% summarise(total_telomere_count = sum(telomere_count))
# 
# telomere_counts_per_chr <- all_asm_chr %>% spread(key = asm, value = total_telomere_count, fill = NA) %>% select(id, ARS_UCD2.0, Wagyu, Tuli)
# names(telomere_counts_per_chr) <- c("Chromosomes", "ARS_UCD2.0", "Wagyu", "Tuli")
# 
# ## clean table
# ARS_tlm_count <- ARS_filt  %>% select(id,forward_repeat_number,reverse_repeat_number,asm) %>% group_by(id) %>% summarise(forward_repeat_number = max(forward_repeat_number), reverse_repeat_number=max(reverse_repeat_number), asm = first(asm)) %>%
#   arrange(as.numeric(id))
# hap1_tlm_count <- hap1_filt  %>% select(id,forward_repeat_number,reverse_repeat_number,asm) %>% group_by(id) %>% summarise(forward_repeat_number = max(forward_repeat_number), reverse_repeat_number=max(reverse_repeat_number), asm = first(asm)) %>%
#   arrange(as.numeric(id))
# hap2_tlm_count <- hap2_filt  %>% select(id,forward_repeat_number,reverse_repeat_number,asm) %>% group_by(id) %>% summarise(forward_repeat_number = max(forward_repeat_number), reverse_repeat_number=max(reverse_repeat_number), asm = first(asm)) %>%
#   arrange(as.numeric(id))
# 
# ### CREATING PLOTS
# ## plot telomere mean
# all_asm_compute_mean <- all_asm_chr %>% na.omit() %>% group_by(asm) %>% summarise_at(vars(total_telomere_count),list(mean = mean, sd = sd)) %>%
#   as.data.frame()
# ggplot(all_asm_compute_mean, aes(x=asm, y=mean)) +
#   geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd)) +
#   geom_point() +
#   labs(title = "Telomere sizes", x = "Genome assembies", y = "Mean (unit)") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# asm <- c("ARS_UCD2.0", "Wagyu", "Tuli")
# colors <- c("ARS_UCD2.0"="#694a30", "Wagyu"="#a4bf01","Tuli"="#257ca3")
# all_asm_size_per_chromosome <- telomere_counts_per_chr %>% gather(key = "Assemblies", value = "Unit", -Chromosomes)
# all_asm_size_per_chromosome_plot <- all_asm_size_per_chromosome %>%
#   ggplot(aes(x = factor(Chromosomes), y = Unit, fill = Assemblies)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(title = "Comparison of Telomere Counts", x = "Chromosomes", y = "Unit") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_fill_manual(values = colors, breaks = asm)
# print(all_asm_size_per_chromosome_plot)

## plot telomere signals  per chromosomes
### notes
# 1. to make a gap between chromosomes in the bedgraph
## a. identify the last line of each group
## b. add another line: start with the last value of the previous number + 40000 with 0 value

gray_colors <- rep(c("#A9A9A9", "#D3D3D3"), length.out = 30)

test <- ARS_bg %>% mutate(row = row_number())
test[156645,4] = 10 #chr14
test[99747,4] = 10 #chr8
test[180576,4] = 10 #chr17
test[120646,4] = 10 #chr10
test[156644,4] = 10 #chr14
test[120647,4] = 10 #chr10
test[77370,4] = 10 #chr6
test[77369,4] = 10 #chr6
test[120648,4] = 10 #chr10
test[207236,4] = 10 #chr21

ARS_bg_edited <- test %>% select(-row)

ARS_bg <- ARS_bg_edited %>%
  group_by(id) %>%
  group_modify(~ .x %>%
                 summarise(start = max(start), end = max(end), value = 0) %>%
                 add_row(.x, .)
  ) %>%
  ungroup()
ARS_bg2 <- left_join(ARS_bg, ARS_size, by = "id")

ARS_bg_final <- ARS_bg2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
ARS_bg_final$id <- factor(ARS_bg_final$id, levels = c(1:29))
ARS_bg_rects <- ARS_size %>% select(id, begin, endloc) %>% group_by(id) %>% summarise(begin = min(begin), endloc = max(endloc)) %>% arrange(as.numeric(id))
ARS_bg_rects$id <- factor(ARS_bg_rects$id, levels = c(1:29))
ARS_bg_plot <- ggplot() +
  geom_rect(data = ARS_bg_rects, aes(xmin = begin, xmax = endloc, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
  geom_step(data = ARS_bg_final, aes(x = start, y = value, color = id)) +
  scale_fill_manual(values = gray_colors) +
  theme_classic() +
  ylim(0, 2000) +
  xlim(0,2700000000) +
  xlab("") +
  ylab("")
print(ARS_bg_plot)

hap1_bg <- hap1_bg %>%
  group_by(id) %>%
  group_modify(~ .x %>%
                 summarise(start = max(start), end = max(end), value = 0) %>%
                 add_row(.x, .)
  ) %>%
  ungroup()
hap1_bg2 <- left_join(hap1_bg, hap1_size, by = "id")
hap1_bg_final <- hap1_bg2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
hap1_bg_final$id <- factor(hap1_bg_final$id, levels = c(1:29))
hap1_bg_rects <- hap1_size %>% select(id, begin, endloc) %>% group_by(id) %>% summarise(begin = min(begin), endloc = max(endloc)) %>% arrange(as.numeric(id))
hap1_bg_rects$id <- factor(hap1_bg_rects$id, levels = c(1:29))
hap1_bg_plot <- ggplot() +
  geom_rect(data = hap1_bg_rects, aes(xmin = begin, xmax = endloc, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
  geom_step(data = hap1_bg_final, aes(x = start, y = value, color = id)) +
  scale_fill_manual(values = gray_colors) +
  theme_classic() +
  ylim(0, 2000) +
  xlim(0,2700000000) +
  xlab("") +
  ylab("")
print(hap1_bg_plot)

# hap2_bg <- hap2_bg %>%
#   group_by(id) %>%
#   group_modify(~ .x %>%
#                  summarise(start = max(start), end = max(end), value = 0) %>%
#                  add_row(.x, .)
#   ) %>%
#   ungroup()
# hap2_bg2 <- left_join(hap2_bg, hap2_size, by = "id")
# hap2_bg_final <- hap2_bg2 %>% mutate(start_new = start+begin, end_new = end+begin) %>% select(id,start=start_new, end=end_new,value)
# hap2_bg_final$id <- factor(hap2_bg_final$id, levels = c(1:29))
# hap2_bg_rects <- hap2_size %>% select(id, begin, endloc) %>% group_by(id) %>% summarise(begin = min(begin), endloc = max(endloc)) %>% arrange(as.numeric(id))
# hap2_bg_rects$id <- factor(hap2_bg_rects$id, levels = c(1:29))
# hap2_bg_plot <- ggplot() +
#   geom_rect(data = hap2_bg_rects, aes(xmin = begin, xmax = endloc, ymin = -Inf, ymax = Inf, fill = id), alpha = 0.4) +
#   geom_step(data = hap2_bg_final, aes(x = start, y = value, color = id)) +
#   scale_fill_manual(values = gray_colors) +
#   theme_classic() +
#   ylim(0, 2000) +
#   xlim(0,2700000000) +
#   xlab("") +
#   ylab("")
# print(hap2_bg_plot)

# ## all in one plot
# ARS_bg_final <- ARS_bg_plot + theme(legend.position = "none") +labs(title = "ARS_UCD2.0")
# hap1_bg_final <- hap1_bg_plot + theme(legend.position = "none") +labs(title = "Wagyu")
# hap2_bg_final <- hap2_bg_plot + theme(legend.position = "none") +labs(title = "Tuli")
# 
# legend_b <- get_legend(ARS_bg_final + theme(legend.position = "right"))
# 
# telomere_bg_plot_paper <- ggarrange(ARS_bg_final, hap1_bg_final, hap2_bg_final,
#                                     hjust=-0.8,
#                                     legend.grob = legend_b, legend = "right",
#                                     ncol = 1, nrow = 3)
# telomere_bg_plot_assemblies <- annotate_figure(telomere_bg_plot_paper, top = text_grob("Telomeric signals in the assemblies", face = "bold", size = 18),
#                                                 left = "Count of TTAGGG/CCCTAA", bottom = "Region")
# print(telomere_bg_plot_assemblies)
# ggsave("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/telomeres/telomere_bg_plot_assemblies.png",telomere_bg_plot_assemblies, width = 12, height = 8, units = "in")

options(scipen = 0)
## Wagyu and ARSUCD plot only
ARS_bg_final <- ARS_bg_plot + theme(legend.position = "none") +labs(title = "ARS_UCD2.0")
hap1_bg_final <- hap1_bg_plot + theme(legend.position = "none") +labs(title = "Wagyu")

legend_b <- get_legend(ARS_bg_final + theme(legend.position = "right"))

telomere_bg_plot_paper <- ggarrange(ARS_bg_final, hap1_bg_final,
                                    hjust=-0.8,
                                    legend.grob = legend_b, legend = "right",
                                    ncol = 1, nrow = 2)
telomere_bg_plot_final_paper <- annotate_figure(telomere_bg_plot_paper, top = text_grob("Telomeric signals in the assemblies", face = "bold", size = 18),
                                                left = "Count of TTAGGG/CCCTAA", bottom = "Region")
print(telomere_bg_plot_final_paper)
ggsave("/Users/polenpineda/Documents/TulixWagyu/genomics_meeting/telomeres/telomere_bg_plot_final.png",telomere_bg_plot_final_paper, width = 12, height = 6, units = "in")
