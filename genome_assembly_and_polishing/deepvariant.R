# Rscript for comparing merqury qv unpolished vs polished with deepvariant

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

dirout <- '/Users/polenpineda/Documents/TulixWagyu/polishing/'
data <- '/Users/polenpineda/Documents/TulixWagyu/polishing/merqury_per_chr_frozen.csv'
df <- read_delim(data, col_names = TRUE)

df_long <- df %>%
  pivot_longer(cols = starts_with("mat") | starts_with("pat"), names_to = "type", values_to = "value") %>%
  mutate(value = as.numeric(value))  %>%
  drop_na()

chr <- c(as.character(1:29), "X", "Y")
df_long$chr <- factor(df_long$chr, levels = chr)

df_mat <- df_long %>%
  filter(type %in% c("mat_QV", "mat_QV_deepv"))
df_pat <- df_long %>%
  filter(type %in% c("pat_QV", "pat_QV_deepv"))

plot_mat <- ggplot(df_mat, aes(x = as.factor(chr), y = value, color = type, group = type)) +
  geom_line() +
  labs(
       x = "Chromosome",
       y = "QV") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(50, 68), breaks = seq(50, 68, by = 1))
  

plot_pat <- ggplot(df_pat, aes(x = as.factor(chr), y = value, color = type, group = type)) +
  geom_line() +
  labs(
       x = "Chromosome",
       y = "QV") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(50, 68), breaks = seq(50, 68, by = 1))

ggsave(filename = "plot_maternal_comparison_QV_withDeepV.png", plot = plot_mat, path = paste0(dirout), width = 8, height = 4, units = "in", dpi = 300)
ggsave(filename = "plot_paternal_comparison_QV_withDeepV.png", plot = plot_pat, path = paste0(dirout), width = 8, height = 4, units = "in", dpi = 300)

plot_mat + plot_pat
