# Combined Sperm Analysis - Concentration and Motility
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Read both datasets
conc_data <- read.csv("sperm_concentration_data.csv")
motility_data <- read.csv("sperm_motility_data.csv")

# Convert concentration to millions for better visualization
conc_data$Concentration_millions <- conc_data$Concentration_spz_ml / 1e6

my_comparisons <- list( c("wt", "unga"))

# Create a comprehensive figure with all parameters
# Concentration plot
p_conc <- ggplot(conc_data, aes(x = factor(Genotype, level=c("wt", "unga")), 
                                 y = Concentration_millions, fill = Genotype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = c("wt" = "#2E86AB", "unga" = "#F24236")) +
  ylim(0,1200) +
  labs(
    title = "Sperm Concentration",
    x = "Genotype",
    y = "Concentration (millions/ml)"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons)

# Motility plots
p_tot_mot <- ggplot(motility_data, aes(x = factor(Genotype, level=c("wt", "unga")), 
                            y = Tot_Motility_pct, fill = Genotype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = c("wt" = "#2E86AB", "unga" = "#F24236")) +
  ylim(0,100) +
  labs(
    title = "Total Motility",
    x = "Genotype",
    y = "Percentage (%)"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons)

p_prog_mot <- ggplot(motility_data, aes(x = factor(Genotype, level=c("wt", "unga")), 
                            y = Prog_Motility_pct, fill = Genotype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = c("wt" = "#2E86AB", "unga" = "#F24236")) +
  ylim(0,90) +
  labs(
    title = "Progressive Motility",
    x = "Genotype",
    y = "Percentage (%)"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons)

p_immotile <- ggplot(motility_data, aes(x = factor(Genotype, level=c("wt", "unga")), 
                            y = Immotile_pct, fill = Genotype)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = c("wt" = "#2E86AB", "unga" = "#F24236")) +
  ylim(0,30) +
  labs(
    title = "Immotile Sperm",
    x = "Genotype",
    y = "Percentage (%)"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons)

# Combine all plots
combined_plot <- ggarrange(p_conc, p_tot_mot, p_prog_mot, p_immotile, 
                             ncol = 4, nrow = 1)

# Save the combined plot
ggsave("combined_sperm_analysis.pdf", combined_plot, units= "cm", width = 15, height = 8, dpi = 300)


