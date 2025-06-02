# Analysis of apoptosis (AcridineOrange foci) following DEB treatment
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Reading data and tidying up a bit for further analysis
ao_data <- read.csv("unga-DEB_AcridineOrange_data.csv", header=TRUE)  %>% 
  unite("GenTreat", Genotype:Treatment, remove = FALSE)

# Relevant statistical comparisons
my_comparisons <- list( c("wt_cntrl", "unga_cntrl"),
                        c("wt_DEB", "unga_DEB"),
                        c("wt_cntrl", "wt_DEB"),
                        c("unga_cntrl", "unga_DEB"))

# Plotting all the data, faceting by different thresholds
ggplot(ao_data, aes(x = factor(GenTreat, level=c("wt_cntrl", "unga_cntrl",
                                                 "wt_DEB", "unga_DEB")), 
                    y = Count, 
                    fill = factor(Genotype, level=c("wt", "unga")))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.8) +
  scale_fill_manual(values = c("wt" = "#2E86AB", "unga" = "#F24236")) +
  facet_wrap(vars(Min_size), scales = "free", nrow=5) +
  labs(
    title = "Apoptotic cells",
    x = "Treatment and genotype",
    y = "Number of AO foci"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")

ggsave("unga_AO_foci_DEB_all.pdf", units = "cm", width=20, height=20)


# Plotting the data with particle size larger than 49
ao_data_49 <- read.csv("unga-DEB_AcridineOrange_data.csv", header=TRUE)  %>% 
  unite("GenTreat", Genotype:Treatment, remove = FALSE) %>%
  filter(Min_size == 49 )


ggplot(ao_data_49, aes(x = factor(GenTreat, level=c("wt_cntrl", "unga_cntrl",
                                                 "wt_DEB", "unga_DEB")), 
                    y = Count, 
                    fill = factor(Genotype, level=c("wt", "unga")))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.15, alpha = 0.8) +
  scale_fill_manual(values = c("wt" = "#2E86AB", "unga" = "#F24236")) +
  labs(
    title = "Apoptotic cells",
    x = "Treatment and genotype",
    y = "Number of AO foci"
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons)

ggsave("unga_AO_foci_DEB_49.pdf", units = "cm", width=8, height=12)
