library(tidyverse)
library(ggplot2)
library(wesanderson)
library(ggpubr)

unga_fert_data <-read.csv("unga_fertility_data.csv", header=TRUE)

unga_fert_data_tidy <- unga_fert_data %>% mutate(cross = str_c(male,"_",female)) %>%
  mutate(fertilized_eggs = total_eggs - unfertilized_eggs) %>%
  select(c(cross, fertilized_eggs, unfertilized_eggs)) %>%
  pivot_longer(cols = !cross,
               names_to = "eggs",
               values_to = "number") 


ggplot(unga_fert_data_tidy, aes(fill=factor(eggs, levels=c('unfertilized_eggs', 'fertilized_eggs')), 
                                x=factor(cross, level=c("wt_wt", "unga_unga", "wt_unga", 
                                                        "unga_wt")),
                                y=number)) + 
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "bottom") +
  labs(x = "Genotype", y="Ratio", fill="Eggs") +
  scale_fill_brewer(palette = "Paired")

unga_fert_data_tidy_2 <- unga_fert_data %>% mutate(cross = str_c(male,"_",female)) %>%
  mutate(fertilized_eggs = total_eggs - unfertilized_eggs) %>%
  mutate(fert_ratio = fertilized_eggs/total_eggs)


my_comparisons <- list( c("wt_wt", "unga_unga"), c("wt_wt", "wt_unga"), c("wt_wt", "unga_wt") )

ggplot(unga_fert_data_tidy_2, aes(x=factor(cross, level=c("wt_wt", "unga_unga", "wt_unga", 
                                                        "unga_wt")),
                                y=fert_ratio,
                                fill = cross)) +
  geom_boxplot() + geom_jitter() +
  theme_classic() +
  labs(x = "Genotype", y="Ratio of fertilized eggs") +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons) +
  scale_fill_manual(values = wes_palette("AsteroidCity3"))

ggsave("20250527-unga_fert_plot.pdf", units = "cm", width=7.5, height=10 )
