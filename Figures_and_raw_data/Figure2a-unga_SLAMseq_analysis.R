library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)


gene_expr_data <- read.csv("SLAMseq_reanalysis/bhat2023-SLAMseq_expr_data.csv")

conversion_data <- read.csv("SLAMseq_reanalysis/bhat2023-SLAMseq-TC_conversion.csv")

ung_expr_data <- gene_expr_data %>% select(-classification) %>% select(-sub.classification) %>%
  filter(Gene.ID %in% c("unga", "ungb")) %>%
  pivot_longer(cols = (-"Gene.ID"), names_to=c("blank", "rep", "time_point"), 
               names_sep="_", 
               values_to="rpm") %>%
  select(-blank)

ung_expression <- ggplot(ung_expr_data, aes(x=time_point, y=rpm, fill=Gene.ID)) +
  geom_boxplot(position=position_dodge(1)) + theme_classic() +
  scale_fill_manual(values = c("lightskyblue2", "indianred3"))

ung_conversion_data <- conversion_data %>% select(-MZ_cluster) %>%
  filter(gene %in% c("unga", "ungb")) %>%
  pivot_longer(cols = (-"gene"), names_to=c("blank1", "blank2", "time_point"), 
               names_sep="_", 
               values_to="fractionTC") %>%
  select(-blank1, -blank2) 

ung_conversion <-ggplot(ung_conversion_data, aes(x=time_point, y=fractionTC, color=gene)) +
  geom_point() + theme_classic() +
  scale_color_manual(values = c("lightskyblue2", "indianred3")) + ylim(0,1)

ggarrange(ung_expression, nrow=1,
          labels=NULL, align = "hv")

ggsave("20250521-ung_SLAMseq.pdf", units = "cm", height = 9, width = 12)

