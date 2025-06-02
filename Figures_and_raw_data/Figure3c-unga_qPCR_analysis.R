#R code for the qPCR analysis

# Load required libraries
library(ggplot2)
library(tidyverse)

# Read the data
unga_qPCR_data <- read_csv("unga_qPCR.csv")


# Filter for relevant genes and create sample type column
unga_qPCR_data_filtered <- unga_qPCR_data %>%
  filter(`Gene Name` %in% c("lsm12b", "unga")) %>%
  mutate(
    Sample_Type = case_when(
      grepl("wt", `Sample Name`, ignore.case = TRUE) ~ "WT",
      grepl("unga", `Sample Name`, ignore.case = TRUE) ~ "unga_mutant",
      TRUE ~ "Other"
    ),
    Cq = as.numeric(Cq)  # Ensure Cq is numeric
  ) %>%
  filter(Sample_Type != "Other")

# Calculate Delta Cq (target gene - housekeeping gene)
# Get housekeeping gene values
housekeeping <- unga_qPCR_data_filtered %>%
  filter(`Gene Name` == "lsm12b") %>%
  select(`Sample Name`, Cq_housekeeping = Cq)

# Get target gene values
target <- unga_qPCR_data_filtered %>%
  filter(`Gene Name` == "unga") %>%
  select(`Sample Name`, Cq, Sample_Type)

# Merge and calculate Delta Cq
delta_cq_data <- target %>%
  inner_join(housekeeping, by = "Sample Name") %>%
  mutate(
    Delta_Cq = Cq - Cq_housekeeping,
    Relative_Expression = 2^(-Delta_Cq)
  )

# Calculate fold change relative to WT
wt_mean_rel_exp <- delta_cq_data %>%
  filter(Sample_Type == "WT") %>%
  pull(Relative_Expression) %>%
  mean()

delta_cq_data <- delta_cq_data %>%
  mutate(Fold_Change = Relative_Expression / wt_mean_rel_exp)

delta_cq_data$Sample_Type <- factor(delta_cq_data$Sample_Type, levels = c("WT", "unga_mutant"))  
df_filtered$Sample_Type <- factor(df_filtered$Sample_Type, levels = c("WT", "unga_mutant"))  

# Display results
print("Delta Cq calculations:")
print(delta_cq_data)

# Summary statistics
print("Relative expression summary:")
delta_cq_data %>%
  group_by(Sample_Type) %>%
  summarise(
    mean_Delta_Cq = mean(Delta_Cq),
    sd_Delta_Cq = sd(Delta_Cq),
    mean_Rel_Exp = mean(Relative_Expression),
    sd_Rel_Exp = sd(Relative_Expression),
    .groups = "drop"
  )

# Create and save the standalone relative expression plot
ggplot(delta_cq_data, aes(x = Sample_Type, y = Relative_Expression)) +
  geom_boxplot(aes(fill = Sample_Type), alpha = 0.7) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
  scale_y_log10() +
  labs(
    #title = "unga relative expression (2^-ΔCq)",
    x = "Strain",
    y = "unga relative expression"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "none"
  ) +
  scale_fill_brewer(type = "qual", palette = "Set2")


# Save relative expression plot in multiple formats
ggsave("unga_relative_expression.pdf", units = "cm", width = 6, height = 8, )
