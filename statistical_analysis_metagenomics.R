





# Loading required packages
library(BiocManager)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(nlme)
library(FSA)
library(ggpubr)
library(rstatix)
library(broom)
library(forcats)
library(ape)
library(plyr)
library(vegan)
library(cowplot)
library(randomcoloR)
library(phyloseq)
library("DESeq2")
library("metacoder")
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(patchwork)
library("writexl")
library(shadowtext)
library(reshape2)
library(magrittr)
library("ggpubr")
library("microbiomeMarker")



# Setup the working directory to the path of the files (otu_table, tables_samples, etc)
# Windows
setwd("C:/Users/andre/OneDrive - IHMT-NOVA/final_data_and_images/scripts/wastewater_metagenomic_analysis_data_and_scripts/script_results")

# MAC
setwd("/Users/andre_santos/Library/CloudStorage/OneDrive-IHMT-NOVA/final_data_and_images/scripts/wastewater_metagenomic_analysis_data_and_scripts/script_results")


# Load up the data
tables_samples_for_paper <- read.xlsx("./table_samples_for_paper.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
tables_samples <- read.xlsx("./table_samples.xlsx", sheet = 1, startRow = 1, colNames = TRUE) # == sample_data
non_phage_otu_table_long <- read.xlsx("./otu_table_non_phage_long.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
non_phage_otu_table <- read.xlsx("./otu_table_non_phage.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
phage_otu_table <- read.xlsx("./otu_table_phage.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
phage_otu_table_long <- read.xlsx("./otu_table_phage_long.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

all_samples_for_R <- read.xlsx("./taxid_df.xlsx", sheet = 1, startRow = 1, colNames = TRUE)

undefined <- read.xlsx("./taxid_df_undefined.xlsx", sheet = 1, startRow = 1, colNames = TRUE)





# Setup the directory for the outputs (plots, excel files, etc)
# Windows
setwd("C:/Users/andre/OneDrive - IHMT-NOVA/final_data_and_images/scripts/wastewater_metagenomic_analysis_data_and_scripts/script_results/R_plots/diversity_analysis_plots_not_normalized")

# MAC
setwd("/Users/andre_santos/Library/CloudStorage/OneDrive-IHMT-NOVA/final_data_and_images/scripts/wastewater_metagenomic_analysis_data_and_scripts/script_results/R_plots/diversity_analysis_plots_not_normalized")


tables_samples_for_paper
colnames(tables_samples_for_paper)

png(file = "./barplot_raw_reads_per_sample.png", bg="white", width = 12, height = 8, units = "in", res = 300)
ggplot(tables_samples_for_paper, aes(x = Sample, y = Number.of.raw.reads)) +
  geom_bar(stat = "identity", show.legend = TRUE) 
dev.off()

png(file = "./barplot_raw_reads_per_sample_2.png", bg="white", width = 12, height = 8, units = "in", res = 300)
ggplot(tables_samples_for_paper, aes(x = Sample.Number, y = Number.of.raw.reads)) +
  geom_bar(stat = "identity", show.legend = TRUE) +  # position="fill", 
  labs(x = "Replicate",
       y = "Number of raw reads") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),         
        axis.title.x = element_text(size = 16), # Adjust the size of the x-axis title
        axis.text.y = element_text(size = 14),  
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18)) + 
  facet_grid(Treament_Plant~Step , scales="fixed", space="fixed")
dev.off()


"
#######################
tables_samples_for_paper.xlsx
#######################
Samples: 32
"
# Normality test for distribution (normal if p-value > 0.05)
shapiro.test(tables_samples_for_paper$Number.of.raw.reads)  # W = 0.74071, p-value = 3.752e-06
shapiro.test(tables_samples_for_paper$Concentration) # W = 0.92013, p-value = 0.02096


min(tables_samples_for_paper$Number.of.raw.reads)
max(tables_samples_for_paper$Number.of.raw.reads)
median(tables_samples_for_paper$Number.of.raw.reads)
iqr(tables_samples_for_paper$Number.of.raw.reads)

min(tables_samples_for_paper$Concentration)
max(tables_samples_for_paper$Concentration)
median(tables_samples_for_paper$Concentration)
iqr(tables_samples_for_paper$Concentration)


"
From the output, the two p-values are smaller than the significance 
level 0.05 implying that the distribution of the data are significantly 
different from normal distribution.
We can't assume the normality.

If the data are not normally distributed, it’s recommended to use the 
non-parametric correlation, including Spearman and Kendall rank-based correlation tests.
"

# Correlation between raw and library concentration
cor.test(tables_samples_for_paper$Number.of.raw.reads, tables_samples_for_paper$Concentration,
         method=c('spearman')) # rho = 0.3465195, p-value = 0.05203



ggplot(tables_samples_for_paper, aes(x = Number.of.raw.reads, y = Concentration)) +
  geom_point(stat = "identity", show.legend = TRUE) +  # position="fill", 
  labs(x = "Replicate",
       y = "Number of raw reads") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),         
        axis.title.x = element_text(size = 16), # Adjust the size of the x-axis title
        axis.text.y = element_text(size = 14),  
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18)) + 
  facet_grid(Treament_Plant~Step , scales="fixed", space="fixed")





colnames(tables_samples_for_paper)[colnames(tables_samples_for_paper) == "Replicate/Pool"] ="Replicate_Pool"

tables_samples_for_paper %>% group_by(Replicate_Pool) %>% 
  dplyr::summarise(median = median(Number.of.raw.reads), IQR = iqr(Number.of.raw.reads))  # Have to use the dplyr package, it was using the plyr, lololololol

"
Is the mean raw reads for the replicates different to the raw reads for the pools?
"
wilcox.test(Number.of.raw.reads ~ Replicate_Pool, data = tables_samples_for_paper)  # might need to be paired, check later
"
The p-value is 0.4037 This means that we do not have enough evidence to reject 
the null hypothesis.
This implies that there is no statistically significant difference 
in the distributions of raw reads between the two sample types based on our data.
"
raw_reads_median_iqr <- tables_samples_for_paper %>%
  group_by(Replicate_Pool) %>%
  dplyr::summarize(
    median = median(Number.of.raw.reads),
    IQR = IQR(Number.of.raw.reads)
  )
raw_reads_median_iqr


# # A tibble: 2 × 3
# Replicate_Pool median      IQR
# <chr>           <dbl>    <dbl>
#   1 Pool           881004 1129542.
# 2 Replicate      906068  978000.

tables_samples_for_paper %>%
  dplyr::summarize(
    median = median(Number.of.raw.reads),
    IQR = IQR(Number.of.raw.reads)
  )

"Correlation between total viral reads vs raw reads (PER SOFTWARE)"
############## NON-PHAGES ##############
# PERFORM PEARSON CORRELATION ANALYSIS
spearman_result <- cor.test(tables_samples$Number_of_all_viral_reads,
                            tables_samples$Raw_reads,
                            method = 'spearman'); spearman_result

# Extract p-values
spearman_p_value <- spearman_result$p.value



tables_samples = tables_samples %>%
  mutate(Software_long = case_when(
    Software == "INSA" ~ "INSaFLU-TELEVIR",
    Software == "GD" ~ "Genome Detective",
    Software == "CZID" ~ "CZ.ID",
    Software == "KR" ~ "Kraken2",
    TRUE ~ as.character(Software)  # Keep the original value if no match
  ))

# Perform Pearson correlation tests for each software and store the results
correlation_results <- tables_samples %>%
  group_by(Software_long) %>%
  dplyr::summarize(spearman_cor = cor(Number_of_all_viral_reads, Raw_reads, method = 'spearman'),
                   spearman_p_value = cor.test(Number_of_all_viral_reads, Raw_reads, method = 'spearman')$p.value)

correlation_results



png(file = "./correlation_all_viral_vs_raw_reads_per_software.png", bg="white", width = 12, height = 8, units = "in", res = 300)
# Create scatter plot with ggplot2
p <- ggplot(tables_samples, aes(x = Number_of_all_viral_reads, y = Raw_reads)) +
  geom_point() +
  labs(x = "Number of all viral reads",
       y = "Number of Raw reads") +
  theme_minimal() +
  facet_wrap(~ Software_long, scales = "free") + 
  theme(axis.text.x = element_text(size = 14),  # Adjust the size of the x-axis text
        axis.title.x = element_text(size = 16), # Adjust the size of the x-axis title
        axis.text.y = element_text(size = 14),  
        axis.title.y = element_text(size = 16),
        strip.text = element_text(size = 18))

# Annotate the plot with the Pearson correlation coefficient
p + geom_text(data = correlation_results, aes(x = Inf, y = Inf, 
                                              label = paste("Spearman r:", format(spearman_cor, digits = 3), 
                                                            "\np-value:", format(spearman_p_value, digits = 3))), 
              hjust = 1.1, vjust = 1.5, size = 6, color = "red")
dev.off()



##############################################
######## Viral identification results ########
##############################################

min(tables_samples$Number_of_all_viral_reads) # 3
max(tables_samples$Number_of_all_viral_reads) # 288464
median(tables_samples$Number_of_all_viral_reads) # 1092.5
iqr(tables_samples$Number_of_all_viral_reads) # 2449.25

# Sample with highest viral reads
tables_samples[order(-tables_samples$Number_of_all_viral_reads), ] %>%
  head(3) %>%
  select(SampleID, Software, Number_of_all_viral_reads)

# 49    D1_GD       GD                    288464
# 34  C1_INSA     INSA                    201198
# 51    D1_KR       KR                    158295

# Sample with lowest viral reads
tables_samples[order(-tables_samples$Number_of_all_viral_reads), ] %>%
  tail(3) %>%
  select(SampleID, Software, Number_of_all_viral_reads)

# 79     E5_KR       KR                         4
# 110    G5_KR       KR                         3
# 118    H2_KR       KR                         3

# Normality test for distribution (normal if p-value > 0.05)
shapiro.test(tables_samples$Number_of_all_viral_reads) # W = 0.25259, p-value < 2.2e-16


# Correlation between raw and number of all_viral_reads
cor.test(tables_samples$Raw_reads, 
         tables_samples$Number_of_all_viral_reads,
         method=c('spearman')) # rho = 0.3916755, p-value = 5.753e-06


# Create the proportions
tables_samples = tables_samples %>%
  mutate(Proportion_raw_vs_viral = Number_of_all_viral_reads / Raw_reads,
         Viral_Proportion_vs_all_reads = Number_of_all_viral_reads / Number_of_all_reads,
         Non_phage_proportion = Number_of_non_phage_reads / Number_of_all_viral_reads,
         Phage_proportion = Number_of_only_phage_reads / Number_of_all_viral_reads)



colnames(tables_samples)
proportions_table = tables_samples[c("SampleID",
                                     "Software",
                                     "Raw_reads",
                                     "Proportion_raw_vs_viral",
                                     "Viral_Proportion_vs_all_reads",
                                     "Non_phage_proportion",
                                     "Phage_proportion")]

write_xlsx(proportions_table, "proportions_table.xlsx")
min(proportions_table$Proportion_raw_vs_viral) # 5.897339e-06
max(proportions_table$Proportion_raw_vs_viral) # 0.3119487

# Sample with highest proportion of raw vs viral
proportions_table[order(-proportions_table$Proportion_raw_vs_viral), ] %>%
  head(3) %>%
  select(SampleID, Software, Proportion_raw_vs_viral)

# Sample with lowest proportion of raw vs viral
proportions_table[order(-proportions_table$Proportion_raw_vs_viral), ] %>%
  tail(3) %>%
  select(SampleID, Software, Proportion_raw_vs_viral)


kruskal.test(tables_samples$Number_of_all_viral_reads, tables_samples$Software)
"
P-value = 1.687e-08, 
chi-squared = 39.059
meaning we reject the null hypothesis.
There is strong evidence to suggest that at least one of the group medians is 
significantly different from the others.


Since the Kruskal-Wallis test indicates that there are significant differences 
among the groups, you may want to conduct post-hoc tests to determine which 
specific groups differ from each other. One common approach is to use pairwise 
comparisons with a correction for multiple testing, such as the Dunn test.
"

library(dunn.test)

# Perform Dunn's test
dunn.test(tables_samples$Number_of_all_viral_reads, tables_samples$Software, method = "bonferroni")

# Calculate median and IQR for the number of all viral reads by software
median_iqr <- tables_samples %>%
  group_by(Software) %>%
  dplyr::summarize(
    median = median(Number_of_all_viral_reads),
    IQR = IQR(Number_of_all_viral_reads)
  )

# Print the result
print(median_iqr)

write_xlsx(median_iqr, "median_iqr_per_software.xlsx")



##############################################
############# Phage vs Non-Phage #############
##############################################
# Phage
min_phage_prop <- min(tables_samples$Number_of_only_phage_reads, na.rm = TRUE); min_phage_prop # 0
max_phage_prop <- max(tables_samples$Number_of_only_phage_reads, na.rm = TRUE); max_phage_prop # 1
median_phage_prop <- median(tables_samples$Number_of_only_phage_reads, na.rm = TRUE); median_phage_prop # 0.4159091
iqr_phage_prop <- iqr(tables_samples$Number_of_only_phage_reads, na.rm = TRUE); iqr_phage_prop # 0.6206303

# Sample with highest phage reads
top_samples_phage <- tables_samples[order(-tables_samples$Phage_proportion), ] %>%
  head(3) %>%
  select(SampleID, Software, Phage_proportion)

top_samples_phage


# Shapiro-Wilk test
shapiro_test_phage<- shapiro.test(tables_samples$Phage_proportion); shapiro_test_phage

# Spearman correlation# Spearman correlationshapiro_test
spearman_corr_phage <- cor.test(tables_samples$Raw_reads, tables_samples$Number_of_only_phage_reads, method = "spearman"); spearman_corr_phage


summary_table_phage <- data.frame(
  Statistic = c("Min Phage Proportion", "Max Phage Proportion", "Median Phage Proportion", "IQR Phage Proportion", 
                "Shapiro-Wilk W", "Shapiro-Wilk p-value", 
                "Spearman Correlation (rho)", "Spearman Correlation p-value"),
  Value = c(min_phage_prop, max_phage_prop, median_phage_prop, iqr_phage_prop, 
            shapiro_test_phage$statistic, shapiro_test_phage$p.value, 
            spearman_corr_phage$estimate, spearman_corr_phage$p.value)
)

# Display top samples separately
top_samples_table_phage <- top_samples_phage

print(summary_table_phage)
print(top_samples_table_phage)

write.xlsx(list(Summary = summary_table_phage, TopSamples = top_samples_table_phage), 
           file = "phage_summary_output.xlsx")


# Non-phage
min_non_phage_prop <- min(tables_samples$Non_phage_proportion, na.rm = TRUE); min_non_phage_prop
max_non_phage_prop <- max(tables_samples$Non_phage_proportion, na.rm = TRUE); max_non_phage_prop
median_non_phage_prop <- median(tables_samples$Non_phage_proportion, na.rm = TRUE); median_non_phage_prop
iqr_non_phage_prop <- iqr(tables_samples$Non_phage_proportion, na.rm = TRUE); iqr_non_phage_prop

# Sample with highest phage reads
top_samples_non_phages <- tables_samples[order(-tables_samples$Non_phage_proportion), ] %>%
  head(3) %>%
  select(SampleID, Software, Non_phage_proportion)

# Shapiro-Wilk test
shapiro_test_non_phage <- shapiro.test(tables_samples$Non_phage_proportion); shapiro_test_non_phage

# Spearman correlation
spearman_corr_non_phage <- cor.test(tables_samples$Raw_reads, tables_samples$Number_of_non_phage_reads, method = "spearman"); spearman_corr_non_phage


summary_table_non_phages <- data.frame(
  Statistic = c("Min Non-Phage Proportion", "Max Non-Phage Proportion", "Median Non-Phage Proportion", "IQR Non-Phage Proportion", 
                "Shapiro-Wilk W", "Shapiro-Wilk p-value", 
                "Spearman Correlation (rho)", "Spearman Correlation p-value"),
  Value = c(min_non_phage_prop, max_non_phage_prop, median_non_phage_prop, iqr_non_phage_prop, 
            shapiro_test_non_phage$statistic, shapiro_test_non_phage$p.value, 
            spearman_corr_non_phage$estimate, spearman_corr_non_phage$p.value)
)

# Display top samples separately
top_samples_table_non_phages <- top_samples_non_phages

print(summary_table_non_phages)
print(top_samples_table_non_phages)

write.xlsx(list(Summary = summary_table_non_phages, TopSamples = top_samples_table_non_phages), 
           file = "non_phage_summary_output.xlsx")

# Create a merged table for phage and non-phage

summary_table_viral <- data.frame(
  Statistic = c("Min Proportion", "Max Proportion", "Median Proportion", "IQR Proportion", 
                "Shapiro-Wilk W", "Shapiro-Wilk p-value", 
                "Spearman Correlation (rho)", "Spearman Correlation p-value"),
  Phage = c(min_phage_prop, max_phage_prop, median_phage_prop, iqr_phage_prop, 
            shapiro_test_phage$statistic, shapiro_test_phage$p.value, 
            spearman_corr_phage$estimate, spearman_corr_phage$p.value),
  
  Non_Phage = c(min_non_phage_prop, max_non_phage_prop, median_non_phage_prop, iqr_non_phage_prop, 
            shapiro_test_non_phage$statistic, shapiro_test_non_phage$p.value, 
            spearman_corr_non_phage$estimate, spearman_corr_non_phage$p.value)
)

summary_table_viral
write.xlsx(list(Summary = summary_table_viral), 
           file = "summary_statistics_phage_and_non_phage.xlsx")



colnames(tables_samples)

# Create the TreatmentPlant_Step columns
tables_samples$Treatment_Step <- paste(tables_samples$treatment_plant, "-", tables_samples$collection_step)

tables_samples = tables_samples %>%
  mutate(Software = case_when(
    Software == "INSA" ~ "INSaFLU-TELEVIR",
    Software == "GD" ~ "Genome Detective",
    Software == "CZID" ~ "CZ.ID",
    Software == "KR" ~ "Kraken2",
    TRUE ~ as.character(Software)  # Keep the original value if no match
  ))


to_melt= c("Number_of_non_phage_reads",
           "Number_of_only_phage_reads")

to_keep_non_melted = colnames(tables_samples) %>% setdiff(., to_melt); to_keep_non_melted

df_total_viral <- melt(tables_samples,  # Reshape data from wide to long format
                       id.vars = to_keep_non_melted, 
                       variable.name = "Total_viral_Proportions", 
                       value.name = "Total_viral_Value")

# Set the factor levels to ensure correct order
df_total_viral$Total_viral_Proportions <- factor(df_total_viral$Total_viral_Proportions, levels = to_melt)


# Define custom color palette
custom_colors_total_viral <- c("Number_of_non_phage_reads" = "#92C5DE",
                               "Number_of_only_phage_reads" = "#D6604D")

custom_labels_total_viral <- c("Number_of_non_phage_reads" = "Proportion of non phage viral reads",
                               "Number_of_only_phage_reads" = "Proportion of phage viral reads")

png(file = "proportion_of_phage_vs_non_phage.png", bg="white",width = 12, height = 17,  units = "in", res = 300)
ggplot(df_total_viral, aes(x = Total_viral_Value, y = Replicate, fill = Total_viral_Proportions)) +
  geom_bar(stat = "identity", show.legend = TRUE, position="fill") +  # position="fill", 
  scale_fill_manual(values = custom_colors_total_viral, labels = custom_labels_total_viral) +
  labs(title = "",
       x = "Proportion",
       y = "Replicate") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 17),  # Change this value to adjust size
        axis.title.y = element_text(size = 17),
        strip.text = element_text(size = 13),) + 
  theme(legend.title=element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1)) +  # Set legend to a single column
  facet_grid(Treatment_Step~Software, scales="fixed", space="fixed") +
  scale_y_discrete(limits=rev)
dev.off()


# Print the result
print(median_iqr)


df_total_viral$Number_of_phage_reads = 
  df_total_viral$Number_of_full_id_phage_reads + 
  df_total_viral$Number_of_partially_identified_phage_reads 


df_total_viral$Number_of_non_phage_reads = 
  df_total_viral$Number_of_full_id_non_phage_reads + 
  df_total_viral$Number_of_partially_identified_non_phage_reads 



##############################################
########### Classification level #############
##############################################
# Transforming tables_samples_df from wide to long format to do the acumulated percentual for each family
to_melt = c("Number_of_full_id_non_phage_reads", 
                       "Number_of_partially_identified_non_phage_reads", 
                       "Number_of_full_id_phage_reads",
                       "Number_of_partially_identified_phage_reads") 

to_keep_non_melted = colnames(tables_samples) %>% setdiff(., to_melt); to_keep_non_melted

# Reshape data from wide to long format
tables_samples_df_for_4_props <- melt(tables_samples, 
                                                id.vars = to_keep_non_melted, 
                                                variable.name = "Proportions",
                                                value.name = "Value")

# Set the factor levels to ensure correct order
tables_samples_df_for_4_props$Proportions <- factor(tables_samples_df_for_4_props$Proportions, 
                                                              levels = to_melt)

# Define custom color palette
custom_colors <- c("Number_of_full_id_non_phage_reads" = "#2166AC",
                             "Number_of_partially_identified_non_phage_reads" = "#92C5DE",
                             "Number_of_full_id_phage_reads" = "#D6604D",
                             "Number_of_partially_identified_phage_reads" = "#FDDBC7"
                             )

# Define custom labels
custom_labels <- c("Number_of_full_id_non_phage_reads" = "Fully Identified Non-Phage Reads",
                  "Number_of_partially_identified_non_phage_reads" = "Partially Identified Non-Phage Reads",
                  "Number_of_full_id_phage_reads" = "Fully Identified Phage Reads",
                  "Number_of_partially_identified_phage_reads" = "Partially Identified Phage Reads"
                  )

png(file = "./plot_sample_by_workflow_by_classification_level.png", bg="white",width = 12, height = 17,  units = "in", res = 300)
ggplot(tables_samples_df_for_4_props, aes(x = Value, y = Replicate, fill = Proportions)) +  # Family, Abundance, .desc=TRUE)
  geom_bar(position="fill", stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = custom_colors, labels = custom_labels) +
  labs(title = "",
       x = "Relative frequency",
       y = "Replicate") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 17),  # Change this value to adjust size
        axis.title.y = element_text(size = 17),
        strip.text = element_text(size = 13),) + 
  theme(legend.title=element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2)) +  # Set legend to a single column
  facet_grid(Treatment_Step~Software, scales="fixed", space="fixed") +
  scale_y_discrete(limits=rev)
dev.off()


colnames(tables_samples)

# Calculate median and IQR for the number of all viral reads by software
median_iqr_phage_vs_non_phage <- tables_samples %>%
  group_by(Software) %>%
  dplyr::summarize(
    median_phage = median(Phage_proportion),
    median_non_phage = median(Non_phage_proportion),
    IQR_phage = IQR(Phage_proportion),
    IQR_non_phage = IQR(Non_phage_proportion)
  )

median_iqr_phage_vs_non_phage


colnames(tables_samples)

# Create the proportions
tables_samples = tables_samples %>%
  mutate(Fully_ID_Non_Phage_Prop = Number_of_full_id_non_phage_reads / Number_of_non_phage_reads,
         Partial_ID_Non_Phage_Prop = Number_of_partially_identified_non_phage_reads / Number_of_non_phage_reads,
         Fully_ID_Phage_Prop = Number_of_full_id_phage_reads / Number_of_only_phage_reads,
         Partial_ID_Phage_Prop = Number_of_partially_identified_phage_reads / Number_of_only_phage_reads)


# Check normality of data
shapiro.test(tables_samples$Fully_ID_Non_Phage_Prop)  # W = 0.92492, p-value = 3.127e-06
shapiro.test(tables_samples$Partial_ID_Non_Phage_Prop)  # W = 0.92492, p-value = 3.127e-06
shapiro.test(tables_samples$Fully_ID_Phage_Prop)  # W = 0.72121, p-value = 1.14e-13
shapiro.test(tables_samples$Partial_ID_Phage_Prop)  # W = 0.72121 p-value = 1.14e-13


# Test for differences in identification power amonst softwares
kruskal_full_id_non_phage = kruskal.test(tables_samples$Number_of_full_id_non_phage_reads, tables_samples$Software); kruskal_full_id_non_phage
# p-value < 0.0001
kruskal_partial_id_non_phage = kruskal.test(tables_samples$Number_of_partially_identified_non_phage_reads, tables_samples$Software); kruskal_partial_id_non_phage
# p-value < 0.0001
kruskal_full_id_phage = kruskal.test(tables_samples$Number_of_full_id_phage_reads, tables_samples$Software); kruskal_full_id_phage
# p-value = 0.0001801
kruskal_partial_id_phage = kruskal.test(tables_samples$Number_of_partially_identified_phage_reads, tables_samples$Software); kruskal_partial_id_phage
# p-value = 0.0001801

summary_table_kruskal_classification_level <- data.frame(
  Viral_Type = c("Non-phage", "", "", "", "", "",
           "Phage", "", "", "", "", ""),
  Classification_level = c("Fully identified", "", "", 
                           "Partially identified", "", "",
                           "Fully identified", "", "", 
                           "Partially identified", "", ""),
  Statistic = c("chi-squared", "df", "p-value"),
  Value = c(kruskal_full_id_non_phage$statistic,
            kruskal_full_id_non_phage$parameter,
            kruskal_full_id_non_phage$p.value,
            kruskal_partial_id_non_phage$statistic,
            kruskal_partial_id_non_phage$parameter,
            kruskal_partial_id_non_phage$p.value,
            kruskal_full_id_phage$statistic,
            kruskal_full_id_phage$parameter,
            kruskal_full_id_phage$p.value,
            kruskal_partial_id_phage$statistic,
            kruskal_partial_id_phage$parameter,
            kruskal_partial_id_phage$p.value)
)

summary_table_kruskal_classification_level

write.xlsx(list(Summary = summary_table_kruskal_classification_level), 
           file = "summary_table_kruskal_classification_level.xlsx")







##############################################
############ Abundance Analysis ##############
##############################################

# Define the function
process_dataframe <- function(df, abundance_threshold = 10) {
  
  # Add TreatmentPlant column
  df <- df %>%
    mutate(TreatmentPlant = case_when(
      substr(sample_code_number, 1, 1) %in% c("A", "B", "C", "D") ~ "WWTP 1",
      substr(sample_code_number, 1, 1) %in% c("E", "F", "G", "H") ~ "WWTP 2",
      TRUE ~ NA_character_
    ))
  
  # Add Step column
  df <- df %>%
    mutate(Step = case_when(
      substr(sample_code_number, 1, 1) == "A" ~ "Step 1",
      substr(sample_code_number, 1, 1) == "B" ~ "Step 2",
      substr(sample_code_number, 1, 1) == "C" ~ "Step 3",
      substr(sample_code_number, 1, 1) == "D" ~ "Step 4",
      substr(sample_code_number, 1, 1) == "E" ~ "Step 1",
      substr(sample_code_number, 1, 1) == "F" ~ "Step 2",
      substr(sample_code_number, 1, 1) == "G" ~ "Step 3",
      substr(sample_code_number, 1, 1) == "H" ~ "Step 4",
      TRUE ~ NA_character_
    ))
  
  # Add the full Step_Code
  df <- df %>%
    mutate(Step_Code = paste(TreatmentPlant, "-", Step))
  
  # Filter based on abundance threshold
  df_filtered <- df %>%
    filter(abundance >= abundance_threshold)
  
  # Return the processed dataframe
  return(df_filtered)
}



# Check the unique sample-workflow initially
length(unique(non_phage_otu_table_long$sample))  # 123
length(grep("INSA", unique(non_phage_otu_table_long$sample))) # 31 # stays the same

# Example usage with df_non_phage dataframe and a threshold of 10
non_phage_otu_table_long <- process_dataframe(non_phage_otu_table_long, abundance_threshold = 10)

length(unique(non_phage_otu_table_long$sample))  # 123 # stayed the same
length(grep("INSA", unique(non_phage_otu_table_long$sample))) # 31 # stayed the same


non_phage_otu_table_long <- non_phage_otu_table_long %>%
  mutate(family = str_replace(family, "^f__", ""),   # Remove "f__" from family
         genus = str_replace(genus, "^g__", ""),
         species = str_replace(species, "^s__", ""))     # Remove "g__" from genus


# Non phage families
length(unique(non_phage_otu_table_long$family)) # 59

# Non phage genus
length(unique(non_phage_otu_table_long$genus)) #  98


# Calculate the relative abundance
non_phage_otu_table_long <- non_phage_otu_table_long %>%
  mutate(Relative_Abundance = abundance / total_reads_per_sample)


# Ensure sample_number is a factor
non_phage_otu_table_long <- non_phage_otu_table_long %>%
  mutate(sample_number = factor(sample_number, levels = c(1, 2, 3, 5)))

non_phage_otu_table_long = non_phage_otu_table_long %>%
  mutate(software = case_when(
    software == "INSA" ~ "INSaFLU-TELEVIR",
    software == "GD" ~ "Genome Detective",
    software == "CZID" ~ "CZ.ID",
    software == "KR" ~ "Kraken2",
    TRUE ~ as.character(software)  # Keep the original value if no match
  ))



# Combine all unique Family values
all_family_ids <- unique(non_phage_otu_table_long$family); all_family_ids # old version

# Use a larger base palette from RColorBrewer for high contrast
color_palette_base <- RColorBrewer::brewer.pal(12, "Set3")  # "Set3" for a larger palette

# Manually add distinct colors
additional_colors <- c("#000000", "#FF00FF", "#00FFFF", "#FFFF00", "#FF0000", "#00FF00", "#0000FF", "#800080", "#808000", "#008080")

# Function to intercalate colors
intercalate_colors <- function(base_colors, additional_colors) {
  max_length <- max(length(base_colors), length(additional_colors))
  intercalated <- c()
  for (i in seq_len(max_length)) {
    if (i <= length(base_colors)) {
      intercalated <- c(intercalated, base_colors[i])
    }
    if (i <= length(additional_colors)) {
      intercalated <- c(intercalated, additional_colors[i])
    }
  }
  return(intercalated)
}

# Intercalate the base and additional colors
intercalated_palette <- intercalate_colors(color_palette_base, additional_colors)

# Ensure the palette has enough colors for all unique Taxonomy_ID values
num_colors_needed <- length(all_family_ids) # Old version
extended_palette_generator <- colorRampPalette(intercalated_palette)
color_palette_total <- extended_palette_generator(num_colors_needed)

non_phage_otu_table_long$family <- fct_reorder(non_phage_otu_table_long$family, non_phage_otu_table_long$Relative_Abundance, .desc=TRUE) # old version


# After reordering, get the reordered levels
reordered_family_levels <- levels(non_phage_otu_table_long$family)

# Create a named vector for color mapping
color_mapping <- setNames(color_palette_total, all_family_ids)


png(file = "./plot_sample_by_workflow_by_rep_number_relative_frequency_vertical_2.png", bg="white",width = 17, height = 15,  units = "in", res = 300)
ggplot(non_phage_otu_table_long, aes(x = sample_number, y = Relative_Abundance, fill = family)) +  
  geom_bar(position="fill", stat = "identity", show.legend = TRUE) +
  scale_fill_manual('Family', values = color_mapping) +
  labs(title = "",
       x = "Replicate",
       y = "Relative frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 20),  # Change this value to adjust size
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 12),) +  # Move legend to the right for vertical space
  guides(fill = guide_legend(ncol = 2)) +  # Set legend to a single column
  facet_grid(software~Step_Code, scales="free", space="free_x") +
  scale_x_discrete(drop = FALSE)
dev.off()



# Relative_Abundance - Horizontal
# png(file = "./plot_sample_by_workflow_by_rep_number_relative_frequency_horizontal.png", bg="white",width = 10, height = 20,  units = "in", res = 300)
# ggplot(non_phage_otu_table_long, aes(x = Relative_Abundance, y = sample_number, fill = family)) +
#   geom_bar(position="fill", stat = "identity", show.legend = TRUE) +
#   scale_fill_manual('Family', values = color_mapping) +
#   labs(title = "",
#        x = "Relative Frequency",
#        y = "Replicate") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), # Rotate x-axis labels for better fit
#         # legend.title = element_blank(),
#         legend.position = "right") +  # Move legend to the right for vertical space
#   guides(fill = guide_legend(ncol = 1)) +  # Set legend to a single column
#   facet_grid(software~Step_Code, scales="free", space="free_x")
# dev.off()



######################################################## Most abundant Families
# Summarize data to find the most abundant families
abundant_families <- non_phage_otu_table_long %>%
  group_by(family) %>%                                # Group by family
  summarise(family, Relative_Abundance, sample) %>%  # Sum relative abundance
  arrange(-Relative_Abundance)                           # Sort in descending order without desc()

# View the results
print(abundant_families)
length(unique(abundant_families$family)) # 59

summary_table_families <- abundant_families %>%
  group_by(family) %>%
  dplyr::summarize(
    max_abundance = max(Relative_Abundance)*100,
    min_abundance = min(Relative_Abundance)*100,
    sample_max = sample[which.max(Relative_Abundance)],
    sample_min = sample[which.min(Relative_Abundance)]
  ) %>%
  arrange(-max_abundance)

print(summary_table_families)
write.xlsx(list(Families_Abundance = summary_table_families), 
           file = "most_abundant_families_output.xlsx")



#### All relative abundances for each family, in each sample (UNGROUPED)
non_phage_otu_table_long[order(-non_phage_otu_table_long$Relative_Abundance), ] %>%
  head(10) %>%
  select(sample_code_number, software, family, abundance, Relative_Abundance)




######################################################## Most abundant Genus
# Summarize data to find the most abundant families
abundant_genus <- non_phage_otu_table_long %>%
  group_by(genus) %>%                                # Group by family
  summarise(genus, Relative_Abundance, sample) %>%  # Sum relative abundance
  arrange(-Relative_Abundance)                           # Sort in descending order without desc()

# View the results# View the resultsfamily
print(abundant_genus)
length(unique(abundant_genus$genus)) # 60

summary_table_genus <- abundant_genus %>%
  group_by(genus) %>%
  dplyr::summarize(
    max_abundance = max(Relative_Abundance),
    min_abundance = min(Relative_Abundance),
    sample_max = sample[which.max(Relative_Abundance)],
    sample_min = sample[which.min(Relative_Abundance)]
  ) %>%
  arrange(-max_abundance)

print(summary_table_genus)
write.xlsx(list(Genus_Abundance = summary_table_genus), 
           file = "most_abundant_genus_output.xlsx")


#### All relative abundances for each genus, in each sample (UNGROUPED)
non_phage_otu_table_long[order(-non_phage_otu_table_long$Relative_Abundance), ] %>%
  head(10) %>%
  select(sample_code_number, software, genus, abundance, Relative_Abundance)


# Combine all unique Family values
all_family_ids_log <- unique(non_phage_otu_table_long$family); all_family_ids

# Use a larger base palette from RColorBrewer for high contrast
color_palette_base_log <- RColorBrewer::brewer.pal(12, "Set3")  # "Set3" for a larger palette

# Manually add distinct colors
additional_colors_log <- c("#000000", "#FF00FF", "#00FFFF", "#FFFF00", "#FF0000", "#00FF00", "#0000FF", "#800080", "#808000", "#008080")

# Function to intercalate colors
intercalate_colors_log <- function(base_colors, additional_colors_log) {
  max_length <- max(length(base_colors), length(additional_colors_log))
  intercalated <- c()
  for (i in seq_len(max_length)) {
    if (i <= length(base_colors)) {
      intercalated <- c(intercalated, base_colors[i])
    }
    if (i <= length(additional_colors_log)) {
      intercalated <- c(intercalated, additional_colors_log[i])
    }
  }
  return(intercalated)
}

# Intercalate the base and additional colors
intercalated_palette_log <- intercalate_colors_log(color_palette_base_log, additional_colors_log)

# Ensure the palette has enough colors for all unique Taxonomy_ID values
num_colors_needed_log <- length(all_family_ids) # Old version
extended_palette_generator_log <- colorRampPalette(intercalated_palette_log)
color_palette_total_log <- extended_palette_generator(num_colors_needed_log)

non_phage_otu_table_long$family <- fct_reorder(non_phage_otu_table_long$family, log(non_phage_otu_table_long$abundance), .desc=TRUE)

# After reordering, get the reordered levels
reordered_family_levels <- levels(non_phage_otu_table_long$family)

# Create a named vector for color mapping
color_mapping_log <- setNames(color_palette_total_log, all_family_ids_log)


# Log e (absolute)
png(file = "./plot_sample_by_workflow_by_rep_number_relative_frequency_vertical_log.png", bg="white",width = 17, height = 15,  units = "in", res = 300)
ggplot(non_phage_otu_table_long, aes(x = sample_number, y = log(abundance), fill = fct_reorder(family, log(abundance), .desc=TRUE))) +
  geom_bar(position="stack", stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = color_mapping_log) +
  labs(title = "",
       x = "Replicate",
       y = "log e (Absolut abundance)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 20),  # Change this value to adjust size
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 12),) +
  guides(fill = guide_legend(ncol= 2)) +  # Set legend to a single row
  facet_grid(software~Step_Code, scales="free_x", space="free_x") +
  scale_x_discrete(drop = FALSE)  # Ensure that 1, 2, 3, and 5 always appear on the x-axis
dev.off()



png(file = "./plot_sample_by_workflow_by_rep_number_relative_frequency_vertical_log_2_2.png", bg="white",width = 17, height = 15,  units = "in", res = 300)
ggplot(non_phage_otu_table_long, aes(x = sample_number, y = log(abundance), fill = family)) +
  geom_bar(position="stack", stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = color_mapping_log) +
  labs(title = "",
       x = "Replicate",
       y = "log e (Absolut abundance)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 20),  # Change this value to adjust size
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 12),) +
  guides(fill = guide_legend(ncol= 2)) +  # Set legend to a single row
  facet_grid(software~Step_Code, scales="free_x", space="free_x") +
  scale_x_discrete(drop = FALSE)  # Ensure that 1, 2, 3, and 5 always appear on the x-axis
dev.off()


non_phage_otu_table_long$log_abundance <- log(non_phage_otu_table_long$abundance)




################################ ADD to final images ###########################
colnames(non_phage_otu_table_long)

head(non_phage_otu_table_long)

testing_facet_bug <- non_phage_otu_table_long %>%
  group_by(family, sample, TreatmentPlant, Step, software, Step_Code, sample_number) %>%
  dplyr::summarise(sum_abundance = sum(abundance))

testing_facet_bug$log_abundance = log(testing_facet_bug$sum_abundance)


png(file = "./plot_sample_by_workflow_by_rep_number_relative_frequency_vertical_log_2_4.png", bg="white",width = 17, height = 15,  units = "in", res = 300)
ggplot(testing_facet_bug, aes(x = sample_number, y = log_abundance, fill = family)) +
  geom_bar(position="stack", stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = color_mapping_log) +
  labs(title = "",
       x = "Replicate",
       y = "log e (Absolut abundance)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1), legend.title = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.title.x = element_text(size = 20),  # Change this value to adjust size
        axis.title.y = element_text(size = 20),
        strip.text = element_text(size = 12),) +
  guides(fill = guide_legend(ncol= 2)) +  # Set legend to a single row
  facet_grid(software~Step_Code, scales="free_x", space="free_x") +
  scale_x_discrete(drop = FALSE)  # Ensure that 1, 2, 3, and 5 always appear on the x-axis
dev.off()

################################################################################












# Number of unique families per fample
num_unique_families_per_sample = non_phage_otu_table_long %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(UniqueFamilies = n_distinct(family))

mean(num_unique_families_per_sample$UniqueFamilies)
max(num_unique_families_per_sample$UniqueFamilies)
min(num_unique_families_per_sample$UniqueFamilies)





##############################################
####### Pathogenic viruses identified ########
##############################################

nrow(non_phage_otu_table_long) # 855 # stayed the same

pathogenic_ids = non_phage_otu_table_long %>%
  filter(pathogenic_or_not == 1)

nrow(pathogenic_ids) #  182 stayed the same

colnames(pathogenic_ids)

dataframe_for_export = pathogenic_ids[c("Step_Code",
                                        "sample",
                                        "abundance",
                                        "family",
                                        "genus",
                                        "species",
                                        "pathogenic_or_not",
                                        "graph_path",
                                        "host_source")]

pathogenic_ids

write.xlsx(list(Potentially_Pathogenic_IDs = pathogenic_ids), 
           file = "potentially_pathogenic_output.xlsx")


# Summarize data to find the most abundant families
abundant_pathogenic_families <- pathogenic_ids %>%
  group_by(family) %>%                                # Group by family
  summarise(family, Relative_Abundance, sample) %>%  # Sum relative abundance
  arrange(-Relative_Abundance)                           # Sort in descending order without desc()

# View the results
print(abundant_pathogenic_families)

# Convert 'family' to character and then sort
abundant_pathogenic_families = sort(as.character(unique(abundant_pathogenic_families$family)))

unique(abundant_pathogenic_families) # 9


# Summarize data to find the most abundant families
abundant_pathogenic_genus <- pathogenic_ids %>%
  group_by(genus) %>%                                # Group by family
  summarise(genus, Relative_Abundance, sample) %>%  # Sum relative abundance
  arrange(-Relative_Abundance)                           # Sort in descending order without desc()

abundant_pathogenic_genus = sort(as.character(unique(abundant_pathogenic_genus$genus)))

# View the results
print(abundant_pathogenic_genus)
unique(abundant_pathogenic_genus)


colnames(pathogenic_ids)

pathogenic_ids

# Families and genera for the potentially pathogenic viruses identified
pathogenic_families_and_genera = pathogenic_ids %>%
  group_by(family, genus) %>%                         
  summarise(family, genus) %>%
  distinct() %>%
  arrange(family, genus)

# Convert 'family' and 'genus' columns to character
pathogenic_families_and_genera$family <- as.character(pathogenic_families_and_genera$family)
pathogenic_families_and_genera$genus <- as.character(pathogenic_families_and_genera$genus)

# Display the updated data frame
pathogenic_families_and_genera


sorted_data <- pathogenic_families_and_genera %>%
  arrange(family, genus)

sorted_data


write.xlsx(list(Potentially_Pathogenic_Tax = sorted_data), 
           file = "potentially_pathogenic_taxonomies_output.xlsx")



##############################################
######## Preparing Phyloseq object ###########
##############################################

# Checking the column names to remove any unwanted columns
colnames(non_phage_otu_table)

# Remove X! and OTU_ID to refactor
non_phage_otu_table <- non_phage_otu_table[,!names(non_phage_otu_table) %in% c("X1", "OTU_ID")]
colnames(non_phage_otu_table)
tables_samples <- tables_samples[,!names(tables_samples) %in% c("X1")]

# Refactor OTU_ID so its correct
rownames(non_phage_otu_table) = paste0("OTU", 1:nrow(non_phage_otu_table)); rownames(non_phage_otu_table)
colnames(non_phage_otu_table)

# Load the OTU table that you read from a file, called otumat
otu_table = as.matrix(non_phage_otu_table[,c(9:131)])

# taxonomy table
tax_table = as.matrix(non_phage_otu_table[,c(1:8)])

# Check if the same samples are kept in both matrices
unique(colnames(otu_table))
length(unique(colnames(otu_table))) # 111
       
unique(tables_samples$SampleID)
length(unique(tables_samples$SampleID)) # 126

# If they are not, we keep only the common ones
tables_samples_filtered <- tables_samples[tables_samples$SampleID %in% colnames(otu_table), ]

# Check again
unique(colnames(otu_table) == tables_samples_filtered$SampleID) # TRUE, proceed


# Convert to phyloseq objects
otu <- otu_table(otu_table, taxa_are_rows = TRUE)
tax <- tax_table(as.matrix(tax_table))
samples <- sample_data(tables_samples_filtered);
# We need to make the SampleID the rownames to match with the otu_table's rownames
rownames(samples) <- samples$SampleID


# Create phyloseq object
physeq <- phyloseq(otu, tax, samples)
physeq


sample_data(physeq)


"
####### 1. Pre-processing
b. Quality Control
"
# Already done previously, min was 10
taxa_sums(physeq)
min(taxa_sums(physeq)) # 1, meaning the taxa with the least reads has 1

sample_sums(physeq)
min(sample_sums(physeq)) # 2, meaning the sample with the least reads has 2

"
####### 2. Data Transformation
# a. Normalization
"
# Not to be done due to various papers advising not to.


###################
"
####### 3. RAREFACTION CURVES
"
# Rarefaction curves are created by randomly re-sampling the pool of N samples multiple times and then plotting the average number of species found in each sample 

tab <- otu_table(physeq)

# Separar as otu tables, uma por software para fazer a curva de rarefação (com a mesma escala do y)
# Filter columns based on sample names
otu_table_czid <- data.frame(otu_table) %>% dplyr::select(contains("CZID"))
otu_table_gd <- data.frame(otu_table) %>% dplyr::select(contains("GD"))
otu_table_kr <- data.frame(otu_table) %>% dplyr::select(contains("KR"))
otu_table_insa <- data.frame(otu_table) %>% dplyr::select(contains("INSA"))

otu_czid <- otu_table(otu_table_czid, taxa_are_rows = TRUE)
otu_gd <- otu_table(otu_table_gd, taxa_are_rows = TRUE)
otu_kr <- otu_table(otu_table_kr, taxa_are_rows = TRUE)
otu_insa <- otu_table(otu_table_insa, taxa_are_rows = TRUE)

# Filter rows based on the "Software" column
samples_data_czid <- tables_samples_filtered %>% dplyr::filter(Software == "CZ.ID")
sample_data_gd <- tables_samples_filtered %>% dplyr::filter(Software == "Genome Detective")
sample_data_kr <- tables_samples_filtered %>% dplyr::filter(Software == "Kraken2")
sample_data_insa <- tables_samples_filtered %>% dplyr::filter(Software == "INSaFLU-TELEVIR")

samples_czid <- sample_data(samples_data_czid)
samples_gd <- sample_data(sample_data_gd)
samples_kr <- sample_data(sample_data_kr)
samples_insa <- sample_data(sample_data_insa)

# SampleID must be rownames, otherwise it does not match
rownames(samples_czid) <- samples_czid$SampleID
rownames(samples_gd) <- samples_gd$SampleID
rownames(samples_kr) <- samples_kr$SampleID
rownames(samples_insa) <- samples_insa$SampleID

# Create phyloseq object
physeq_czid <- phyloseq(otu_czid, samples_czid)
physeq_gd <- phyloseq(otu_gd, samples_gd)
physeq_kr <- phyloseq(otu_kr, samples_kr)
physeq_insa <- phyloseq(otu_insa, samples_insa)




png(file = "./combined_rarefaction_curves.png", bg="white", width = 12, height = 10, units = "in", res = 600)
par(mfrow = c(2, 2))  # 2x2 layout for four plots
par(mar = c(5, 6, 4, 2) + 0.1)  # Set margins
par(cex.axis = 1.5, cex.lab = 1.5)  # Set axis and label sizes

# Plot 1: CZID
tab_czid <- otu_table(physeq_czid)
class(tab_czid) <- "matrix"
vegan::rarecurve(t(tab_czid), step = 1000, lwd = 2, ylab = "OTU", cex = 1.0, main = "CZ.ID") # ylim = ylim_values
# axis(2, at = yticks, labels = yticks)

tab_gd <- otu_table(physeq_gd)
class(tab_gd) <- "matrix"
vegan::rarecurve(t(tab_gd), step = 1000, lwd = 2, ylab = "OTU",cex = 1.0, main = "Genome Detective")
# axis(2, at=yticks, labels=yticks)

tab_kr <- otu_table(physeq_kr)
class(tab_kr) <- "matrix"
vegan::rarecurve(t(tab_kr),step = 1000, lwd = 2, ylab="OTU",cex = 1.0, main = "Kraken2")
# axis(2, at=yticks, labels=yticks)

tab_insa <- otu_table(physeq_insa)
class(tab_insa) <- "matrix"
vegan::rarecurve(t(tab_insa), step = 1000, lwd = 2, ylab = "OTU", cex = 1.0, ylim = c(0, 500), main = "INSaFLU-TELEVIR")
yticks <- seq(0, 500, by = 200)
axis(2, at = yticks, labels = yticks)

# Close the PNG device
dev.off()


# # Set the y-axis limits and labels
# ylim_values <- c(0, 200)
# yticks <- c(0, 50, 100, 150, 200)
# 
# tab_czid <- otu_table(physeq_czid)
# class(tab_czid) <- "matrix"
# png(file = "./plot_rarefaction_curve_czid.png", bg="white",width = 12, height = 10,  units = "in", res = 300)
# par(mar = c(5, 6, 4, 2) + 0.1)  # Increasing the second value (left margin)
# par(cex.axis = 1.5, cex.lab = 1.5) # Adjust size of axis labels and title 
# vegan::rarecurve(t(tab_czid), step = 1000, lwd = 2, ylab = "OTU", cex = 1.0, ylim = ylim_values)
# axis(2, at = yticks, labels = yticks) # Add axis ticks, ensuring tick labels are appropriately sized
# dev.off()


# tab_gd <- otu_table(physeq_gd)
# class(tab_gd) <- "matrix"
# png(file = "./plot_rarefaction_curve_gd.png", bg="white",width = 12, height = 10,  units = "in", res = 300)
# par(mar = c(5, 6, 4, 2) + 0.1)  # Increasing the second value (left margin)
# par(cex.axis = 1.5, cex.lab = 1.5) # Adjust size of axis labels and title 
# vegan::rarecurve(t(tab_gd), step = 1000, lwd = 2, ylab = "OTU",cex = 1.0, ylim = ylim_values)
# axis(2, at=yticks, labels=yticks)
# dev.off()
# 
# 
# tab_kr <- otu_table(physeq_kr)
# class(tab_kr) <- "matrix"
# png(file = "./plot_rarefaction_curve_kr.png", bg="white",width = 12, height = 10,  units = "in", res = 300)
# par(mar = c(5, 6, 4, 2) + 0.1)  # Increasing the second value (left margin)
# par(cex.axis = 1.5, cex.lab = 1.5) # Adjust size of axis labels and title vegan::rarecurve(t(tab_kr),step = 1000, lwd = 2, ylab = "OTU",cex = 0.7, ylim = ylim_values)
# vegan::rarecurve(t(tab_kr),step = 1000, lwd = 2, ylab="OTU",cex = 1.0, ylim = ylim_values)
# axis(2, at=yticks, labels=yticks)
# dev.off()
# 
# 
# tab_insa <- otu_table(physeq_insa)
# class(tab_insa) <- "matrix"
# png(file = "./plot_rarefaction_curve_insa.png", bg="white", width = 12, height = 10, units = "in", res = 300)
# par(mar = c(5, 6, 4, 2) + 0.1)  # Increasing the second value (left margin)
# par(cex.axis = 1.5, cex.lab = 1.5) # Adjust size of axis labels and title 
# vegan::rarecurve(t(tab_insa), step = 1000, lwd = 2, ylab = "OTU", cex = 1.0, ylim = c(0, 500))
# yticks <- seq(0, 500, by = 200)
# axis(2, at = yticks, labels = yticks)
# dev.off()





"
####### 3. Abundance Estimation
a. Alpha Diversity
"
sample_data(physeq)$Software <- as.factor(sample_data(physeq)$Software)
sample_data(physeq)$Treatment_Step <- as.factor(sample_data(physeq)$Treatment_Step)

# png(file = "./plot_richness_diversity_boxplot_taxon.png", bg="white",width = 10, height = 7,  units = "in", res = 300)
# plot_richness(physeq, x="Treatment_Step", color="Software", measures=c("Observed","Shannon","Simpson","InvSimpson")) + geom_boxplot()
# dev.off()

richness_estimate <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson"))

richness_df_for_excel = richness_estimate; richness_df_for_excel
richness_df_for_excel$Sample = rownames(richness_estimate)

# Reorder column by name manually
new_order = c("Sample", "Observed", "Shannon", "Simpson")
richness_df_for_excel <- richness_df_for_excel[, new_order]

# Write to excel table
write_xlsx(richness_df_for_excel, "richness_estimates.xlsx")


# Add a column for SampleID in richness_estimate based on rownames
richness_estimate$SampleID <- rownames(richness_estimate)

# Merge the two data frames based on multiple columns
columns_to_merge <- c("SampleID", "Replicate", "Software", "Treatment_Plant", "Step", "Treatment_Step")
merged_data = dplyr::left_join(richness_estimate, tables_samples_filtered, by = "SampleID")

colnames(merged_data)

# To change labels in facet
software_labs <- c("CZID", "KR","GD", "INSA")
names(software_labs) <- c("CZ.ID","Kraken2","GenomeDetective", "INSaFLU")

treatment_step_labs <- c("1 - 1", 
                         "1 - 2",
                         "1 - 3", 
                         "1 - 4",
                         "2 - 1", 
                         "2 - 2",
                         "2 - 3", 
                         "2 - 4")
names(treatment_step_labs) = c("WWTP 1 - Step  1",
                               "WWTP 1 - Step  2",
                               "WWTP 1 - Step  3",
                               "WWTP 1 - Step  4",
                               "WWTP 2 - Step  1",
                               "WWTP 2 - Step  2",
                               "WWTP 2 - Step  3",
                               "WWTP 2 - Step  4")
dot_size_in_plot = 3

colnames(merged_data)


png(file = "./plot_indexes_for_each_sample_software.png", bg="white",width = 10, height = 10,  units = "in", res = 300)
p1 <- ggplot(merged_data, aes(x= Software, y = Observed)) +
  geom_point(aes(col=as.character(Replicate)), size=dot_size_in_plot) +
  geom_line(aes()) +
  scale_colour_manual(name="Replicate", labels=c("1","2","3","5"),values = c("red","green","blue", "purple")) +
  theme_bw() +
  labs(title = "",
       x = "",
       y = "Observed index",tag="A") +
  facet_grid(~Treatment_Step, scales="free", space="free_x",labeller = labeller(Software = software_labs, Treatment_Step = treatment_step_labs)) +
  scale_y_continuous(limits = c(min(merged_data$Observed, na.rm = TRUE), max(merged_data$Observed, na.rm = TRUE))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme (axis.text.x = element_text(angle = 45, hjust = 1)) + # theme(axis.text.x = element_blank()) +  # OLD VERSION  
  theme(legend.position="none") 


p2 <- ggplot(merged_data, aes(x= Software, y = Shannon)) +
  geom_point(aes(col=as.character(Replicate)), size=dot_size_in_plot) +
  geom_line(aes()) +
  scale_colour_manual(name="Replicate", labels=c("1","2","3","5"),values = c("red","green","blue", "purple")) +
  theme_bw() +
  labs(title = "",
       x = "",
       y = "Shannon index",tag="B") +
  facet_grid(~Treatment_Step, scales="free", space="free_x",labeller = labeller(Software = software_labs, Treatment_Step = treatment_step_labs)) +
  scale_y_continuous(limits = c(min(merged_data$Shannon, na.rm = TRUE), max(merged_data$Shannon, na.rm = TRUE))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme (axis.text.x = element_text(angle = 45, hjust = 1)) + # theme(axis.text.x = element_blank()) +  # OLD VERSION
  theme(legend.position="none") 

p3 <- ggplot(merged_data, aes(x= Software, y = Simpson)) +
  geom_point(aes(col=as.character(Replicate)), size=dot_size_in_plot) +
  geom_line(aes()) +
  scale_colour_manual(name="Replicate", labels=c("1","2","3","5"),values = c("red","green","blue", "purple")) +
  theme_bw() +
  labs(title = "",
       x = "",
       y = "Simpson index",tag="C") +
  facet_grid(~Treatment_Step, scales="free", space="free_x",labeller = labeller(Software = software_labs, Treatment_Step = treatment_step_labs)) +
  scale_y_continuous(limits = c(min(merged_data$Simpson, na.rm = TRUE), max(merged_data$Simpson, na.rm = TRUE))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme (axis.text.x = element_text(angle = 45, hjust = 1)) + # theme(axis.text.x = element_blank()) +  # OLD VERSION  
  theme(legend.position="none") 


# Create grid
ggpubr::ggarrange(
  p1, p2, p3, # list of plots
  labels = "AUTO", # labels
  common.legend = TRUE, # COMMON LEGEND
  legend = "bottom", # legend position
  align = "hv", # Align them both, horizontal and vertical
  nrow = 3 # number of rows
)

dev.off()



# Test for normality
shapiro.test(merged_data$Observed) # p-value < 0.0001, not normal distribution
shapiro.test(merged_data$Shannon) # p-value < 0.0001, not normal distribution
shapiro.test(merged_data$Simpson) # p-value < 0.0001, not normal distribution
"
The Shapiro-Wilk test returned p-values less than 0.0001, 
indicating that none of the variables follow a normal distribution.
"

"
Levene's test to check for variance homogeneity across 
different software groups for each diversity index.
"
library(car)
leveneTest(merged_data$Observed ~ merged_data$Software) # p < 0.0004775
leveneTest(merged_data$Shannon ~ merged_data$Software) # p = 0.9292, homogenous
leveneTest(merged_data$Simpson ~ merged_data$Software) # p = 0.1693, homogenous
"
# normality is violated but variance is homogeneous:
- Shannon
- Simpson
"
kruskal.test(Observed ~ Software, data = merged_data) # p-value = 3.758e-08
kruskal.test(Shannon ~ Software, data = merged_data) # p-value = 9.194e-09
kruskal.test(Simpson ~ Software, data = merged_data) # p-value = 7.793e-12
"
Which indicates that there are significant differences in the 
two indexes between at least two of the software 
groups.

Follow up with post-hoc tests to see which specific groups are
different from one another.
"

dunn.test(merged_data$Observed, merged_data$Software, method = "holm")
"
Diff:
GD vs INSA
CZID vs KR
INSA vs KR
"


dunn.test(merged_data$Shannon, merged_data$Software, method = "holm")
"
Diff:
CZID vs INSA
INSA vs GD
INSA vs KR
"

dunn.test(merged_data$Simpson, merged_data$Software, method = "holm")
"
Diff:
INSA vs CZID
INSA vs GD
INSA vs KR
"


colnames(richness_estimate)


################################################################################
"Undefined"

undefined_data_grouped <- undefined %>% 
  dplyr::group_by(SampleID) %>% 
  dplyr::summarise(Sum_Reads = sum(Reads),
            .groups = 'drop'); undefined_data_grouped


min(undefined_data_grouped$Sum_Reads)
median(undefined_data_grouped$Sum_Reads)
max(undefined_data_grouped$Sum_Reads)




"
####### 3. Abundance Estimation
b. Beta Diversity
"
distance_matrix <- phyloseq::distance(physeq, method = "bray", weighted=F)
ordination <- ordinate(physeq, method="PCoA", distance=distance_matrix)

pcoa_loadings = ordination$values
write_xlsx(pcoa_loadings, "pcoa_loadings.xlsx")

pcoa_loadings
# 1.301727e-01 = 13.01 %
# 1.116493e-01 = 11.16 %

png(file = "./plot_ordination.png", bg="white",width = 10, height = 7,  units = "in", res = 300)
plot_ordination(physeq, ordination, color="Software") + theme(aspect.ratio=1)
dev.off()




# Principal Coordinates Analysis
pcoa_bray <- pcoa(distance_matrix, correction="none", rn=NULL)
pcoa_bray
# 1.301727e-01 = 13.01 %
# 1.116493e-01 = 11.16 %

pco_1_percent = round(pcoa_bray$values[1:2,2][1]*100, digits=2);pco_1_percent
pco_2_percent = round(pcoa_bray$values[1:2,2][2]*100, digits=2);pco_2_percent



pcoa_bray_df <- data.frame(pcoa_bray$vectors[,1:3])
colnames(pcoa_bray_df) <- c("PCo1", "PCo2","PCo3")

# Add a column for SampleID in pcoa_bray_df based on rownames
pcoa_bray_df$SampleID <- rownames(pcoa_bray_df)

# Merge the two data frames based on multiple columns
columns_to_merge <- c("SampleID", "Replicate", "Software", "Treatment_Plant", "Step", "Treatment_Step")
merged_data_pcoa = dplyr::left_join(pcoa_bray_df, tables_samples_filtered, by = "SampleID")

colnames(merged_data_pcoa)


png(file = "./plot_pcoa_bray_software_treatment_plant.png", bg="white",width = 8, height = 8,  units = "in", res = 300)
ggplot(merged_data_pcoa, aes(x = PCo1, y = PCo2)) + 
  geom_point(aes(shape=treatment_plant, color=Software), size = 4) +
  scale_shape_manual(values= c(16, 17, 18, 19, 20, 21, 22, 23)) +  # Adjust the shape values as needed for each treatment plant
  scale_colour_manual(values = c("blue", "red", "green", "yellow")) +  # Adjust the colors for each software
  labs(title = "",
       x = paste("PCo 1 (", pco_1_percent, "%)", sep = ""), 
       y = paste("PCo 2 (", pco_2_percent, "%)", sep = "")) +
  theme(legend.position = "top") +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  theme_minimal()
dev.off()




# Add a column for "Sample_Letter" and "Replicate_Or_Not" to check if reps of same sample group or not
merged_data_pcoa$Sample_Letter = substr(merged_data_pcoa$Sample, 1, 1)

# Add a new column to your data frame
merged_data_pcoa$Pool_Replicate <- ifelse(merged_data_pcoa$Replicate == 5, "Pool", "Replicate")


png(file = "./plot_pcoa_bray_sample_per_rep_or_pool.png", bg="white",width = 8, height = 8,  units = "in", res = 300)
ggplot(merged_data_pcoa, aes(x = PCo1, y = PCo2)) + 
  geom_point(aes(shape = Pool_Replicate, color = Sample_Letter), size = 4) +
  scale_shape_manual(values = c("Pool" = 16, "Replicate" = 17)) +  # Assign shapes
  # Optional: Customize colors for each Sample
  # scale_colour_manual(values = c("Sample1" = "blue", "Sample2" = "red", ...)) +
  labs(title = "",
       x = paste("PCo 1 (", pco_1_percent, "%)", sep = ""), 
       y = paste("PCo 2 (", pco_2_percent, "%)", sep = "")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")
dev.off()







png(file = "./plot_distances_side_by_side.png", bg="white",width = 16, height = 10,  units = "in", res = 300)

p1 = ggplot(merged_data_pcoa, aes(x = PCo1, y = PCo2)) + 
  geom_point(aes(shape=treatment_plant, color=Software), size = 4) +
  scale_shape_manual(values= c(16, 17, 18, 19, 20, 21, 22, 23)) +  # Adjust the shape values as needed for each treatment plant
  scale_colour_manual(values = c("blue", "red", "green", "yellow")) +  # Adjust the colors for each software
  labs(title = "",
       x = paste("PCo 1 (", pco_1_percent, "%)", sep = ""), 
       y = paste("PCo 2 (", pco_2_percent, "%)", sep = "")) +
  theme(legend.position = "bottom") +
  geom_hline(yintercept=0, linetype="dashed", color = "gray") +
  geom_vline(xintercept=0, linetype="dashed", color = "gray") +
  theme_minimal()
  
p2 = ggplot(merged_data_pcoa, aes(x = PCo1, y = PCo2)) + 
  geom_point(aes(shape = Pool_Replicate, color = Sample_Letter), size = 4) +
  scale_shape_manual(values = c("Pool" = 16, "Replicate" = 17)) +  # Assign shapes
  # Optional: Customize colors for each Sample
  # scale_colour_manual(values = c("Sample1" = "blue", "Sample2" = "red", ...)) +
  labs(title = "",
       x = paste("PCo 1 (", pco_1_percent, "%)", sep = ""), 
       y = paste("PCo 2 (", pco_2_percent, "%)", sep = "")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")


# Create grid
ggpubr::ggarrange(
  p1, p2, # list of plots
  labels = "AUTO", # labels
  common.legend = FALSE, # COMMON LEGEND
  legend = "bottom", # legend position
  align = "hv", # Align them both, horizontal and vertical
  ncol = 2 # number of rows
)
dev.off()








########################
# PERMANOVA
# tested with permutational multivariate analysis of variance using distance matrices (PERMANOVA) with the adonis() function from the vegan package
# group homogeneity of variances was tested by implementing the betadisper() function
# https://rstudio-pubs-static.s3.amazonaws.com/343284_cbadd2f3b7cd42f3aced2d3f42dc6d33.html#introduction

#Run PERMANOVA on distances Bray-Curtis
perm<-adonis2(distance_matrix ~ tables_samples_filtered$Software +
                tables_samples_filtered$treatment_plant  +
                tables_samples_filtered$collection_step ,  # sample_data(physeq)$Treatment_Step, 
              permutations = 1000)
perm







##############################################
################# Crassphage #################
##############################################
"
No filter, taken from all_samples_for_R.xlsx
"
# colnames(all_samples_for_R)
# length(unique(all_samples_for_R$SampleID)) # 126
# 
only_crass <- all_samples_for_R[all_samples_for_R$order %in% "Crassvirales", ]
# 
# sort(unique(only_crass$family));length(sort(unique(only_crass$family)))
# # 4 (5 with undefined)
# 
# sort(unique(only_crass$genus)); length(sort(unique(only_crass$genus)))
# # 27 (28 with undefined)
# 
# sort(unique(only_crass$species)); length(sort(unique(only_crass$species)))
# # 60 (62 with undefined and uncultured)
# 
# sort(unique(only_crass$SampleID)); length(sort(unique(only_crass$SampleID)))
# # 66 sample-software pairs
# 
# sort(unique(only_crass$Sample)); length(sort(unique(only_crass$Sample)))
# # 27 (there were 28 in total)
# 
# 
# "Correlation between CrAssphage and various other metrics, see where to fit this"
# # Correlation between crass sum count and number of pathogenic identifications
# colnames(tables_samples)
# 
# head(tables_samples)
# 
# 
# # Shapiro-Wilk test
# shapiro.test(tables_samples$Number_of_viral_pathogenic_IDs) # p-value = 5.306e-13
# shapiro.test(tables_samples$Number_of_crassphages_reads) # p-value < 2.2e-16
# shapiro.test(tables_samples$Number_of_poly_reads) # p-value < 2.2e-16
# shapiro.test(tables_samples$Number_of_adeno_reads) # p-value < 2.2e-16
# shapiro.test(tables_samples$Number_of_boca_reads) # p-value < 2.2e-16
# 
# 
# # Spearman correlation
# cor.test(tables_samples$Number_of_crassphages_reads, tables_samples$Number_of_viral_pathogenic_IDs, method = "spearman")
# "
# rho = 0.7104539
# p-value < 2.2e-16
# "
# cor.test(tables_samples$Number_of_crassphages_reads, tables_samples$Number_of_poly_reads, method = "spearman")
# "
# rho = p-value = 9.517e-07
# p-value 
# "
# cor.test(tables_samples$Number_of_crassphages_reads, tables_samples$Number_of_boca_reads, method = "spearman")
# "
# rho = 0.5024344
# p-value = 1.056e-08 
# "
# cor.test(tables_samples$Number_of_crassphages_reads, tables_samples$Number_of_adeno_reads, method = "spearman")
# "
# rho = 0.2637531
# p-value = 0.004397
# "


"
Filter for 10 Reads
"
filter_dataframe_by_col_numeric_value <- function(df, col, abundance_threshold = 10) {
  # Filter based on abundance threshold
  df_filtered <- df %>%
    filter({{col}} >= abundance_threshold)
  
  # Return the processed dataframe
  return(df_filtered)
}


only_crass_filtered <- filter_dataframe_by_col_numeric_value(only_crass, Reads ,abundance_threshold = 10)


sort(unique(only_crass_filtered$family));length(sort(unique(only_crass_filtered$family)))
# 4 (5 with undefined)

sort(unique(only_crass_filtered$genus)); length(sort(unique(only_crass_filtered$genus)))
# 21 (22 with undefined)

sort(unique(only_crass_filtered$species)); length(sort(unique(only_crass_filtered$species)))
# 38 (40 with undefined and uncultured)

sort(unique(only_crass_filtered$SampleID)); length(sort(unique(only_crass_filtered$SampleID)))
# 32 sample-software pairs

sort(unique(only_crass_filtered$Sample)); length(sort(unique(only_crass_filtered$Sample)))
# 12 (there were 28 in total)

# Families and genera for the potentially pathogenic viruses identified
crassphage_families_and_genera = only_crass_filtered %>%
  group_by(family, genus) %>%                         
  summarise(family, genus) %>%
  distinct() %>%
  arrange(family, genus)


crassphage_families_and_genera

crassphage_families_and_genera = crassphage_families_and_genera[!grepl("undefined",crassphage_families_and_genera$family),]

crassphage_families_and_genera


# Convert 'family' and 'genus' columns to character
crassphage_families_and_genera$family <- as.character(crassphage_families_and_genera$family)
crassphage_families_and_genera$genus <- as.character(crassphage_families_and_genera$genus)

# Display the updated data frame
crassphage_families_and_genera


sorted_data_crass <- crassphage_families_and_genera %>%
  arrange(family, genus)

sorted_data_crass




write.xlsx(list(Crassphage_Families_and_Genera = crassphage_families_and_genera), 
           file = "crassphage_families_and_genera_output.xlsx")


# Families and genera for the potentially pathogenic viruses identified
crassphage_full_tax = only_crass_filtered %>%
  group_by(family, genus, species) %>%                         
  summarise(family, genus, species) %>%
  distinct() %>%
  arrange(family, genus, species)

crassphage_full_tax

crassphage_full_tax = crassphage_full_tax[!grepl("undefined",crassphage_full_tax$family),]

crassphage_full_tax

write.xlsx(list(Crassphage_Full_tax = crassphage_full_tax), 
           file = "crassphage_full_taxonomy_output.xlsx")


"Correlation between CrAssphage and various other metrics, see where to fit this"
# Correlation between crass sum count and number of pathogenic identifications
colnames(tables_samples)

head(tables_samples)

tables_samples_crass_filtered = filter_dataframe_by_col_numeric_value(tables_samples, Number_of_crassphages_reads ,abundance_threshold = 10)


# Shapiro-Wilk test
shapiro.test(tables_samples_crass_filtered$Number_of_viral_pathogenic_IDs) # p-value = 1.634e-07
shapiro.test(tables_samples_crass_filtered$Number_of_crassphages_reads) # p-value = 7.493e-15
shapiro.test(tables_samples_crass_filtered$Number_of_poly_reads) # p-value = 1.211e-10
shapiro.test(tables_samples_crass_filtered$Number_of_adeno_reads) # p-value = 1.932e-10
shapiro.test(tables_samples_crass_filtered$Number_of_boca_reads) # p-value = 3.176e-12


# Spearman correlation
cor.test(tables_samples_crass_filtered$Number_of_crassphages_reads, tables_samples_crass_filtered$Number_of_viral_pathogenic_IDs, method = "spearman")
"
rho = 0.3920924  
p-value = p-value = 0.006415
"
cor.test(tables_samples_crass_filtered$Number_of_crassphages_reads, tables_samples_crass_filtered$Number_of_poly_reads, method = "spearman")
"
rho = 0.3377514 
p-value = 0.02024
"
cor.test(tables_samples_crass_filtered$Number_of_crassphages_reads, tables_samples_crass_filtered$Number_of_adeno_reads, method = "spearman")
"
rho = 0.3023269  
p-value = p-value = 0.03888
"
cor.test(tables_samples_crass_filtered$Number_of_crassphages_reads, tables_samples_crass_filtered$Number_of_boca_reads, method = "spearman")
"
rho = 0.1809979 
p-value = p-value = 0.2234
"




