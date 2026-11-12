# ============================================================================
# STEP 1: Load Required Libraries
# ============================================================================
# Install any missing packages: install.packages(c("tidyverse", "viridis", "ggthemes"))

library(tidyverse)
library(viridis)
library(ggthemes)

# ============================================================================
# STEP 2: Define File Paths
# ============================================================================
# Adjust these paths if your working directory structure is different
data_dir <- "."
mapping_file_fp <- file.path(data_dir, "20260313_Mali_PERCH_mapping_file.tsv")
ab2_fp <- file.path(data_dir, "PERCH_antibody_panel2.tsv")
auto_fp <- file.path(data_dir, "PERCHTotal_AutoantigenPanel.tsv")

# ============================================================================
# STEP 3: Define the Exact Upper and Lower Order Lists
# ============================================================================

upper_antibodies <- c(
  "U1-snRNP", "U1-snRNP68", "U1-snRNA", "U1-A", "VEGFB", "SRP54", 
  "SSA1", "SSB", "Scl-70", "RNP", "Ro60", "Ro52", "Ribp P2", 
  "Ribp P1", "Ribp P0", "Proteinase 3", "PM-Scl75", "PM-Scl100", 
  "PL-7", "PL-12"
)

lower_antibodies <- c(
  "PDC-E2", "PDC-E1alpha", "PDC-E1beta", "MIP-1alpha", "MDA5", 
  "La/SSB", "Ku", "Jo-1", "IL-1RA", "IL-6", "IL-5", "IL-33", 
  "IL-32", "IL-31", "IL-21", "IL-22", "IL-17A", "IL-13",  # Removed the duplicate "IL-17A" here
  "IL-12", "IL-10", "IL-4", "IL-1", "IFN-omega", "IFN-alpha2", 
  "IFN-lambda3", "IFN-lambda2", "IFN-lambda1", "IFN-gamma", "IFN-beta", 
  "IFN-epsilon", "IFN-alpha17", "IFN-alpha14", "IFN-alpha8", 
  "IFN-alpha4", "IFN-alpha16", "IFN-alpha13", "IFN-alpha7", "IFN-alpha", 
  "IL-3", "GM-CSF", "HGF", "PDGF", "CXCL5", "CXCL3", "CXCL8", 
  "CXCL11", "CXCL12", "CXCL4", "CXCL13", "CXCL6", "CXCL10", "CXCL1", 
  "CCL2", "CCL25", "CCL26", "CCL13", "CCL15", "CCL19", "CCL5", 
  "CCL17", "CCL4", "CCL11", "CCL1", "CCL3", "ACE2"
)

# Create a reference dataframe for ordering and faceting
antibody_order_df <- tibble(
  antibody = c(upper_antibodies, lower_antibodies),
  portion = c(rep("Upper", length(upper_antibodies)), rep("Lower", length(lower_antibodies)))
) %>%
  mutate(antibody_ordered = factor(antibody, levels = rev(antibody))) %>%
  mutate(portion = factor(portion, levels = c("Upper", "Lower")))

# ============================================================================
# STEP 4: Load and Merge Data (Replicating the Preamble Logic)
# ============================================================================

# 1. Load Mapping File (for the 'death' column)
mapping_data <- read.delim(mapping_file_fp, stringsAsFactors = F) %>%
  # Replicate the subject_id creation logic from the preamble
  mutate(subject_id = paste0("p_", tolower(ParticipantID))) %>%
  select(subject_id, death)

# 2. Load Panel 2 (specifically to get the 'Group' / study_group for each sample)
ab2_data <- read_delim(ab2_fp, show_col_types = FALSE) %>%
  rename(age_group = "Age Range", Sample = "Location") %>%
  select(Sample, age_group, Group) %>%
  unique()

# 3. Load Autoantigen Panel and Merge
auto_long <- read_delim(auto_fp, show_col_types = FALSE) %>%
  rename(age_group = "Age Range") %>%
  # Join Group from Panel 2
  left_join(ab2_data, by = c("Sample", "age_group")) %>%
  # Reshape to long format
  pivot_longer(cols = -c(Sample, age_group, Group), names_to = "antibody", values_to = "value_orj") %>%
  # Replicate subject_id creation (e.g. M00231-1B -> p_m00231)
  separate(Sample, into = c("subject_id", "extra"), sep = "-", remove = FALSE) %>%
  mutate(subject_id = paste0("p_", tolower(subject_id))) %>%
  # Recode Group into study_group factors matching the original script
  rename(study_group = Group) %>%
  mutate(study_group = factor(study_group)) %>%
  mutate(study_group = fct_recode(study_group, 
                                  `Control without RTI` = "without RTI",
                                  `Control with RTI` = "with RTI",
                                  `Severe case` = "Severe",
                                  `Very severe case` = "Very Severe")) %>%
  # Join mapping data to get the 'death' column
  left_join(mapping_data, by = "subject_id") %>%
  # Engineer the exact 'study_group_death' column used for the heatmap x-axis
  mutate(study_group_collapse_control = fct_collapse(study_group, Control = c("Control without RTI", "Control with RTI"))) %>%
  mutate(study_group_death = ifelse(death %in% "Yes", "Death", as.character(study_group_collapse_control))) %>%
  mutate(study_group_death = factor(study_group_death, levels = c("Control", "Severe case", "Very severe case", "Death"))) %>%
  # Standardize Age Group factor
  mutate(age_group = factor(age_group, levels = c("0-5 months", "6-11 months", "12-23 months", "24-59 months")))

# ============================================================================
# STEP 5: Filter, Scale, and Summarize Data
# ============================================================================

global_antibodies <- c(
  "Anti-human IgG Fc fragment specific", 
  "Anti-human IgG F (ab') fragment specific", 
  "Anti-human IgG (H+L)", 
  "Human IgG from serum"
)

summary_data <- auto_long %>%
  # Filter out global and BPI
  filter(!antibody %in% c(global_antibodies, "BPI")) %>%
  # Keep only the antibodies in our defined Upper/Lower lists
  filter(antibody %in% c(upper_antibodies, lower_antibodies)) %>%
  # Filter out samples with missing clinical data 
  filter(!is.na(study_group_death) & !is.na(age_group)) %>%
  
  # Calculate Z-score per antibody
  group_by(antibody) %>%
  mutate(value_z = as.numeric(scale(value_orj))) %>%
  ungroup() %>%
  
  # Calculate Mean Z-score per clinical group and age group
  group_by(antibody, study_group_death, age_group) %>%
  summarize(mean_z = mean(value_z, na.rm = TRUE), .groups = 'drop') %>%
  
  # Join the ordering dataframe
  left_join(antibody_order_df, by = "antibody") %>%
  filter(!is.na(portion))

# ============================================================================
# STEP 6: Plot the Heatmap
# ============================================================================

heatmap_plot <- ggplot(summary_data, aes(x = study_group_death, y = antibody_ordered, fill = mean_z)) +
  geom_tile() +
  facet_grid(portion ~ age_group, scales = "free_y", space = "free_y") +
  scale_fill_viridis(option = "magma", name = "Mean Z-score") +
  theme_clean() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom",
    legend.key.width = unit(1.2, "cm")
  ) +
  labs(
    x = "",
    y = ""
  )

# ============================================================================
# STEP 7: View and Save
# ============================================================================

print(heatmap_plot)

ggsave("PERCH_Autoantigen_Heatmap_Split.pdf", 
       plot = heatmap_plot, 
       height = 14, 
       width = 7, 
       useDingbats = FALSE)
