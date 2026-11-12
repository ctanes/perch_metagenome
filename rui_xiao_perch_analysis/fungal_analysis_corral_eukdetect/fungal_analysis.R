# =====================================================================
# 0. Setup
# =====================================================================

library(tidyverse)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(VennDiagram)
library(gridExtra)
library(ggpubr)
library(rstatix)
library(DescTools)
library(DT)

# Set group colors
group_colors <- c(
  Control  = "#BA616C",
  Survived = "#2767B0",
  Died     = "#BEC100"
)

# =====================================================================
# 1. Read files
# =====================================================================

# Mapping file
mapping <- read_tsv(
  "20260313_Mali_PERCH_mapping_file.tsv",
  na = c("", "NA"),
  show_col_types = FALSE
)

# BUSCO taxonomy link
busco <- read_tsv(
  "busco_taxid_link.txt",
  col_names = c("taxon_name", "tax_id"),
  show_col_types = FALSE
) %>%
  mutate(tax_id = as.character(tax_id)) %>%
  distinct(tax_id, .keep_all = TRUE) %>%
  mutate(
    taxon_label = str_remove(taxon_name, "^Fungi-") %>%
      str_replace_all("_", " ")
  )

# CPM matrix
cpm <- read_tsv(
  "cpm.matrix.tsv",
  na = c("", "NA"),
  show_col_types = FALSE
)

# Rename first column to tax_id
colnames(cpm)[1] <- "tax_id"

cpm <- cpm %>%
  mutate(tax_id = as.character(tax_id))

# =====================================================================
# 2. Create grouping variable
# =====================================================================

grouping <- mapping %>%
  select(SampleID, study_group_collapse, study_group, death) %>%
  mutate(
    group = case_when(
      study_group_collapse == "Community control" ~ "Control",
      study_group == "Very severe case" & death == "No" ~ "Survived",
      study_group == "Very severe case" & death == "Yes" ~ "Died",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group, levels = c("Control", "Survived", "Died")))

# Check group counts
print(table(grouping$group, useNA = "ifany"))

# =====================================================================
# 3. Subset CPM matrix to relevant samples and convert blanks to 0
# =====================================================================

sample_keep <- intersect(colnames(cpm)[-1], grouping$SampleID)

cat("Samples in mapping groups:", length(unique(grouping$SampleID)), "\n")
cat("Samples in CPM matrix:", length(colnames(cpm)[-1]), "\n")
cat("Samples kept:", length(sample_keep), "\n")

cpm_sub <- cpm %>%
  select(tax_id, all_of(sample_keep)) %>%
  mutate(across(-tax_id, ~ as.numeric(replace_na(., 0))))

# Create long format with group and taxonomy names
cpm_long <- cpm_sub %>%
  pivot_longer(-tax_id, names_to = "SampleID", values_to = "CPM") %>%
  left_join(grouping %>% select(SampleID, group), by = "SampleID") %>%
  left_join(busco, by = "tax_id") %>%
  filter(!is.na(group))

# Check object
cat("Long CPM rows:", nrow(cpm_long), "\n")
cat("Taxa:", n_distinct(cpm_long$tax_id), "\n")
cat("Samples:", n_distinct(cpm_long$SampleID), "\n")

# =====================================================================
# 4. Analysis 1: Detection stacked bar chart
# =====================================================================

detection <- cpm_long %>%
  group_by(SampleID, group) %>%
  summarise(total_CPM = sum(CPM), .groups = "drop") %>%
  mutate(Detection = if_else(total_CPM > 0, "Detected", "Undetected"))

detection_counts <- detection %>%
  count(group, Detection) %>%
  complete(group, Detection, fill = list(n = 0)) %>%
  group_by(group) %>%
  mutate(pct = n / sum(n) * 100)

# Fill colors: detected = group color, undetected = grey
detection_fill_vals <- c(
  "Detected.Control"   = "#BA616C",
  "Detected.Survived"  = "#2767B0",
  "Detected.Died"      = "#BEC100",
  "Undetected.Control" = "grey70",
  "Undetected.Survived"= "grey70",
  "Undetected.Died"    = "grey70"
)

# Percentage stacked bar
p_detection_pct <- ggplot(
  detection_counts,
  aes(x = group, y = pct, fill = interaction(Detection, group))
) +
  geom_col(position = "stack", width = 0.7, color = "black") +
  geom_text(
    aes(label = paste0(round(pct, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values = detection_fill_vals) +
  labs(
    title = "Fungal detection status by group",
    x = NULL,
    y = "Percentage of samples",
    fill = "Detection"
  ) +
  theme_bw() +
  theme(legend.position = "right")

print(p_detection_pct)


# Absolute count stacked bar
p_detection_count <- ggplot(
  detection_counts,
  aes(x = group, y = n, fill = interaction(Detection, group))
) +
  geom_col(position = "stack", width = 0.7, color = "black") +
  geom_text(
    aes(label = n),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values = detection_fill_vals) +
  labs(
    title = "Fungal detection counts by group",
    x = NULL,
    y = "Number of samples",
    fill = "Detection"
  ) +
  theme_bw() +
  theme(legend.position = "right")

print(p_detection_count)


# Stats: 3 x 2 Fisher exact test
det_table <- detection %>%
  count(group, Detection) %>%
  pivot_wider(names_from = Detection, values_from = n, values_fill = 0)

det_matrix <- det_table %>%
  select(Detected, Undetected) %>%
  as.matrix()

rownames(det_matrix) <- det_table$group

cat("\nDetection 3x2 Fisher exact test:\n")
print(fisher.test(det_matrix))

# Pairwise Fisher tests
pairwise_fisher_detection <- function(g1, g2) {
  tab <- table(
    Detection = detection$Detection[detection$group %in% c(g1, g2)],
    Group     = detection$group[detection$group %in% c(g1, g2)]
  )
  ft <- fisher.test(tab)
  tibble(
    comparison = paste(g1, "vs", g2),
    p_value = ft$p.value
  )
}

pairwise_detection_stats <- bind_rows(
  pairwise_fisher_detection("Control", "Survived"),
  pairwise_fisher_detection("Control", "Died"),
  pairwise_fisher_detection("Survived", "Died")
) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

cat("\nPairwise Fisher tests for detection:\n")
print(pairwise_detection_stats)

datatable(pairwise_detection_stats)

# =====================================================================
# 5. Analysis 2: Venn diagram of shared/unique taxa
# =====================================================================

# Group-level taxa presence: at least one sample in group has CPM > 0
taxa_by_group <- cpm_long %>%
  filter(CPM > 0) %>%
  distinct(group, tax_id) %>%
  left_join(busco, by = "tax_id")

taxa_sets <- split(taxa_by_group$tax_id, taxa_by_group$group)

venn_plot <- venn.diagram(
  x = taxa_sets,
  category.names = c("Control", "Survived", "Died"),
  filename = NULL,
  output = TRUE,
  col = "black",
  fill = c("#BA616C", "#2767B0", "#BEC100"),
  alpha = 0.45,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("#BA616C", "#2767B0", "#BEC100"),
  main = "Shared and unique fungal taxa"
)

grid.newpage()
grid.draw(venn_plot)


# Create table of unique and shared taxa
A <- taxa_sets[["Control"]]
B <- taxa_sets[["Survived"]]
C <- taxa_sets[["Died"]]

only_A <- setdiff(A, union(B, C))
only_B <- setdiff(B, union(A, C))
only_C <- setdiff(C, union(A, B))
AB_only <- setdiff(intersect(A, B), C)
AC_only <- setdiff(intersect(A, C), B)
BC_only <- setdiff(intersect(B, C), A)
all3 <- intersect(intersect(A, B), C)

partition_list <- list(
  "Control only"             = only_A,
  "Survived only"            = only_B,
  "Died only"                = only_C,
  "Control & Survived only"  = AB_only,
  "Control & Died only"      = AC_only,
  "Survived & Died only"     = BC_only,
  "All three groups"         = all3
)

venn_table <- bind_rows(lapply(names(partition_list), function(cat) {
  ids <- partition_list[[cat]]
  if (length(ids) == 0) {
    tibble(category = cat, tax_id = NA_character_)
  } else {
    tibble(category = cat, tax_id = ids)
  }
})) %>%
  left_join(busco, by = "tax_id") %>%
  select(category, tax_id, taxon_label, everything())

print(head(venn_table, 20))

datatable(venn_table)

# =====================================================================
# 6. Analysis 3: Total CPM per sample / fungal burden
# =====================================================================

sample_burden <- cpm_long %>%
  group_by(SampleID, group) %>%
  summarise(total_CPM = sum(CPM), .groups = "drop")

# Stats
kruskal_burden <- sample_burden %>%
  kruskal_test(total_CPM ~ group)

dunn_burden <- sample_burden %>%
  dunn_test(total_CPM ~ group, p.adjust.method = "BH")

cat("\nKruskal-Wallis test for total CPM:\n")
print(kruskal_burden)

cat("\nDunn post-hoc tests for total CPM:\n")
print(dunn_burden)

# Box + violin plot
p_burden <- ggplot(sample_burden, aes(x = group, y = total_CPM)) +
  geom_violin(aes(fill = group), alpha = 0.3, trim = TRUE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(aes(color = group), width = 0.1, alpha = 0.45) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  stat_compare_means(
    comparisons = list(
      c("Control", "Survived"),
      c("Control", "Died"),
      c("Survived", "Died")
    ),
    method = "wilcox.test",
    label = "p.adj",
    p.adjust.method = "BH"
  ) +
  labs(
    title = "Total fungal CPM by group",
    x = NULL,
    y = "Total CPM per sample"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p_burden


# =====================================================================
# Section 7. Analysis 4: Comprehensive heatmap
# =====================================================================

# Load viridis for color scale
library(viridis)

# Prepare column annotation, ordering columns by group
annot_col <- grouping %>%
  filter(SampleID %in% colnames(cpm_sub)[-1]) %>%   # exclude tax_id column
  arrange(group) %>%
  column_to_rownames("SampleID") %>%
  select(group)

# Reorder CPM matrix columns to match annotation order
cpm_ordered <- cpm_sub %>%
  select(tax_id, all_of(rownames(annot_col))) %>%
  column_to_rownames("tax_id") %>%
  as.matrix()

# Log10 transform with pseudocount 1
log_cpm <- log10(cpm_ordered + 1)

# Create mapping from tax_id to unique fungal name
tax_labels <- busco %>%
  mutate(unique_label = make.unique(taxon_label)) %>%   # avoid duplicate row names
  select(tax_id, unique_label) %>%
  deframe()

# Replace row names (tax_id) with fungal names
row_names <- tax_labels[rownames(log_cpm)]
row_names[is.na(row_names)] <- rownames(log_cpm)[is.na(row_names)]  # fallback for missing IDs
rownames(log_cpm) <- row_names

# Annotation colors
annot_colors <- list(group = group_colors)

# Heatmap with viridis color palette
pheatmap(
  log_cpm,
  annotation_col = annot_col,
  annotation_colors = annot_colors,
  show_rownames = TRUE,          # show fungal names
  show_colnames = FALSE,         # hide sample names for clarity
  cluster_rows = TRUE,           # cluster taxa by abundance pattern
  cluster_cols = FALSE,          # do NOT cluster samples; keep grouped by group
  color = viridis::viridis(100), # viridis color scale
  main = "Log10(CPM+1) fungal taxa by sample",
  fontsize_row = 6,              # small row font so many taxa can be visible
  border_color = NA
)



# =====================================================================
# Analysis 5: Richness per sample (fixed)
# =====================================================================

richness <- cpm_long %>%
  group_by(SampleID, group) %>%
  summarise(richness = sum(CPM > 0), .groups = "drop")

rich_summary <- richness %>%
  group_by(group) %>%
  summarise(
    mean = mean(richness),
    sd = sd(richness),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  )

# Stats
kruskal_richness <- richness %>%
  kruskal_test(richness ~ group)

dunn_richness <- richness %>%
  dunn_test(richness ~ group, p.adjust.method = "BH")

cat("\nKruskal-Wallis test for richness:\n")
print(kruskal_richness)

cat("\nDunn post-hoc tests for richness:\n")
print(dunn_richness)

# Bar chart with mean +/- SE and raw points
p_richness <- ggplot(richness, aes(x = group, y = richness, fill = group)) +
  geom_col(
    data = rich_summary,
    aes(y = mean),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    data = rich_summary,
    aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),   # fixed: x = group added
    width = 0.2,
    inherit.aes = FALSE
  ) +
  geom_jitter(width = 0.1, alpha = 0.45, show.legend = FALSE) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(
    comparisons = list(
      c("Control", "Survived"),
      c("Control", "Died"),
      c("Survived", "Died")
    ),
    method = "wilcox.test",
    label = "p.adj",
    p.adjust.method = "BH"
  ) +
  labs(
    title = "Fungal taxa richness by group",
    x = NULL,
    y = "Number of detected taxa per sample"
  ) +
  theme_bw() +
  theme(legend.position = "none")

p_richness





# =====================================================================
# 9. Analysis 6: Pairwise differential abundance, top 15 by log2FC
# =====================================================================

# Helper function to compare two groups for every taxon
pairwise_de <- function(data, g1, g2, pseudocount = 1) {
  comp_data <- data %>%
    filter(group %in% c(g1, g2)) %>%
    mutate(group = factor(group, levels = c(g1, g2)))
  
  taxa <- unique(comp_data$tax_id)
  
  res <- bind_rows(lapply(taxa, function(t) {
    vals <- comp_data %>% filter(tax_id == t)
    x <- vals$CPM[vals$group == g1]
    y <- vals$CPM[vals$group == g2]
    
    if (length(x) == 0 | length(y) == 0) return(NULL)
    
    mean_x <- mean(x)
    mean_y <- mean(y)
    
    log2FC <- log2((mean_x + pseudocount) / (mean_y + pseudocount))
    
    p_value <- tryCatch(
      wilcox.test(x, y, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    
    prev_x <- mean(x > 0) * 100
    prev_y <- mean(y > 0) * 100
    
    tibble(
      tax_id = t,
      mean_g1 = mean_x,
      mean_g2 = mean_y,
      prevalence_g1 = prev_x,
      prevalence_g2 = prev_y,
      log2FC = log2FC,
      p_value = p_value
    )
  })) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(desc(abs(log2FC)))
  
  return(res)
}

# Compare all pairs
de_cs <- pairwise_de(cpm_long, "Control", "Survived")
de_cd <- pairwise_de(cpm_long, "Control", "Died")
de_sd <- pairwise_de(cpm_long, "Survived", "Died")

# Add taxonomy names
de_cs <- de_cs %>% left_join(busco, by = "tax_id") %>%
  select(tax_id, taxon_label, everything())
de_cd <- de_cd %>% left_join(busco, by = "tax_id") %>%
  select(tax_id, taxon_label, everything())
de_sd <- de_sd %>% left_join(busco, by = "tax_id") %>%
  select(tax_id, taxon_label, everything())

# Save all DE results
write_tsv(de_cs, "DE_Control_vs_Survived_all.tsv")
write_tsv(de_cd, "DE_Control_vs_Died_all.tsv")
write_tsv(de_sd, "DE_Survived_vs_Died_all.tsv")

# Top 15 by absolute log2FC
top15_cs <- de_cs %>% slice_max(abs(log2FC), n = 15)
top15_cd <- de_cd %>% slice_max(abs(log2FC), n = 15)
top15_sd <- de_sd %>% slice_max(abs(log2FC), n = 15)

cat("\nTop 15 taxa by log2FC: Control vs Survived\n")
print(top15_cs)

cat("\nTop 15 taxa by log2FC: Control vs Died\n")
print(top15_cd)

cat("\nTop 15 taxa by log2FC: Survived vs Died\n")
print(top15_sd)

write_tsv(top15_cs, "top15_Control_vs_Survived.tsv")
write_tsv(top15_cd, "top15_Control_vs_Died.tsv")
write_tsv(top15_sd, "top15_Survived_vs_Died.tsv")


# =====================================================================
# Analysis 6b: Pairwise differential abundance plots
# =====================================================================

# Add direction annotation: which group is higher?
add_direction <- function(de_table, g1, g2) {
  de_table %>%
    mutate(
      higher_in = if_else(log2FC > 0, g1, g2),
      direction_color = if_else(log2FC > 0, group_colors[g1], group_colors[g2])
    )
}

de_cs_plot <- add_direction(de_cs, "Control", "Survived")
de_cd_plot <- add_direction(de_cd, "Control", "Died")
de_sd_plot <- add_direction(de_sd, "Survived", "Died")

# ---------------------------------------------------------------------
# Function: top 20 differential taxa bar plot
# ---------------------------------------------------------------------

plot_top20_log2FC <- function(de_plot, title) {
  top20 <- de_plot %>%
    slice_max(abs(log2FC), n = 20) %>%
    mutate(
      taxon_label = if_else(is.na(taxon_label), paste0("TaxID ", tax_id), taxon_label),
      taxon_label = fct_reorder(taxon_label, log2FC)   # sort by log2FC
    )
  
  # Determine which group names appear for caption
  g_pos <- unique(de_plot$higher_in[de_plot$log2FC > 0])
  g_neg <- unique(de_plot$higher_in[de_plot$log2FC < 0])
  cap <- paste0("Positive = higher in ", paste(g_pos, collapse = "/"),
                "; Negative = higher in ", paste(g_neg, collapse = "/"))
  
  ggplot(top20, aes(x = taxon_label, y = log2FC, fill = direction_color)) +
    geom_col(color = "black", width = 0.7) +
    coord_flip() +
    scale_fill_identity() +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    labs(
      title = title,
      x = NULL,
      y = "log2 Fold Change",
      caption = cap
    ) +
    theme_bw() +
    theme(
      plot.caption = element_text(size = 9, face = "italic"),
      axis.text.y = element_text(size = 8)
    )
}

# Generate and save top20 plots
p_top_cs <- plot_top20_log2FC(de_cs_plot, "Top 20 differential taxa: Control vs Survived")
p_top_cd <- plot_top20_log2FC(de_cd_plot, "Top 20 differential taxa: Control vs Died")
p_top_sd <- plot_top20_log2FC(de_sd_plot, "Top 20 differential taxa: Survived vs Died")

p_top_cs
p_top_cd
p_top_sd


# =====================================================================
# Analysis 6b (revised): Volcano plots for pairwise DE
# =====================================================================

plot_volcano_revised <- function(de_table, g1, g2, title) {
  
  de_vol <- de_table %>%
    mutate(
      # Handle NA p-values as non-significant
      p_adj_safe = if_else(is.na(p_adj), 1, p_adj),
      neg_log10_padj = pmin(-log10(p_adj_safe), 5),   # cap at 5
      significant = p_adj_safe < 0.05,
      label = if_else(is.na(taxon_label), paste0("TaxID ", tax_id), taxon_label)
    )
  
  ggplot(de_vol, aes(x = log2FC, y = neg_log10_padj, color = significant)) +
    geom_point(size = 6, alpha = 0.8) +
    
    # Color scheme: red = significant, grey = not significant
    scale_color_manual(
      values = c("FALSE" = "grey70", "TRUE" = "red"),
      breaks = c(FALSE, TRUE),
      labels = c("Not significant", "Significant")
    ) +
    
    # Reference lines
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    
    # Label only significant points
    geom_text(
      data = filter(de_vol, significant),
      aes(label = label),
      size = 3,
      vjust = -0.5,
      hjust = 0.5,
      check_overlap = TRUE,        # avoid overlapping labels
      color = "black"
    ) +
    
    # Custom x-axis with direction labels
    scale_x_continuous(
      limits = c(-5, 5),
      breaks = c(-5, 0, 5),
      labels = c(paste0("Higher in ", g2), "0", paste0("Higher in ", g1)),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    # Fixed y-axis 0 to 5
    scale_y_continuous(
      limits = c(0, 5),
      breaks = seq(0, 5, 1),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    
    labs(
      title = title,
      x = NULL,
      y = "-log10(adjusted p-value)",
      color = NULL
    ) +
    
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold")
    )
}

p_vol_cs <- plot_volcano_revised(de_cs_plot, "Control", "Survived", 
                                     "Volcano: Control vs Survived")
p_vol_cd <- plot_volcano_revised(de_cd_plot, "Control", "Died", 
                                     "Volcano: Control vs Died")
p_vol_sd <- plot_volcano_revised(de_sd_plot, "Survived", "Died", 
                                     "Volcano: Survived vs Died")

p_vol_cs
p_vol_cd
p_vol_sd




# =====================================================================
# 10. Additional analyses
# =====================================================================

# ----------
# 10a. Shannon diversity index
# ----------

shannon <- cpm_long %>%
  group_by(SampleID, group) %>%
  summarise(
    total = sum(CPM),
    shannon = if_else(
      total == 0,
      0,
      -sum((CPM[CPM > 0] / total) * log(CPM[CPM > 0] / total))
    ),
    .groups = "drop"
  ) %>%
  select(SampleID, group, shannon)

kruskal_shannon <- shannon %>%
  kruskal_test(shannon ~ group)

dunn_shannon <- shannon %>%
  dunn_test(shannon ~ group, p.adjust.method = "BH")

cat("\nKruskal-Wallis test for Shannon diversity:\n")
print(kruskal_shannon)

cat("\nDunn post-hoc tests for Shannon diversity:\n")
print(dunn_shannon)

p_shannon <- ggplot(shannon, aes(x = group, y = shannon, fill = group)) +
  geom_violin(alpha = 0.3, trim = TRUE) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.45, show.legend = FALSE) +
  scale_fill_manual(values = group_colors) +
  stat_compare_means(
    comparisons = list(
      c("Control", "Survived"),
      c("Control", "Died"),
      c("Survived", "Died")
    ),
    method = "wilcox.test",
    label = "p.adj",
    p.adjust.method = "BH"
  ) +
  labs(
    title = "Shannon diversity by group",
    x = NULL,
    y = "Shannon index"
  ) +
  theme_bw()

p_shannon



# ----------
# 10b. PCA of log10(CPM+1)
# ----------

pca_result <- prcomp(t(log_cpm), center = TRUE, scale. = FALSE)

pca_scores <- pca_result$x %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(grouping, by = "SampleID")

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4, alpha = 0.65) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "PCA of log10(CPM+1) fungal profiles",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)"),
    color = "Group"
  ) +
  theme_bw()

p_pca



# =====================================================================
# Section 10c: Top 20 most abundant taxa per group (fixed)
# =====================================================================

# Identify top 20 taxa by mean CPM within each group
top_taxa <- cpm_long %>%
  group_by(group, tax_id, taxon_label) %>%
  summarise(mean_CPM = mean(CPM), .groups = "drop") %>%
  group_by(group) %>%
  slice_max(mean_CPM, n = 20) %>%
  ungroup() %>%
  mutate(
    taxon_label = if_else(
      is.na(taxon_label),
      paste0("TaxID ", tax_id),
      taxon_label
    )
  )

# Create a unique label per group to allow independent ordering within facets
top_taxa <- top_taxa %>%
  # Sort by group and descending mean_CPM to set factor levels correctly
  arrange(group, desc(mean_CPM)) %>%
  mutate(
    facet_label = paste(group, taxon_label, sep = "__"),
    facet_label = factor(facet_label, levels = unique(facet_label))
  )

# Plot with facets and independent ordering
p_top_taxa <- ggplot(
  top_taxa,
  aes(x = facet_label, y = mean_CPM, fill = group)
) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ group, scales = "free") +
  coord_flip() +
  # Remove the group prefix from x-axis labels
  scale_x_discrete(labels = function(x) gsub(".*__", "", x)) +
  labs(
    title = "Top 20 most abundant fungal taxa by group",
    x = NULL,
    y = "Mean CPM"
  ) +
  theme_bw()

p_top_taxa


