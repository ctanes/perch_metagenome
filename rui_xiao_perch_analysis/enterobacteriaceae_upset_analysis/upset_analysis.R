# script for making custom upset plot show enterobacter detected

suppressMessages({
  library(ComplexUpset)
  library(ggplot2)
  library(dplyr)
})

# ------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------
n_top <- 20   # number of top species to use as sets

# Group colors (as requested)
group_colors <- c(
  "Control"          = "#CC6677",
  "Severe case"      = "#88CCEE",
  "Very severe case" = "#0077BB",
  "Death"            = "#BEC100",
  "Unknown"          = "#bdc3c7"
)

# ------------------------------------------------------------
# STEP 1. Load the species presence matrix (species x samples)
# ------------------------------------------------------------
mat <- read.delim("enterobacteriaceae_presence_matrix_474.tsv",
                  row.names = 1, check.names = FALSE)

# ------------------------------------------------------------
# STEP 2. Load abundance proportions (per-sample unassigned_entero)
# ------------------------------------------------------------
ab <- read.delim("perch_kraken_npop_entero_proportions.tsv",
                 stringsAsFactors = FALSE)

# ------------------------------------------------------------
# STEP 3. Load metadata (group + death)
# ------------------------------------------------------------
meta <- read.delim("20260313_Mali_PERCH_mapping_file.tsv",
                   stringsAsFactors = FALSE)
meta <- meta[, c("SampleID", "study_group", "death")]

# ------------------------------------------------------------
# STEP 4. Helper to build one UpSet plot on a FILTERED matrix
# ------------------------------------------------------------
make_filtered_upset <- function(thresh, suffix, title) {

  # --- 4a. Identify qualifying sample IDs ---
  ids <- ab$SampleID[ab$unassigned_entero > thresh]

  # --- 4b. Filter the binary matrix to those samples (columns) ---
  sub <- mat[, intersect(colnames(mat), ids), drop = FALSE]
  cat(sprintf("[%s] threshold %.4f -> %d samples, %d species\n",
              title, thresh, ncol(sub), nrow(sub)))

  # --- 4c. Transpose to samples x species, ensure 0/1 ---
  df <- as.data.frame(t(sub))
  df[] <- lapply(df, as.integer)
  df$SampleID <- rownames(df)

  # --- 4d. Merge metadata ---
  df <- merge(df, meta, by = "SampleID", all.x = TRUE)

  # --- 4e. Group factor for stacked color bars ---
  df$grp_col <- case_when(
    is.na(df$study_group) ~ "Unknown",
    df$study_group %in% c("Control without RTI", "Control with RTI") ~ "Control",
    df$study_group == "Severe case"      ~ "Severe case",
    df$study_group == "Very severe case" ~ "Very severe case",
    TRUE ~ df$study_group
  )
  df$grp_col[!is.na(df$death) & df$death == "Yes"] <- "Death"
  df$grp_col <- factor(df$grp_col,
                       levels = c("Control", "Severe case",
                                  "Very severe case", "Death", "Unknown"))

  # --- 4f. Top-N species sets by prevalence in the FILTERED samples ---
  present <- sort(rowSums(sub), decreasing = TRUE)
  sets <- names(present)[present > 0][1:min(n_top, sum(present > 0))]

  # --- 4g. Sort sets by size descending (so set-size bars descend) ---
  set_sizes <- sapply(sets, function(s) sum(df[[s]]))
  sets <- sets[order(set_sizes, decreasing = TRUE)]

  # --- 4h. Build the UpSet plot ---
  #     base annotation: per-sample stacked bars with BLACK borders
  p <- upset(
    df, sets,
    base_annotations = list(
      'Intersection size' = ggplot(mapping = aes(fill = grp_col)) +
        geom_bar(stat = "count", position = "stack",
                 aes(group = SampleID),      # each sample = its own unit
                 colour = "black",           # black border on every sample
                 linewidth = 0.25) +
        ylab("Intersection size") +
        scale_fill_manual(values = group_colors, name = "Group") +
        theme(legend.position = "bottom")
    ),
    set_sizes = upset_set_size(),
    name = "Enterobacteriaceae species",
    width_ratio = 0.25,
    min_size = 1,
    sort_sets = "descending",              # largest set on the left
    sort_intersections = "ascending",      # least complex first
    sort_intersections_by = "degree"       # complexity = # species in set
  )

  # --- 4i. Save as pdf / png / svg ---
  ggsave(paste0("upset_filter_v2_", suffix, ".pdf"), p, width = 14, height = 9)
  ggsave(paste0("upset_filter_v2_", suffix, ".png"), p, width = 14, height = 9, dpi = 300)
  ggsave(paste0("upset_filter_v2_", suffix, ".svg"), p, width = 14, height = 9)

  cat(sprintf("[%s] Saved upset_filter_v2_%s.pdf / .png / .svg\n",
              title, suffix))
  return(p)
}

# ------------------------------------------------------------
# STEP 5. Produce both filtered plots
# ------------------------------------------------------------
# Plot 1: filter to samples with unassigned_entero > 1%
p1 <- make_filtered_upset(thresh = 0.01, suffix = "1pct",
                          title = ">1%")

# Plot 2: filter to samples with unassigned_entero > 0.01%
p2 <- make_filtered_upset(thresh = 0.0001, suffix = "001pct",
                          title = ">0.01%")

# ------------------------------------------------------------
# STEP 6. (Optional) view in RStudio
# ------------------------------------------------------------
# print(p1)
# print(p2)
cat("Done.\n")
