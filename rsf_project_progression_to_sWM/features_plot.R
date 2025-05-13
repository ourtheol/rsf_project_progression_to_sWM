library(ComplexHeatmap)
library(circlize)

# training cohort
d <- read.table("./data/training_cohort_initial64patients.tsv", header = TRUE)
title <- "Heatmap of the features used for training the Random Survival Forest model"
cohort <- "train_set"

# # validation cohort 
# d <- read.table("./data/validation_cohort.tsv", header = TRUE)
# title <- "Heatmap of the features in the validation cohort"
# cohort <- "validation_set"



clinical <- c("Months",	"Progression",	"Baseline_IgM",	"Hb",	"B2M",	"Albumin",	"BMI")

# Separate the data
# Matrix for the main heatmap (binary features)
binary_matrix <- d[, !(colnames(d) %in% clinical)] # Select all columns except the first one

# Data frame for row annotations (continuous features)
annotation_df <- d[, (colnames(d) %in% clinical)]


# Define colors for annotations (ComplexHeatmap syntax)

# Continuous IgM
col_fun_igm = colorRamp2(
  range(annotation_df$Baseline_IgM), # Calculate range automatically
  c("lightblue", "darkblue")         # Low and high colors
)

# Continuous Hb
col_fun_hb = colorRamp2(
  range(annotation_df$Hb), # Calculate range automatically
  c("darkorchid1", "darkorchid4")         # Low and high colors
)

# Continuous B2M
col_fun_b2m = colorRamp2(
  range(annotation_df$B2M), # Calculate range automatically
  c("darksalmon", "coral3")         # Low and high colors
)

# Continuous Albumin
col_fun_albumin = colorRamp2(
  range(annotation_df$Albumin), # Calculate range automatically
  c("aquamarine", "aquamarine4")         # Low and high colors
)

# Continuous BMI
col_fun_bmi = colorRamp2(
  range(annotation_df$BMI), # Calculate range automatically
  c("brown1","brown4")         # Low and high colors
)

# Continuous Months
col_fun_months = colorRamp2(
  range(annotation_df$Months), # Calculate range automatically
  c("yellow","darkorange3")         # Low and high colors
)

# # Discrete Group
col_group = c("0" = "grey", "1" = "red")

# Define the main heatmap colors (binary data)
col_binary = c("0" = "white", "1" = "black") # Name colors by the values


# Create HeatmapAnnotation objects
# Specify 'which = "row"' for row annotations

# Annotation for the LEFT side 
left_anno <- HeatmapAnnotation(
  df = data.frame(Baseline_IgM = annotation_df$Baseline_IgM,
                  Hb = annotation_df$Hb,
                  B2M = annotation_df$B2M,
                  Albumin = annotation_df$Albumin,
                  BMI = annotation_df$BMI
  ),
  col = list(Baseline_IgM = col_fun_igm,
             Hb = col_fun_hb,
             B2M = col_fun_b2m,
             Albumin = col_fun_albumin,
             BMI = col_fun_bmi
  ),
  which = "row",
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 8)
)

# Annotation for the RIGHT side
right_anno <- HeatmapAnnotation(
  df = data.frame(Progression = annotation_df$Progression,
                  Months = annotation_df$Months),
  col = list(Progression = col_group,
             Months = col_fun_months),
  which = "row",
  show_legend = TRUE,
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 8)
)


# Create the Heatmap object
# Note: Input matrix needs to be numeric for ComplexHeatmap coloring usually
binary_matrix_numeric <- as.matrix(binary_matrix) # Ensure it's a numeric matrix

ht <- Heatmap(
  binary_matrix_numeric,
  name = "Presence",               # Name for the main heatmap legend
  col = col_binary,                # Colors for the main heatmap
  
  # --- Assign annotations ---
  left_annotation = left_anno,   # Annotation(s) for the left
  right_annotation = right_anno, # Annotation(s) for the right
  
  # --- Appearance & Clustering ---
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_side = "left",         # Keep original row names on the left
  column_title = title, 
  row_title = NULL,                # Hide default row title if desired
  
  # --- Legend specific for binary data ---
  heatmap_legend_param = list(
    title = "Feature",
    at = c(0, 1),
    labels = c("Absent", "Present")
  )
  # show_row_names = TRUE, # default
  # show_column_names = TRUE, # default
)

# Draw the heatmap
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


pdf(paste0("./results/", cohort, ".pdf"), 20, 20)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
