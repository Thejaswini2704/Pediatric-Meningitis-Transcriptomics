rm(list = ls())
gc()

library(glmnet)
library(caret)
library(dplyr)

# ---------------------------------------------------------
# Load expression matrix
# ---------------------------------------------------------

x <- read.csv(
  "results/biomarkers/Bacterial_vs_Viral_Top100_log2CPM.csv",
  row.names = 1,
  check.names = FALSE
)

x <- as.matrix(x)

# Samples are columns in expression matrix
x <- t(x)

# ---------------------------------------------------------
# Sample groups
# ---------------------------------------------------------

y <- factor(
  c(
    rep("Viral", 5),
    rep("Healthy", 5),
    rep("Bacterial", 5)
  )
)

stopifnot(nrow(x) == length(y))

cat("Samples:", nrow(x), "\n")
cat("Features:", ncol(x), "\n")
print(table(y))

# ---------------------------------------------------------
# Stability selection
# ---------------------------------------------------------

set.seed(123)

n_bootstrap <- 100

selection_count <- setNames(
  numeric(ncol(x)),
  colnames(x)
)

successful_models <- 0

for (i in seq_len(n_bootstrap)) {

  # Stratified subsampling
  idx <- createDataPartition(
    y,
    p = 0.8,
    list = FALSE
  )

  x_sub <- x[idx, , drop = FALSE]
  y_sub <- droplevels(y[idx])

  # Make sure all 3 groups are represented
  if (length(unique(y_sub)) < 3) {
    next
  }

  fit <- tryCatch(
    {
      cv.glmnet(
        x = x_sub,
        y = y_sub,
        family = "multinomial",
        alpha = 0.5,
        type.measure = "class",
        nfolds = 3
      )
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(fit)) {
    next
  }

  successful_models <- successful_models + 1

  # Extract coefficients at lambda.min
  co <- coef(
    fit,
    s = "lambda.min"
  )

  selected <- unique(
    unlist(
      lapply(
        co,
        function(m) {

          vals <- as.matrix(m)

          ids <- rownames(vals)

          ids[
            ids != "(Intercept)" &
            rowSums(abs(vals)) != 0
          ]
        }
      )
    )
  )

  selected <- intersect(
    selected,
    names(selection_count)
  )

  selection_count[selected] <-
    selection_count[selected] + 1
}

# ---------------------------------------------------------
# Stability results
# ---------------------------------------------------------

cat(
  "\nSuccessful models:",
  successful_models,
  "/",
  n_bootstrap,
  "\n"
)

stability_results <- data.frame(
  isoform_id = names(selection_count),
  Selection_Count = as.numeric(selection_count),
  Frequency = as.numeric(selection_count) /
    successful_models * 100
)

stability_results <- stability_results %>%
  arrange(desc(Frequency))

# ---------------------------------------------------------
# Save all features
# ---------------------------------------------------------

write.csv(
  stability_results,
  "results/biomarkers/Bacterial_vs_Viral_Stability_Selection.csv",
  row.names = FALSE
)

# ---------------------------------------------------------
# Top 20
# ---------------------------------------------------------

top20 <- stability_results %>%
  slice_head(n = 20)

write.csv(
  top20,
  "results/biomarkers/Bacterial_vs_Viral_Stable_Top20.csv",
  row.names = FALSE
)

cat("\nTop 20 stable transcripts:\n")
print(top20)
