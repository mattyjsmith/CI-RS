################################################################################
# Merge the per-replicate results (results/<run_type>/sim_*_n<n>_rep<j>.rds) into one file per sample size,
# in the format read by 09a_performance_bias.R and 09_performance_measures.R.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "main"))

run_type <- "bias"     # "bias" or "coverage"
sample_sizes <- if (run_type == "bias") c(1000, 3000, 5000) else c(1000)

for (n in sample_sizes) {
  files <- list.files(file.path("results", run_type), pattern = sprintf("_n%d_rep[0-9]+\\.rds$", n), full.names = TRUE)
  if (length(files) == 0) next
  res <- do.call(rbind, lapply(files, readRDS))
  res$n <- n
  saveRDS(res, sprintf("merged_sim_results_%s_n%d.rds", run_type, n))
  cat(sprintf("%s, n = %d: %d rows\n", run_type, n, nrow(res)))
  print(table(res$Setting))
}
