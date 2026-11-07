################################################################################
# Title: Simulation of weighting approach for CI & RS
# Subtitle: TRUE values from the counterfactual outcomes of a very large simulated population
# Both counterfactual worlds (everyone treated / everyone untreated) are simulated for every patient, and the true
# 2-year cumulative incidences of cancer death and other-cause death are read from them directly.
# The population is simulated in chunks to limit memory use; progress is saved after each chunk, so if the script
# is interrupted just run it again.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "main"))
rm(list = ls())
library(Epi)
library(plyr)
library(cmprsk)

source(here::here("R", "02_load_lifetable.R"))
source(here::here("R", "03_data_design.R"))
source(here::here("R", "04_simulate_survival.R"))

n_total    <- 10000000   # patients in total (use 1000000 for a quick test)
chunk_size <- 1000000    # patients per chunk (reduce if memory is a problem)
horizon    <- 2          # years since diagnosis
# Chunks are recorded with their size (new file ..._v2.csv); chunks of another size are ignored, so changing n_total or chunk_size can never silently reuse old chunks.
chunk_file <- "truth_chunks_v2.csv"

done <- if (file.exists(chunk_file)) read.csv(chunk_file) else NULL
n_chunks <- ceiling(n_total / chunk_size)

for (cc in seq_len(n_chunks)) {
  if (!is.null(done) && any(done$chunk == cc & done$chunk_size == chunk_size)) next
  cat(sprintf("\nChunk %d of %d, started %s\n", cc, n_chunks, format(Sys.time(), "%H:%M")))
  
  set.seed(999 + cc)
  invisible(capture.output(base_data <- cDataDesignOptim(
    n = chunk_size, cens.admin = 2015, ydiagmin = 2010, ydiagmax = 2010,
    betaagecr = -0.5, betasex = 0.5, betadep = -0.5, betastage = 0
  )))
  invisible(capture.output(sim_data <- cdatasimulationT1WeibOptim.v2(
    mydat = base_data, dist.T1 = "Weibull", lambda.weib = 0.25, rho.weib = 0.7,
    BetaAge = 0.1, Betasex = -0.1, Betatrt = 0.5, Betadep = 0.4, Betastage = 0,
    myexpmort = expratesASD
  )))
  
  cif0 <- cuminc(ftime = sim_data$TimeT0, fstatus = sim_data$CauseT0)   # everyone untreated
  cif1 <- cuminc(ftime = sim_data$TimeT1, fstatus = sim_data$CauseT1)   # everyone treated
  p1 <- timepoints(cif1, times = horizon)$est
  p0 <- timepoints(cif0, times = horizon)$est
  
  row <- data.frame(chunk = cc, chunk_size = chunk_size,
                    cif_c1 = p1[1, 1], cif_c0 = p0[1, 1],    # cancer death: treated, untreated
                    cif_o1 = p1[2, 1], cif_o0 = p0[2, 1])    # other-cause death: treated, untreated
  write.table(row, chunk_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(chunk_file), append = file.exists(chunk_file))
  rm(base_data, sim_data, cif0, cif1); invisible(gc())
}

# ---- Combine the chunks ----
ch <- read.csv(chunk_file)
ch <- ch[ch$chunk_size == chunk_size, ]   # only chunks of the current size
m  <- colMeans(ch[, c("cif_c1", "cif_c0", "cif_o1", "cif_o0")])
true_values <- data.frame(
  RD         = unname(m["cif_c1"] - m["cif_c0"]),
  logRR      = unname(log(m["cif_c1"] / m["cif_c0"])),
  RD_comp    = unname(m["cif_o1"] - m["cif_o0"]),
  logRR_comp = unname(log(m["cif_o1"] / m["cif_o0"]))
)
print(true_values)
write.csv(true_values, "true_values.csv", row.names = FALSE)

# Monte Carlo error of the true risk differences (needs at least 2 chunks)
if (nrow(ch) > 1) {
  cat("\nMonte Carlo SE of the true RD (cancer):      ", sd(ch$cif_c1 - ch$cif_c0) / sqrt(nrow(ch)), "\n")
  cat("Monte Carlo SE of the true RD (other cause): ", sd(ch$cif_o1 - ch$cif_o0) / sqrt(nrow(ch)), "\n")
}
