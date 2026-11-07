################################################################################
# TRUE values for each scenario, from the counterfactual outcomes of 10,000,000 simulated patients (in chunks).
# Progress is saved after every chunk (truth_chunks_othercause.csv); if interrupted, just run it again.
# Also prints the 2-year cumulative incidences, so you can see how common other-cause death is in each scenario.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "other_cause"))
rm(list = ls())
library(Epi)
library(plyr)
library(cmprsk)

source("00_scenarios.R")
source(here::here("R", "02_load_lifetable.R"))
source(here::here("R", "03_data_design_extended.R"))
source(here::here("R", "04_simulate_survival_extended.R"))

n_total    <- 1000000   # patients per scenario (use 1000000 for a quick test)
chunk_size <- 100000
horizon    <- 2          # years since diagnosis
# Chunks are recorded with their size (new file ..._v2.csv); chunks of another size are ignored, so changing n_total or chunk_size can never silently reuse old chunks.
chunk_file <- "truth_chunks_othercause_v2.csv"

done <- if (file.exists(chunk_file)) read.csv(chunk_file, stringsAsFactors = FALSE) else NULL
n_chunks <- ceiling(n_total / chunk_size)

for (s in seq_len(nrow(scenarios))) {
  sc <- scenarios[s, ]
  for (cc in seq_len(n_chunks)) {
    if (!is.null(done) && any(done$scenario == sc$scenario & done$chunk == cc & done$chunk_size == chunk_size)) next
    cat(sprintf("\nScenario %s, chunk %d of %d, started %s\n", sc$scenario, cc, n_chunks, format(Sys.time(), "%H:%M")))
    
    set.seed(999 + cc)
    invisible(capture.output(base_data <- cDataDesignOptim(
      n = chunk_size, cens.admin = 2015, ydiagmin = 2010, ydiagmax = 2010,
      betaagecr = -0.5, betasex = 0.5, betadep = -0.5, betastage = 0, age_shift = sc$age_shift)))
    invisible(capture.output(sim_data <- cdatasimulationT1WeibOptim.v2(
      mydat = base_data, dist.T1 = "Weibull", lambda.weib = 0.25, rho.weib = 0.7,
      BetaAge = 0.1, Betasex = -0.1, Betatrt = 0.5, Betadep = 0.4, Betastage = 0,
      BetaOther = log(sc$HR_other), myexpmort = expratesASD)))
    
    cif0 <- cuminc(ftime = sim_data$TimeT0, fstatus = sim_data$CauseT0)   # everyone untreated
    cif1 <- cuminc(ftime = sim_data$TimeT1, fstatus = sim_data$CauseT1)   # everyone treated
    p1 <- timepoints(cif1, times = horizon)$est
    p0 <- timepoints(cif0, times = horizon)$est
    
    row <- data.frame(scenario = sc$scenario, chunk = cc, chunk_size = chunk_size,
                      cif_c1 = p1[1, 1], cif_c0 = p0[1, 1], cif_o1 = p1[2, 1], cif_o0 = p0[2, 1])
    write.table(row, chunk_file, sep = ",", row.names = FALSE,
                col.names = !file.exists(chunk_file), append = file.exists(chunk_file))
    rm(base_data, sim_data, cif0, cif1); invisible(gc())
  }
}

# ---- Combine the chunks ----
ch <- read.csv(chunk_file, stringsAsFactors = FALSE)
ch <- ch[ch$chunk_size == chunk_size, ]   # only chunks of the current size
out <- do.call(rbind, lapply(scenarios$scenario, function(sn) {
  d <- ch[ch$scenario == sn, ]
  if (nrow(d) == 0) return(NULL)
  m <- colMeans(d[, c("cif_c1", "cif_c0", "cif_o1", "cif_o0")])
  data.frame(Scenario = sn, n_chunks = nrow(d),
             RD = unname(m["cif_c1"] - m["cif_c0"]), logRR = unname(log(m["cif_c1"] / m["cif_c0"])),
             RD_comp = unname(m["cif_o1"] - m["cif_o0"]), logRR_comp = unname(log(m["cif_o1"] / m["cif_o0"])),
             cancer_CIF_treated = unname(m["cif_c1"]), cancer_CIF_untreated = unname(m["cif_c0"]),
             other_CIF_treated = unname(m["cif_o1"]), other_CIF_untreated = unname(m["cif_o0"]))
}))
write.csv(out, "true_values_othercause.csv", row.names = FALSE)
print(out, digits = 4)
