################################################################################
# TRUE values for each scenario, from the counterfactual outcomes of simulated patients (in chunks; progress is saved after every
# chunk, so if interrupted just run it again). Chunks of a different size from the current setting are ignored.
# FIRST run it with n_total <- 1000000 (about 3 minutes per scenario) and look at the printed cumulative incidences: the untreated cancer
# incidence should be near 35% (aggressive), 12% (medium) and 3% (indolent). If not, change lambda in 00_scenarios.R.
# Then run with n_total <- 10000000 (about 35 minutes per scenario); change chunk_size at the same time as n_total.
################################################################################

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
setwd(here::here("simulations", "cancer_type"))
rm(list = ls())
library(Epi); library(plyr); library(cmprsk)

source("00_scenarios.R")
source(here::here("R", "02_load_lifetable.R"))
source(here::here("R", "03_data_design_extended.R"))
source(here::here("R", "04_simulate_survival_extended.R"))

n_total    <- 10000000    # patients per scenario: 1000000 for the quick check, 10000000 for the final truth
chunk_size <- 1000000     # 100000 with n_total 1000000; 1000000 with n_total 10000000
horizon    <- 2
chunk_file <- "truth_chunks_cancertype.csv"

done <- if (file.exists(chunk_file)) read.csv(chunk_file, stringsAsFactors = FALSE) else NULL
n_chunks <- ceiling(n_total / chunk_size)

for (s in seq_len(nrow(scenarios))) {
  sc <- scenarios[s, ]
  for (cc in seq_len(n_chunks)) {
    if (!is.null(done) && any(done$scenario == sc$scenario & done$chunk == cc & done$lambda == sc$lambda & done$chunk_size == chunk_size)) next
    cat(sprintf("\nScenario %s, chunk %d of %d, started %s\n", sc$scenario, cc, n_chunks, format(Sys.time(), "%H:%M")))
    set.seed(999 + cc)
    invisible(capture.output(base_data <- cDataDesignOptim(
      n = chunk_size, cens.admin = 2015, ydiagmin = 2010, ydiagmax = 2010,
      betaagecr = -0.5, betasex = 0.5, betadep = -0.5, betastage = 0, age_shift = sc$age_shift)))
    invisible(capture.output(sim_data <- cdatasimulationT1WeibOptim.v2(
      mydat = base_data, dist.T1 = "Weibull", lambda.weib = sc$lambda, rho.weib = 0.7,
      BetaAge = 0.1, Betasex = -0.1, Betatrt = 0.5, Betadep = 0.4, Betastage = 0,
      BetaOther = 0, myexpmort = expratesASD)))
    cif0 <- cuminc(ftime = sim_data$TimeT0, fstatus = sim_data$CauseT0)   # everyone untreated
    cif1 <- cuminc(ftime = sim_data$TimeT1, fstatus = sim_data$CauseT1)   # everyone treated
    p1 <- timepoints(cif1, times = horizon)$est
    p0 <- timepoints(cif0, times = horizon)$est
    row <- data.frame(scenario = sc$scenario, lambda = sc$lambda, chunk = cc, chunk_size = chunk_size,
                      cif_c1 = p1[1, 1], cif_c0 = p0[1, 1], cif_o1 = p1[2, 1], cif_o0 = p0[2, 1])
    write.table(row, chunk_file, sep = ",", row.names = FALSE,
                col.names = !file.exists(chunk_file), append = file.exists(chunk_file))
    rm(base_data, sim_data, cif0, cif1); invisible(gc())
  }
}

# ---- Combine the chunks of the current size ----
ch <- read.csv(chunk_file, stringsAsFactors = FALSE)
out <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(s) {
  sc <- scenarios[s, ]
  d <- ch[ch$scenario == sc$scenario & ch$lambda == sc$lambda & ch$chunk_size == chunk_size, ]
  if (nrow(d) == 0) return(NULL)
  m <- colMeans(d[, c("cif_c1", "cif_c0", "cif_o1", "cif_o0")])
  data.frame(Scenario = sc$scenario, lambda = sc$lambda, patients = nrow(d) * chunk_size,
             RD = unname(m["cif_c1"] - m["cif_c0"]), logRR = unname(log(m["cif_c1"] / m["cif_c0"])),
             RD_comp = unname(m["cif_o1"] - m["cif_o0"]), logRR_comp = unname(log(m["cif_o1"] / m["cif_o0"])),
             cancer_CIF_treated = unname(m["cif_c1"]), cancer_CIF_untreated = unname(m["cif_c0"]),
             other_CIF_treated = unname(m["cif_o1"]), other_CIF_untreated = unname(m["cif_o0"]),
             cancer_share_of_deaths_untreated = unname(m["cif_c0"] / (m["cif_c0"] + m["cif_o0"])))
}))
write.csv(out, "true_values_cancertype.csv", row.names = FALSE)
print(out, digits = 4)
