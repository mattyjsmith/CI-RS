################################################################################
# Runs worker_fun(id) for every id in parallel (forked workers), but stops any worker that runs for more than
# timeout_secs, so that one stuck replicate cannot block the whole run (a plain mclapply waits for it forever).
# Returns the ids that did NOT finish properly (stopped, failed or died); they can simply be run again.
################################################################################

run_with_watchdog <- function(ids, worker_fun, n_workers = 4, timeout_secs = 1800, poll_secs = 2) {
  queue    <- ids
  active   <- list()
  not_done <- c()
  
  while (length(queue) > 0 || length(active) > 0) {
    
    # start new workers while there is room
    while (length(queue) > 0 && length(active) < n_workers) {
      id  <- queue[1]
      queue <- queue[-1]
      job <- parallel::mcparallel({
        tryCatch({ worker_fun(id); TRUE }, error = function(e) paste("error:", conditionMessage(e)))
      }, silent = TRUE)
      active[[length(active) + 1]] <- list(job = job, id = id, start = Sys.time())
    }
    
    Sys.sleep(poll_secs)
    
    # which workers have finished?
    got  <- parallel::mccollect(lapply(active, function(a) a$job), wait = FALSE)
    done <- if (is.null(got)) character(0) else names(got)
    
    still <- list()
    for (a in active) {
      pid <- as.character(a$job$pid)
      if (pid %in% done) {
        r <- got[[pid]]
        if (!isTRUE(r)) {
          message("Replicate ", a$id, " did not finish properly: ", paste(r, collapse = " "))
          not_done <- c(not_done, a$id)
        }
      } else if (as.numeric(difftime(Sys.time(), a$start, units = "secs")) > timeout_secs) {
        message("Replicate ", a$id, " ran for more than ", timeout_secs, " seconds and was stopped")
        tools::pskill(a$job$pid, tools::SIGKILL)
        invisible(parallel::mccollect(a$job, wait = TRUE))
        not_done <- c(not_done, a$id)
      } else {
        still[[length(still) + 1]] <- a
      }
    }
    active <- still
  }
  not_done
}
