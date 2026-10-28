# Reproduce all results of the paper. Run from the repository root:
#   Rscript run_all.R          # results from the saved simulation output (results/sim)
#   Rscript run_all.R --sim    # rerun the simulations first (about 15 hours on 16 cores)
# The simulation output (results/sim/*.rds) is not stored in the repository because
# of its size; without it, run with --sim.
run <- function(script, args = character()) {
  message("== ", paste(script, paste(args, collapse = " ")))
  status <- system2("Rscript", c(script, args))
  if (status != 0) stop(script, " failed")
}

if ("--sim" %in% commandArgs(trailingOnly = TRUE)) {
  run("sim/sim_PH.R")
  for (shape in c("bathub", "decreasing", "increasing")) run("sim/sim_TVE.R", shape)
  run("sim/sim_NLE.R")
}

run("analysis/res_PH.R")      # Tables 2-3, Appendix A
run("analysis/res_TVE.R")     # Figures 1-3, Appendix B
run("analysis/res_NLE.R")     # Figures 4-6, Table C1, Appendix C
run("case_study/case_study_age.R")      # Figure 7
run("case_study/case_study_surgery.R")  # Figure 8
