################################################################################
####  example_gl_to_ly.R
####
####  Simulate the LY arm from the GL arm reference using only
####  LY arm summary statistics (no individual-level LY data in simulation).
####
####  Workflow:
####    1. Load example_new.csv and split by treatment arm
####    2. Define column types
####    3. Extract LY arm summary statistics (means, SDs, mins, maxs)
####       -- these are the only LY arm inputs to the simulation
####    4. Run run_simulation() using GL arm (reference) + LY summaries (target)
####    5. Diagnostic: compare simulated output vs LY target summaries
####       Optional post-hoc validation vs actual LY individual data
####
####  Key point: run_simulation() receives only GL arm individual data
####  (dat_ref = dat.c.gl) plus four named vectors of LY summary statistics.
####  dat.c.ly is NOT passed to run_simulation().
####
####  Author: Emma "Yao" Chen
####  R >= 4.1.2 | EmpiricalSim package required
################################################################################

library(EmpiricalSim)

#------------------------------------------------------------------------------
# Step 1: Load data and split by treatment arm
#------------------------------------------------------------------------------

dat.c    <- read.csv("example_new.csv", header = TRUE)
dat.c.gl <- subset(dat.c, TRT01P == "GL", select = -TRT01P)
dat.c.ly <- subset(dat.c, TRT01P == "LY", select = -TRT01P)

cat("Arm sizes:\n")
cat("  GL (reference):", nrow(dat.c.gl), "patients,", ncol(dat.c.gl), "variables\n")
cat("  LY (target):   ", nrow(dat.c.ly), "patients,", ncol(dat.c.ly), "variables\n\n")

#------------------------------------------------------------------------------
# Step 2: Define column types
# Column order matches CSV: AGE, SEX, HBA1CBL, BFSGMGDL, BTRGMGDL,
#   NHypoe_BL, THypoe_BL, HBA1C_wk4..TRG_wk52 (15 cols), NHypoe_wk12..THypoe_wk52
#------------------------------------------------------------------------------

type_vec <- c(
  "continuous",          # AGE
  "binary",              # SEX
  rep("continuous", 3),  # HBA1CBL, BFSGMGDL, BTRGMGDL
  "binary",              # NHypoe_BL
  "ordinal",             # THypoe_BL
  rep("continuous", 15), # HBA1C_wk4, HBA1C_wk12, HBA1C_wk26, HBA1C_wk39, HBA1C_wk52,
                         # FSG_wk4, FSG_wk12, FSG_wk26, FSG_wk39, FSG_wk52,
                         # TRG_wk4, TRG_wk12, TRG_wk26, TRG_wk39, TRG_wk52
  rep("ordinal", 6)      # NHypoe_wk12, NHypoe_wk26, NHypoe_wk52,
                         # THypoe_wk12, THypoe_wk26, THypoe_wk52
)
stopifnot(length(type_vec) == ncol(dat.c.gl))

#------------------------------------------------------------------------------
# Step 3: LY arm summary statistics
#
# extract_target_summaries() computes the four vectors that run_simulation()
# needs for the non-survival "range" scaling: means, SDs, mins, maxs.
#
# In a real-world scenario these values would come from a published table or
# a one-time computation on trial data. From this point onward, only these
# four named vectors are used — not individual LY patient records.
#------------------------------------------------------------------------------

summ_ly      <- extract_target_summaries(dat.c.ly, types = type_vec)
target_means <- summ_ly$means
target_sds   <- summ_ly$sds
target_mins  <- summ_ly$mins
target_maxs  <- summ_ly$maxs

cat("=== LY arm target summary statistics (inputs to simulation) ===\n")
cat("Means:\n");  print(round(target_means, 3))
cat("SDs:\n");    print(round(target_sds,   3))
cat("Mins:\n");   print(round(target_mins,  3))
cat("Maxs:\n");   print(round(target_maxs,  3))
cat("\n")

#------------------------------------------------------------------------------
# Step 4: Simulate LY arm from GL arm reference
#
# dat_ref        = dat.c.gl  -- GL arm individual data (correlation structure)
# target_means/sds/mins/maxs -- LY arm summaries only (no individual LY data)
#
# dat.c.ly is NOT used here.
#------------------------------------------------------------------------------

set.seed(42)
result <- run_simulation(
  dat_ref        = dat.c.gl,
  types          = type_vec,
  N_sim          = 5000,
  target_means   = target_means,
  target_sds     = target_sds,
  target_mins    = target_mins,
  target_maxs    = target_maxs,
  scaling_method = "range",
  verbose        = TRUE
)

sim_data <- result$sim_data

cat("\nSimulated data: ", nrow(sim_data), "rows x", ncol(sim_data), "columns\n\n")

#------------------------------------------------------------------------------
# Step 5a: Primary diagnostic — simulated vs LY target summaries
# This comparison uses only the target_means/sds vectors, NOT dat.c.ly.
#------------------------------------------------------------------------------

cat("=== Simulated vs LY target summary statistics ===\n")
diag <- data.frame(
  column         = colnames(dat.c.gl),
  type           = type_vec,
  LY_target_mean = round(target_means, 3),
  Sim_mean       = round(colMeans(sim_data), 3),
  mean_diff      = round(colMeans(sim_data) - target_means, 3),
  LY_target_sd   = round(target_sds, 3),
  Sim_sd         = round(apply(sim_data, 2, sd), 3),
  stringsAsFactors = FALSE
)
print(diag, row.names = FALSE)
cat("\n")

#------------------------------------------------------------------------------
# Step 5b: Post-hoc validation — simulated vs actual LY individual data
# (Uses dat.c.ly for reference only; this step is optional and purely
#  informational. The simulation above did not use dat.c.ly.)
#------------------------------------------------------------------------------

cat("=== Post-hoc validation: Simulated vs actual LY arm ===\n")
cat("(dat.c.ly used here for comparison only — NOT used in simulation)\n\n")

valid <- data.frame(
  column          = colnames(dat.c.gl),
  type            = type_vec,
  LY_actual_mean  = round(colMeans(dat.c.ly,              na.rm = TRUE), 3),
  Sim_mean        = round(colMeans(sim_data),                            3),
  LY_actual_sd    = round(apply(dat.c.ly, 2, sd,          na.rm = TRUE), 3),
  Sim_sd          = round(apply(sim_data, 2, sd),                        3),
  stringsAsFactors = FALSE
)
print(valid, row.names = FALSE)
cat("\n")

cat("=== Example complete ===\n")
cat("sim_data is available for downstream analysis.\n")
