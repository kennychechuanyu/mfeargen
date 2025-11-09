## Code to prepare experimental datasets from Yu et al. (2023)
## This script processes the raw reference data and creates package-ready datasets

# Load required packages
library(tools)

# Define base path
base_path <- "c:/Users/u0142671/OneDrive - KU Leuven/u0142671_Kenny_Yu/Research/pacakge - mfeargen/mfeargen"

# ============================================================================
# Experiment 1: Single CS+ Design
# ============================================================================

message("Processing Experiment 1 data...")

# Read raw data
exp1_raw <- readRDS(file.path(base_path, "reference/Data_s1.rds"))

# Clean and validate Experiment 1 data
exp1_data <- exp1_raw

# Rename columns to be more R-friendly (snake_case)
names(exp1_data)[names(exp1_data) == "block number"] <- "block_number"
names(exp1_data)[names(exp1_data) == "trial number"] <- "trial_number"
names(exp1_data)[names(exp1_data) == "Startle_Circle"] <- "startle_circle"
names(exp1_data)[names(exp1_data) == "Startle_ITI"] <- "startle_iti"
names(exp1_data)[names(exp1_data) == "Per_size"] <- "per_size"
names(exp1_data)[names(exp1_data) == "US_expect"] <- "us_expect"
names(exp1_data)[names(exp1_data) == "Total_trial_nr"] <- "total_trial_nr"
names(exp1_data)[names(exp1_data) == "CStrials"] <- "cs_trials"
names(exp1_data)[names(exp1_data) == "CS_phy"] <- "cs_phy"
names(exp1_data)[names(exp1_data) == "Phy_size"] <- "phy_size"

# Convert participant to numeric for consistency
exp1_data$participant <- as.numeric(exp1_data$participant)

# Add metadata as attributes
attr(exp1_data, "experiment") <- "Experiment 1"
attr(exp1_data, "design") <- "Single CS+ (1CS)"
attr(exp1_data, "n_participants") <- length(unique(exp1_data$participant))
attr(exp1_data, "n_trials") <- nrow(exp1_data)
attr(exp1_data, "source") <- "Yu et al. (2023)"
attr(exp1_data, "doi") <- "10.1038/s44271-023-00005-0"

# Basic validation (NA values allowed - they occur during ITI trials)
stopifnot(
  !any(is.na(exp1_data$participant)),
  all(exp1_data$us_expect[!is.na(exp1_data$us_expect)] >= 0),
  all(exp1_data$us_expect[!is.na(exp1_data$us_expect)] <= 10),  # Scale is 0-10
  all(exp1_data$dphy[!is.na(exp1_data$dphy)] >= 0),
  all(exp1_data$dper[!is.na(exp1_data$dper)] >= 0)
)

message(sprintf("  - %d participants", attr(exp1_data, "n_participants")))
message(sprintf("  - %d observations", attr(exp1_data, "n_trials")))

# ============================================================================
# Experiment 2: Two CS Design (CS+m and CS+p)
# ============================================================================

message("Processing Experiment 2 data...")

# Read raw data
exp2_raw <- readRDS(file.path(base_path, "reference/Data_s2.rds"))

# Clean and validate Experiment 2 data
exp2_data <- exp2_raw

# Rename columns to be more R-friendly (snake_case)
names(exp2_data)[names(exp2_data) == "block number"] <- "block_number"
names(exp2_data)[names(exp2_data) == "trial number"] <- "trial_number"
names(exp2_data)[names(exp2_data) == "Startle_Circle"] <- "startle_circle"
names(exp2_data)[names(exp2_data) == "Startle_ITI"] <- "startle_iti"
names(exp2_data)[names(exp2_data) == "Per_size"] <- "per_size"
names(exp2_data)[names(exp2_data) == "US_expect"] <- "us_expect"
names(exp2_data)[names(exp2_data) == "Total_trial_nr"] <- "total_trial_nr"
names(exp2_data)[names(exp2_data) == "CSptrials"] <- "csp_trials"
names(exp2_data)[names(exp2_data) == "CSmtrials"] <- "csm_trials"
names(exp2_data)[names(exp2_data) == "USp"] <- "us_p"
names(exp2_data)[names(exp2_data) == "USm"] <- "us_m"
names(exp2_data)[names(exp2_data) == "CS_phy_p"] <- "cs_phy_p"
names(exp2_data)[names(exp2_data) == "CS_phy_m"] <- "cs_phy_m"
names(exp2_data)[names(exp2_data) == "Phy_size"] <- "phy_size"

# Convert participant to numeric for consistency
exp2_data$participant <- as.numeric(exp2_data$participant)

# Add metadata as attributes
attr(exp2_data, "experiment") <- "Experiment 2"
attr(exp2_data, "design") <- "Two CS+ design (2CS)"
attr(exp2_data, "n_participants") <- length(unique(exp2_data$participant))
attr(exp2_data, "n_trials") <- nrow(exp2_data)
attr(exp2_data, "source") <- "Yu et al. (2023)"
attr(exp2_data, "doi") <- "10.1038/s44271-023-00005-0"

# Basic validation (NA values allowed - they occur during ITI trials)
stopifnot(
  !any(is.na(exp2_data$participant)),
  all(exp2_data$us_expect[!is.na(exp2_data$us_expect)] >= 0),
  all(exp2_data$us_expect[!is.na(exp2_data$us_expect)] <= 10),  # Scale is 0-10
  all(exp2_data$d_phy_m[!is.na(exp2_data$d_phy_m)] >= 0),
  all(exp2_data$d_phy_p[!is.na(exp2_data$d_phy_p)] >= 0),
  all(exp2_data$d_per_m[!is.na(exp2_data$d_per_m)] >= 0),
  all(exp2_data$d_per_p[!is.na(exp2_data$d_per_p)] >= 0)
)

message(sprintf("  - %d participants", attr(exp2_data, "n_participants")))
message(sprintf("  - %d observations", attr(exp2_data, "n_trials")))

# ============================================================================
# Save datasets
# ============================================================================

message("Saving datasets...")

# Save to data/ directory using usethis-style compression
usethis::use_data(exp1_data, overwrite = TRUE, compress = "bzip2")
usethis::use_data(exp2_data, overwrite = TRUE, compress = "bzip2")

message("Dataset preparation complete!")
message("Files created:")
message("  - data/exp1_data.rda")
message("  - data/exp2_data.rda")

# ============================================================================
# Summary statistics
# ============================================================================

message("\n=== Dataset Summary ===")
message("\nExperiment 1:")
message(sprintf("  Participants: %d", length(unique(exp1_data$participant))))
message(sprintf("  Observations: %d", nrow(exp1_data)))
message(sprintf("  Variables: %d", ncol(exp1_data)))
message(sprintf("  US expectancy range: [%d, %d]",
                min(exp1_data$us_expect), max(exp1_data$us_expect)))

message("\nExperiment 2:")
message(sprintf("  Participants: %d", length(unique(exp2_data$participant))))
message(sprintf("  Observations: %d", nrow(exp2_data)))
message(sprintf("  Variables: %d", ncol(exp2_data)))
message(sprintf("  US expectancy range: [%d, %d]",
                min(exp2_data$us_expect), max(exp2_data$us_expect)))
