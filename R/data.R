#' Data Processing for mfeargen
#'
#' Data processing functions are defined in 02_data_preparation.R
#' This file is kept for package structure compatibility.
#'
#' @name data-processing
#' @keywords internal
NULL

#' Experiment 1 Data from Yu et al. (2023)
#'
#' Fear generalization data from Experiment 1 of Yu et al. (2023), which used
#' a single CS+ design (1CS). Participants underwent a conditioning phase where
#' one stimulus (CS+) was paired with an aversive outcome (US), followed by a
#' generalization phase where similar stimuli were presented without reinforcement.
#'
#' @format A data frame with 7,520 observations and 17 variables:
#' \describe{
#'   \item{participant}{Participant identifier (numeric, 1-40)}
#'   \item{trials}{Trial number within participant}
#'   \item{dphy}{Physical distance from CS+ (Euclidean distance in pixel space)}
#'   \item{dper}{Perceptual distance from CS+ (measured through MDS)}
#'   \item{block_number}{Experimental block number}
#'   \item{trial_number}{Trial number within block}
#'   \item{US}{US presence indicator (0 = absent, 1 = present)}
#'   \item{startle_circle}{Startle probe during stimulus presentation (0/1)}
#'   \item{startle_iti}{Startle probe during inter-trial interval (0/1)}
#'   \item{per_size}{Perceptual size of stimulus (MDS-derived)}
#'   \item{us_expect}{US expectancy rating (0-10 scale, NA for ITI trials)}
#'   \item{total_trial_nr}{Total trial number across experiment}
#'   \item{stimulus}{Stimulus label (CS+ or generalization stimuli)}
#'   \item{stimulus_phy}{Physical stimulus identifier (S1-S10)}
#'   \item{cs_trials}{Indicator for CS+ trials (1) vs generalization trials (0)}
#'   \item{cs_phy}{Physical size of CS+ stimulus}
#'   \item{phy_size}{Physical size of presented stimulus}
#' }
#'
#' @details
#' The experiment used circle stimuli varying in size. The CS+ was reinforced
#' during conditioning, and generalization was measured using stimuli of different
#' sizes. Distance metrics include both physical (pixel-based) and perceptual
#' (MDS-derived) measures.
#'
#' NA values in \code{us_expect}, \code{dphy}, and \code{dper} occur during
#' inter-trial intervals (ITI) when no stimulus was presented.
#'
#' @source Yu, K., Tuerlinckx, F., Vanpaemel, W., & Zaman, J. (2023).
#'   A multi-fidelity computational approach to understand human fear learning.
#'   \emph{Communications Psychology}, \emph{1}(1), 5.
#'
#' @references
#' \doi{10.1038/s44271-023-00005-0}
#'
#' @examples
#' \dontrun{
#' # Load the data
#' data(exp1_data)
#'
#' # View structure
#' str(exp1_data)
#'
#' # Number of participants
#' length(unique(exp1_data$participant))
#'
#' # Prepare for analysis
#' # (Rename columns to match expected format)
#' exp1_data$US_expect <- exp1_data$us_expect
#' exp1_data$Per_size <- exp1_data$per_size
#' exp1_data$Phy_size <- exp1_data$phy_size
#'
#' model_data <- prepare_mfeargen_data(
#'   data = exp1_data,
#'   model_type = "1cs"
#' )
#' }
"exp1_data"

#' Experiment 2 Data from Yu et al. (2023)
#'
#' Fear generalization data from Experiment 2 of Yu et al. (2023), which used
#' a two CS+ design (2CS). Participants underwent conditioning with two distinct
#' CS+ stimuli: one close to the perceptual prototype (CS+p) and one close to
#' the physical prototype (CS+m). This design allows assessment of whether
#' generalization is based on physical or perceptual similarity.
#'
#' @format A data frame with 7,200 observations and 24 variables:
#' \describe{
#'   \item{participant}{Participant identifier (numeric, 1-40)}
#'   \item{trials}{Trial number within participant}
#'   \item{d_phy_m}{Physical distance from CS+m (more physically prototypical)}
#'   \item{d_phy_p}{Physical distance from CS+p (more perceptually prototypical)}
#'   \item{d_per_m}{Perceptual distance from CS+m}
#'   \item{d_per_p}{Perceptual distance from CS+p}
#'   \item{block_number}{Experimental block number}
#'   \item{trial_number}{Trial number within block}
#'   \item{US}{US presence indicator (0 = absent, 1 = present)}
#'   \item{startle_circle}{Startle probe during stimulus presentation (0/1)}
#'   \item{startle_iti}{Startle probe during inter-trial interval (0/1)}
#'   \item{per_size}{Perceptual size of stimulus (MDS-derived)}
#'   \item{us_expect}{US expectancy rating (0-10 scale, NA for ITI trials)}
#'   \item{total_trial_nr}{Total trial number across experiment}
#'   \item{group}{Experimental group assignment}
#'   \item{stimulus}{Stimulus label (CS+, CS-, or generalization stimuli)}
#'   \item{stimulus_phy}{Physical stimulus identifier (S1-S10)}
#'   \item{csp_trials}{Indicator for CS+p trials}
#'   \item{csm_trials}{Indicator for CS+m trials}
#'   \item{us_p}{US presence for CS+p (1 = present, 0 = absent, -1 = other trial)}
#'   \item{us_m}{US presence for CS+m (1 = present, 0 = absent, -1 = other trial)}
#'   \item{cs_phy_p}{Physical size of CS+p stimulus}
#'   \item{cs_phy_m}{Physical size of CS+m stimulus}
#'   \item{phy_size}{Physical size of presented stimulus}
#' }
#'
#' @details
#' This experiment used a 2CS design to distinguish between physical and perceptual
#' generalization patterns. Two CS+ stimuli were reinforced: CS+p (perceptually
#' prototypical) and CS+m (physically prototypical). Generalization stimuli varied
#' systematically in both physical and perceptual space.
#'
#' The dual distance metrics (\code{d_phy_m}/\code{d_phy_p} and
#' \code{d_per_m}/\code{d_per_p}) allow assessment of whether participants
#' generalize based on physical similarity, perceptual similarity, or both.
#'
#' NA values in \code{us_expect} and distance variables occur during inter-trial
#' intervals (ITI) when no stimulus was presented.
#'
#' @source Yu, K., Tuerlinckx, F., Vanpaemel, W., & Zaman, J. (2023).
#'   A multi-fidelity computational approach to understand human fear learning.
#'   \emph{Communications Psychology}, \emph{1}(1), 5.
#'
#' @references
#' \doi{10.1038/s44271-023-00005-0}
#'
#' @examples
#' \dontrun{
#' # Load the data
#' data(exp2_data)
#'
#' # View structure
#' str(exp2_data)
#'
#' # Number of participants
#' length(unique(exp2_data$participant))
#'
#' # Prepare for analysis
#' # (Rename columns to match expected format)
#' exp2_data$US_expect <- exp2_data$us_expect
#' exp2_data$Per_size <- exp2_data$per_size
#' exp2_data$Phy_size <- exp2_data$phy_size
#' exp2_data$USp <- exp2_data$us_p
#' exp2_data$USm <- exp2_data$us_m
#'
#' model_data <- prepare_mfeargen_data(
#'   data = exp2_data,
#'   model_type = "2cs"
#' )
#' }
"exp2_data"