#' Experiment 1 Data from Yu et al. (2023)
#'
#' Trial-level data from the single-CS+ experiment in Yu et al. (2023).
#' Ratings and distances use the original study's variable names and units.
#'
#' @format A data frame with 7,520 rows and 17 variables:
#' \describe{
#'   \item{participant}{Participant identifier, 1-40.}
#'   \item{trials}{Trial order within participant.}
#'   \item{dphy, dper}{Physical and perceptual distances from CS+.}
#'   \item{block_number, trial_number}{Experimental block and trial within block.}
#'   \item{US}{Observed aversive outcome, 0 or 1.}
#'   \item{startle_circle, startle_iti}{Startle probe indicators during stimulus
#'     presentation and the inter-trial interval.}
#'   \item{per_size}{Perceptual stimulus measurement.}
#'   \item{us_expect}{US expectancy rating on the original 0-10 scale.}
#'   \item{total_trial_nr}{Trial number across the experiment.}
#'   \item{stimulus, stimulus_phy}{Stimulus label and physical stimulus identifier.}
#'   \item{cs_trials}{CS+ trial indicator, 0 or 1.}
#'   \item{cs_phy, phy_size}{Physical CS+ value and presented stimulus value.}
#' }
#' @details
#' The stored data contain missing ratings and distances. Use
#' [load_paper_data()] to inspect the original data and
#' [replicate_paper_analysis()] for the supported reanalysis workflow.
#' @source Yu, K., Tuerlinckx, F., Vanpaemel, W., & Zaman, J. (2023).
#'   Humans display interindividual differences in the latent mechanisms
#'   underlying fear generalization behaviour. *Communications Psychology*, 1, 5.
#'   \doi{10.1038/s44271-023-00005-0}
#' @examples
#' data(exp1_data)
#' head(exp1_data)
"exp1_data"

#' Experiment 2 Data from Yu et al. (2023)
#'
#' Trial-level data from the differential conditioning experiment in
#' Yu et al. (2023), with CS+ and CS- cues and distances from each cue.
#'
#' @format A data frame with 7,200 rows and 24 variables:
#' \describe{
#'   \item{participant}{Participant identifier, 1-40.}
#'   \item{trials}{Trial order within participant.}
#'   \item{d_phy_m, d_phy_p}{Physical distances from CS- and CS+, respectively.}
#'   \item{d_per_m, d_per_p}{Perceptual distances from CS- and CS+, respectively.}
#'   \item{block_number, trial_number}{Experimental block and trial within block.}
#'   \item{US}{Observed aversive outcome, 0 or 1.}
#'   \item{startle_circle, startle_iti}{Startle probe indicators during stimulus
#'     presentation and the inter-trial interval.}
#'   \item{per_size}{Perceptual stimulus measurement.}
#'   \item{us_expect}{US expectancy rating on the original 0-10 scale.}
#'   \item{total_trial_nr}{Trial number across the experiment.}
#'   \item{group}{Experimental group assignment.}
#'   \item{stimulus, stimulus_phy}{Stimulus label and physical stimulus identifier.}
#'   \item{csp_trials, csm_trials}{CS+ and CS- trial indicators, respectively.}
#'   \item{us_p}{Excitatory learning target: 0 or 1.}
#'   \item{us_m}{Inhibitory learning target: -1 on CS- trials and 0 otherwise.}
#'   \item{cs_phy_p, cs_phy_m}{Physical CS+ and CS- reference values.}
#'   \item{phy_size}{Physical value of the presented stimulus.}
#' }
#' @details
#' `us_m = -1` is an inhibitory target, not a missing-value code.
#' The stored data retain a legacy `design` attribute reading
#' `"Two CS+ design (2CS)"`; the cue indicators and outcomes identify the
#' CS+/CS- design described here. The raw dataset is preserved unchanged.
#'
#' Use [replicate_paper_analysis()] for the supported reanalysis workflow.
#' @source Yu, K., Tuerlinckx, F., Vanpaemel, W., & Zaman, J. (2023).
#'   Humans display interindividual differences in the latent mechanisms
#'   underlying fear generalization behaviour. *Communications Psychology*, 1, 5.
#'   \doi{10.1038/s44271-023-00005-0}
#' @examples
#' data(exp2_data)
#' with(exp2_data, table(us_m, csm_trials))
"exp2_data"
