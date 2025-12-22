# globals.R
# Global variable declarations for R CMD check
# These objects are created dynamically (NSE, model internals, data columns)

utils::globalVariables(c(

  ## tidy evaluation / NSE
  ".data",

  ## Common data and plotting variables
  "time", "time_var", "id", "y", "predY",
  "CI.inf", "CI.sup", "CV", "xilim",
  "st.1", "wk.1",

  ## Cumulative incidence / risk estimates
  "Cum_01_2.5", "Cum_01_97.5", "Cum_01_est",
  "Cum_02_2.5", "Cum_02_97.5", "Cum_02_est",
  "Cum_12_2.5", "Cum_12_97.5", "Cum_12_est",

  ## Baseline hazard / parametric survival components
  "Gompertz", "Weibull", "G_W",
  "Gompert.1_01", "Gompert.1_01.se", "Gompert.1_01.name",
  "Gompert.2_01", "Gompert.2_01.se", "Gompert.2_01.name",
  "Gompert.1_02", "Gompert.1_02.se", "Gompert.1_02.name",
  "Gompert.2_02", "Gompert.2_02.se", "Gompert.2_02.name",
  "Gompert.1_12", "Gompert.1_12.se", "Gompert.1_12.name",
  "Gompert.2_12", "Gompert.2_12.se", "Gompert.2_12.name",

  ## Longitudinal / association parameters
  "alpha_b", "alpha_b_01", "alpha_b_02",
  "alpha_y_slope_var",
  "param_mean", "random.effects",

  ## Internal model objects and data structures
  "Objectlsjm", "object",
  "data.long", "data.long.Case1",
  "data.id", "data.id.1",
  "U_base", "sigma_long", "corr_intra_inter",

  ## Loop indices and case-specific helpers
  "id.boucle",
  "id.Case1_boucle", "id.Case1bis_boucle",
  "id.Case2_boucle", "id.Case3_boucle",
  "index_b_slope", "nbCase1", "left_trunc",

  ## Model formulas and grouping structures
  "formGroup", "formFixedVar", "formRandomVar",
  "hazard_baseline_12"
))
