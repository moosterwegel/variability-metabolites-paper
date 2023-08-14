source('code/utils.R')
source('code/brms_functions.R')

#' Wrapper function that facilitates running the models reported in the paper by predefining them.
#'
#' @param experiment string of either 'three-level-model', 'two-level-model', 'three-level-model-fixef', 'three-level-model-fixef-traf', 'two-level-model-stratified-by-center'
#' @param nr_threads_arg number of CPU threads to pass on to future::plan()
#' @param strat_arg a string that describes by what variable the analysis needs to be stratified (NULL if none). With stratification we mean fitting separate models for the different levels of the STRATIFICATION_VAR 
#' @param annotated_only_arg TRUE if only the annotated compounds need to be analyzed
#' @param pre_comp_arg do the brms models need to be precompiled?
#'
#' @return
#' @export
#'
#' @examples see 1_runs_paper.R for examples
run_experiment_wrapper <- function(experiment, nr_threads_arg, strat_arg = NULL, annotated_only_arg = FALSE, pre_comp_arg = FALSE) {
  OPTIONS <- list(
    STRATIFICATION_VAR = strat_arg,
    ANNOTATED_ONLY = annotated_only_arg,
    NR_THREADS = nr_threads_arg,
    pre_compiled = pre_comp_arg
  )
  
  if(experiment == 'three-level-model') {
    PARAMS <- list(
      SEED = 1305,
      CORES_BRMS = 1,
      NR_MODEL_ITERATIONS = 10000,
      ADAPT_DELTA = 0.99,
      BRMS_BACKEND = 'cmdstanr',
      RESULTS_PATH = 'results/unadjusted/',
      RESULTS_FILE_NAME = 'three_level_model.csv'
    )
    
    model <- brms::brmsformula(log(y1) | cens(cen1) ~ 1 + (1 | centre / subjectid))
    
  } else if(experiment == 'two-level-model') {
    PARAMS <- list(
      SEED = 1305,
      CORES_BRMS = 1,
      NR_MODEL_ITERATIONS = 4000,
      ADAPT_DELTA = 0.95,
      BRMS_BACKEND = 'cmdstanr',
      RESULTS_PATH = 'results/unadjusted/',
      RESULTS_FILE_NAME = 'two_level_model.csv'
    )
    
    model <- brms::brmsformula(log(y1) | cens(cen1) ~ 1 + (1 | subjectid))
    
  } else if(experiment == 'three-level-model-fixef') {
    
    PARAMS <- list(
      SEED = 1305,
      CORES_BRMS = 4,
      NR_MODEL_ITERATIONS = 10000,
      ADAPT_DELTA = 0.99,
      BRMS_BACKEND = 'cmdstanr',
      RESULTS_PATH = 'results/adjusted/',
      RESULTS_FILE_NAME = 'three_level_model_all_w_fixef_smooth.csv'
    )
    OPTIONS$ANNOTATED_ONLY <- TRUE
    
    model <- brms::brmsformula(log(y1) | cens(cen1) ~ 1 + s(age_pem) + sq_sex + bmi + (1 | centre / subjectid))
  
    } else if(experiment == 'three-level-model-fixef-traf') {
      
      PARAMS <- list(
        SEED = 1305,
        CORES_BRMS = 4,
        NR_MODEL_ITERATIONS = 10000,
        ADAPT_DELTA = 0.99,
        BRMS_BACKEND = 'cmdstanr',
        RESULTS_PATH = 'results/adjusted/',
        RESULTS_FILE_NAME = 'three_level_model_all_w_fixef_smooth_traf.csv'
      )
      OPTIONS$ANNOTATED_ONLY <- TRUE
      
      model <- brms::brmsformula(log(y1) | cens(cen1) ~ 1 + s(age_pem) + sq_sex + bmi + traf + (1 | centre / subjectid))
      
    } else if(experiment == 'two-level-model-stratified-by-center') {
    
    PARAMS <- list(
      SEED = 1305,
      CORES_BRMS = 1,
      NR_MODEL_ITERATIONS = 4000,
      ADAPT_DELTA = 0.95,
      BRMS_BACKEND = 'cmdstanr',
      RESULTS_PATH = 'results/unadjusted/',
      RESULTS_FILE_NAME = 'two_level_model_strat_by_centre.csv'
    )
    
    model <- brms::brmsformula(log(y1) | cens(cen1) ~ 1 + (1 | subjectid)) # two-level model
    OPTIONS$STRATIFICATION_VAR <- 'centre'
    
  } else {
     print("This experiment was not recognized")
     stop()
  }
  
  MODEL_SPEC <- list(
    FAMILY = stats::gaussian(),
    FORMULA = model
  )

  run_experiment(MODEL_SPEC, OPTIONS, PARAMS)
}

#' Runs the models with the settings defined in its arguments. `run_experiment_wrapper` contains pre-defined experiments.
#' 
#' @param model_spec list specifying family (FAMILY) and formula (FORMULA) of model.
#' @param options list specifying the following:
#' pre_compiled: do the brms models need to be precompiled?,
#' NR_THREADS: number of CPU threads to pass on to future::plan(), 
#' ANNOTATED_ONLY: TRUE if only the annotated compounds need to be analyzed 
#' STRATIFICATION_VAR: a string that describes by what variable the analysis needs to be stratified (NULL if none). With stratification we mean fitting separate models for the different levels of the STRATIFICATION_VAR 
#' @param params list specifying the settings for the brms model.
#'
#' @return
#' writes results to params$RESULTS_PATH}{params$RESULTS_FILE_NAME}
#'
#' @examples see `run_experiment_wrapper`
run_experiment <- function(model_spec, options, params) {
  future::plan(future::multisession, workers = options$NR_THREADS)
  set.seed(params$SEED)
  
  data <- get_processed_data()
  df <- data$covariates |> dplyr::left_join(data$lcms) |> dplyr::mutate(subjectid = as.factor(subjectid))
  
  annotations <- get_annotations()
  
  ########### calling functions with set parameters ########### 
  if (!is.null(options$STRATIFICATION_VAR)) {
    results <- df |>
      dplyr::filter(!is.na(.data[[options$STRATIFICATION_VAR]])) |> # prevent NA groups from forming
      dplyr::group_by(.data[[options$STRATIFICATION_VAR]]) |>
      dplyr::group_modify(~ fit_brms_model_wrapper(.x, annotations, model_spec, options, params))
  } else {
    results <- df |>
      dplyr::group_by(NULL) |> 
      dplyr::group_modify(~ fit_brms_model_wrapper(.x, annotations, model_spec, options, params))
  }
  
  ########### write results ########### 
  results |> readr::write_csv(here::here(glue::glue('{params$RESULTS_PATH}{params$RESULTS_FILE_NAME}')))
}
