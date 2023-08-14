source('code/icc_functions.R')
source('code/utils.R')

#' Runs a single brms model defined by the arguments.
#'
#' @param meta list with non-compound variables that are needed to fit the model, e.g. subject, center, age etc.
#' @param compound_vector vector that contains the data of a compound
#' @param compound_string string that describes name of compound in compound_vector
#' @param compiled_model fitted brms model. If options$pre_compiled priors from this compiled will be used to fit the model
#' @param model_spec see run_experiment function
#' @param params see run_experiment function
#' @param options see run_experiment function
#'
#' @return summary statistics of the model and calculated statistics (i.e. ICC and its components)
#' write summary statistics to {params$RESULTS_PATH}/temp/{compound_string}.RDS
#'
#' @examples see run_experiment function
fit_brms_model <- function(meta, compound_vector, compound_string, compiled_model, model_spec, params, options) {
  data <- meta |> dplyr::bind_cols("y" = compound_vector)
  processed_compound <- process_compound(data)
  
  # print(glue::glue("{compound_string}: {processed_compound$lod}"))
  
  if(options$pre_compiled) {
    fit <- brms::update(compiled_model,
                  newdata = processed_compound$data,
                  backend = params$BRMS_BACKEND,
                  cores = params$CORES_BRMS,
                  iter = params$NR_MODEL_ITERATIONS,
                  seed = params$SEED,
                  control = list(adapt_delta = params$ADAPT_DELTA))
    
    nr_of_levels <- params$nr_of_levels 
  } else {
    fit <- brms::brm(data = processed_compound$data, 
               family = model_spec$FAMILY, 
               formula = model_spec$FORMULA,
               backend = params$BRMS_BACKEND,
               cores = params$CORES_BRMS,
               iter = params$NR_MODEL_ITERATIONS,
               seed = params$SEED,
               control = list(adapt_delta = params$ADAPT_DELTA))
    
    nr_of_levels <- length(brms::ranef(fit)) + 1
  }
  
  p_samples <- brms::as_draws_df(fit)
  
  icc_draws <- icc_wrapper(fit, p_samples, nr_of_levels)
  
  summary_stats <- 
    posterior::summarise_draws(icc_draws, "mean","median","sd","mad",
                               ~quantile(.x, probs = c(0.025, 0.975)),
                               "rhat", "ess_bulk", "ess_tail") |> 
    dplyr::bind_cols("compound" = compound_string, 
              "lod" = processed_compound$lod)
  
  summary_stats |> saveRDS(here::here(glue::glue("{params$RESULTS_PATH}/temp/{compound_string}.RDS"))) 
  
  return(summary_stats)
}

#' Function that facilitates in calling fit_brms_model many times. 
#'
#' @param df dataframe with covariates and lcms data
#' @param annotations dataframe with the annotations
#' @param model_spec see run_experiment function
#' @param options see run_experiment function
#' @param params see run_experiment function
#'
#' @return results of experiment
#' @export
#'
#' @examples see run_experiment
fit_brms_model_wrapper <- function(df, annotations, model_spec, options, params) {
  compounds <- df |> dplyr::select(tidyselect::starts_with('X'))
  compounds <- compounds[,sample(ncol(compounds))] # shuffle, maybe there's a bias in ordering and model fitting becomes harder later on etc.
  compound_freq <- compounds |> purrr::map_df(~ perc_detected(.x)) |> tidyr::gather() |> dplyr::rename("compound" = "key", "perc_present" = value)
  
  meta_vars <- insight::find_predictors(model_spec$FORMULA)$conditional
  meta <- df |> dplyr::select(tidyselect::all_of(meta_vars))
  
  annotated_masses_strings <- annotations |> dplyr::pull(compound)
  annotated_compounds <- compounds |> dplyr::select(tidyselect::any_of(annotated_masses_strings))
  
  if (options$ANNOTATED_ONLY) {
    compounds <- annotated_compounds
  }
  
  if (options$pre_compiled) {
    #### initial fit to avoid recompilation (= 15 times faster) ####
    random_compound <- compounds[sample(1:ncol(compounds), 1)] |> dplyr::pull()
    df_random_compound <- meta |> dplyr::bind_cols("y" = random_compound)
    df_random_compound_processed <- process_compound(df_random_compound)
    
    compiled_model <- brms::brm(data = df_random_compound_processed$data, 
                          family = model_spec$FAMILY, 
                          formula = model_spec$FORMULA,
                          backend = params$BRMS_BACKEND,
                          cores = params$CORES_BRMS,
                          iter = params$NR_MODEL_ITERATIONS,
                          seed = params$SEED,
                          control = list(adapt_delta = params$ADAPT_DELTA))
    
    param$NR_OF_LEVELS <- length(brms::ranef(compiled_model)) + 1
  }
 
  # create temp directory so you can gauge progress
  temp_directory <- glue::glue("{params$RESULTS_PATH}/temp")
  unlink(here::here(temp_directory), recursive = TRUE) # remove if it exists from previous experiment
  dir.create(here::here(temp_directory))
  
  results <- compounds |>
    furrr::future_imap(~ fit_brms_model(meta, .x, .y, compiled_model, model_spec, params, options), 
                .options = furrr::furrr_options(seed = params$SEED)) |>
    dplyr::bind_rows() |> 
    dplyr::left_join(compound_freq)
  
  return(results)
}