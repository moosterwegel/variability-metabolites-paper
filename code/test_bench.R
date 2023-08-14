source('code/utils.R')
source('code/icc_functions.R')
set.seed(1305)

model <- brms::brmsformula(log(y1) | cens(cen1) ~ 1 + (1 | centre / subjectid))
model_spec <- list(
  FAMILY = stats::gaussian(),
  FORMULA = model
)

data <- get_processed_data()
df <- data$covariates |> dplyr::left_join(data$lcms)
annotations <- get_annotations()

compounds <- df |> dplyr::select(starts_with('X'))
compounds <- compounds[,sample(ncol(compounds))]
random_compound <- compounds[sample(1:ncol(compounds), 1)] |> dplyr::pull()
# or pick one yourself, below is 2-Furoylglycine
# random_compound <- compounds['X169.0374_2.1930275'] |> dplyr::pull()

meta_vars <- insight::find_predictors(model_spec$FORMULA)$conditional
meta <- df |> dplyr::select(tidyselect::all_of(meta_vars))

# assemble dataframe to pass to brms
df_random_compound <- meta |> dplyr::bind_cols("y" = random_compound)
# label censored observations correctly
df_random_compound_processed <- process_compound(df_random_compound)

params <- list(
  SEED = 1305,
  CORES_BRMS = 4,
  NR_MODEL_ITERATIONS = 2000,
  ADAPT_DELTA = 0.99,
  BRMS_BACKEND = 'cmdstanr')

fit <- brms::brm(data = df_random_compound_processed$data, 
                 family = model_spec$FAMILY, 
                 formula = model_spec$FORMULA,
                 backend = params$BRMS_BACKEND,
                 cores = params$CORES_BRMS,
                 iter = params$NR_MODEL_ITERATIONS,
                 seed = params$SEED,
                 control = list(adapt_delta = params$ADAPT_DELTA))

nr_of_levels <- length(brms::ranef(fit)) + 1

p_samples <- brms::as_draws_df(fit)
icc_draws <- icc_wrapper(fit, p_samples, nr_of_levels)
posterior::summarise_draws(icc_draws, "mean","median","sd","mad",
                           ~quantile(.x, probs = c(0.025, 0.975)),
                           "rhat", "ess_bulk", "ess_tail")
