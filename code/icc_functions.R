#' Function that selects the right ICC calculation depending on the number of levels in the model (only 2 or 3 supported)
#'
#' @param fit brms fitted model
#' @param p_samples posterior samples from fitted model (result of brms::as_draws_df)
#' @param nr_of_levels 
#'
#' @return posterior draws of the ICC statistic
#' @export
#'
#' @examples
icc_wrapper <- function(fit, p_samples, nr_of_levels) {
  if (nr_of_levels == 2) {
    level_2 <- glue::glue('sd_{fit$ranef$group[1]}__Intercept')
    
    icc_draws <- purrr::pmap_df(
      list(
        # follow notation from formulas in paper
        B0 = p_samples$b_Intercept,
        sigma_2 = p_samples[[level_2]],
        sigma_1 = p_samples$sigma
      ),
      .f = get_iccs_brms
    )
  }
  
  if (nr_of_levels == 3) {
    level_2 <- glue::glue('sd_{fit$ranef$group[2]}__Intercept')
    level_3 <- glue::glue('sd_{fit$ranef$group[1]}__Intercept')
    
    icc_draws <- purrr::pmap_df(
      list(
        # follow notation from formulas in paper
        B0 = p_samples$b_Intercept,
        sigma_3 = p_samples[[level_3]],
        sigma_2 = p_samples[[level_2]],
        sigma_1 = p_samples$sigma
      ),
      .f = get_iccs_three_levels_brms
    )
    
  }
  return(icc_draws)
}

#' Calculate ICC implied by parameters from three-level brms model.
#'
#' @param B0 intercept
#' @param sigma_3 sigma third level (e.g. center level)
#' @param sigma_2 sigma second level (e.g. subject level)
#' @param sigma_1 sigma first level (e.g. measurement level)
#'
#' @return ICC and its components, only on log scale.
#' @export
#'
#' @examples
get_iccs_three_levels_brms <- function(B0, sigma_3, sigma_2, sigma_1) {
  icc_log_level_2 <- (sigma_2^2 + sigma_3^2) / (sigma_1^2 + sigma_2^2 + sigma_3^2)
  icc_log_level_3 <- (sigma_3^2) / (sigma_1^2 + sigma_2^2 + sigma_3^2)
  icc_log_minus_centre_nominator <- (sigma_2^2) / (sigma_1^2 + sigma_2^2 + sigma_3^2)
  icc_log_minus_centre <- (sigma_2^2) / (sigma_1^2 + sigma_2^2)
  
  return(
    list(
      "icc_log_level_2" = icc_log_level_2,
      "icc_log_level_3" = icc_log_level_3,
      "icc_log_minus_centre_nominator" = icc_log_minus_centre_nominator,
      "icc_log_minus_centre" = icc_log_minus_centre,
      "wi_sd_log" = sigma_1,
      "bw_sd_log_level_2" = sigma_2,
      "bw_sd_log_level_3" = sigma_3
    ))
}

#' Calculate ICC implied by parameters from two-level brms model.
#'
#' @param B0 intercept
#' @param sigma_2 sigma second level (e.g. subject level)
#' @param sigma_1 sigma first level (e.g. measurement level)
#' @param ... 
#'
#' @return ICC and its components, on both log and data scale
#' @export
#'
#' @examples
get_iccs_brms <- function(B0, sigma_2, sigma_1, ...) {
  params_data_scale <- get_params_data_scale_brms(B0, sigma_2, sigma_1)
  
  icc_data <- params_data_scale$sigma_2^2 / (params_data_scale$sigma_2^2 + params_data_scale$sigma_1^2)
  icc_log <- sigma_2^2/(sigma_2^2 + sigma_1^2)
  
  return(
    list(
      "icc_data" = icc_data,
      "icc_log" = icc_log,
      "wi_sd_data" = params_data_scale$sigma_1,
      "bw_sd_data" = params_data_scale$sigma_2,
      "wi_sd_log" = sigma_1,
      "bw_sd_log" = sigma_2
    ))
}

#' Convert parameters brms model from log scale to data scale. 
#' Based on method described on adapted from https://rpsychologist.com/GLMM-part1-lognormal 
#'
#' @param B0 intercept
#' @param sigma_2 sigma second level
#' @param sigma_1 sigma first level
#' @param ... 
#'
#' @return parameters on data scale
#' @export
#'
#' @examples
get_params_data_scale_brms <- function(B0, sigma_2, sigma_1, ...) {
  mu <-  stats::integrate(
    f = function(x) {
      exp(x + sigma_1^2 / 2) * stats::dnorm(x, B0, sd = sigma_2)
    },
    lower = B0 - 10 * sigma_2,
    upper = B0 + 10 * sigma_2
  )$value
  
  sigma2_bw <-  stats::integrate(
    f = function(x) {
      (exp(x + sigma_1^2 / 2) - mu)^2  * stats::dnorm(x, B0, sd = sigma_2)
    },
    lower = B0 - 10 * sigma_2,
    upper = B0 + 10 * sigma_2
  )$value
  
  sigma2_wi <-  stats::integrate(
    f = function(x) {
      (exp(sigma_1^2) - 1) * exp(2 * x + sigma_1^2)  * stats::dnorm(x, B0, sd = sigma_2)
    },
    lower = B0 - 10 * sigma_2,
    upper = B0 + 10 * sigma_2
  )$value
  
  return(
    list(
      "sigma_1" = sqrt(sigma2_wi),
      "sigma_2" = sqrt(sigma2_bw),
      "mu" = mu
    ))
}