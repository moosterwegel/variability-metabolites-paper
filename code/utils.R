BASE_PATH_DATA <- 'data'

perc_missing <- \(x) sum(x == 1) / length(x) # calculates proportion of samples where compound was missing
perc_detected <- \(x) 1 - perc_missing(x) # calculates proportion of samples where compound was detected

#' Processes vector of compound such that missing variables (values below limit of detection (lod) are correctly labelled as censored
#'
#' @param data compound data
#'
#' @return list with lowest detected value (LOD) and vector that labels the censored labels
#' @export
#'
#' @examples
process_compound <- function(data) {
  lod <- data |> dplyr::filter(y > 1) |> dplyr::pull(y) |> min()
  
  compound_processed <- data |> 
    dplyr::mutate(cen1 = dplyr::if_else(y < lod, 'left', 'none'),
           detected = dplyr::if_else(y >= lod, TRUE, FALSE), 
           y1 = dplyr::if_else(detected, y, lod))
  
  return(list("lod" = lod,
              "data" = compound_processed))
}


#' Just a wrapper for reading the data with one line of code
#'
#' @return list with covariate and lcms data
#' @export
#'
#' @examples
get_processed_data  <- function() {
  covariates <- utils::read.csv(glue::glue("{BASE_PATH_DATA}/processed_covariate_data.csv")) 
  lcms <- utils::read.csv(glue::glue("{BASE_PATH_DATA}/processed_lcms_data.csv"))
  
  return(list("covariates" = covariates, 
              "lcms" = lcms))
}

#' Adds ChEBI class to compounds based on methodology described in paper. The selected_classes is result of said methodology.
#'
#' @return dataframe with compound and its class
#' @export
#'
#' @examples
get_processed_ancestors <- function() {
  selected_classes <- c('glycerophospholipid',
                        'O-acylcarnitine', 'amino acid', 'phosphatidylcholine', 'steroid')
  
  df <- readxl::read_excel(here::here('data', 'ancestors_annotations.xlsx'))
  
  df_processed <- df |> dplyr::filter(chebi_ancestors %in% selected_classes) |> 
    dplyr::rename(class = chebi_ancestors)
  
  return(df_processed)
}

#' Get the annotations and clean data a bit.
#'
#' @return dataframe with annotated compounds
#' @export
#'
#' @examples
get_annotations <- function() {
  annotations <- readxl::read_excel(here::here(BASE_PATH_DATA, 'annotations.xlsx')) |> 
    tidyr::extract(feature_belongs_to_same_compound, c("base", "detail"),
            regex = "(\\d+)([[:lower:]])?", 
            remove = FALSE, convert = TRUE) |> 
    tidyr::unite(compound, c(mass, retention_time), sep = '_', remove = FALSE) |> 
    dplyr::mutate(compound = paste("X", compound, sep = ""), 
           base = as.character(base)) |> 
    janitor::clean_names()
  
  return(annotations)
}

#' Wrapper to read the annotations and their KEGG pathways
#'
#' @return dataframe with annotations and their KEGG pathways
#' @export
#'
#' @examples
get_annotations_plus_kegg <- function() {
  annotations <- utils::read.csv(here::here(BASE_PATH_DATA, 'annotations_plus_kegg_pathways.csv'))
  return(annotations)
}


#' Adds ChEBI classes to result dataframe
#'
#' @param results dataframe with results (see 2. summarize_results.R)
#'
#' @return result dataframe with ChEBI class added
#' @export
#'
#' @examples
results_w_classes <- function(results) {
  classes <- get_processed_ancestors()  |> dplyr::select(mass, retention_time, class)
  
  return(results |> dplyr::left_join(classes, dplyr::join_by(mass, retention_time), multiple = 'all')) 
}

#' Adds KEGG pathways to results
#'
#' @param results 
#'
#' @return
#' @export
#'
#' @examples
results_w_pathways <- function(results) {
  annotations_incl_pathways <- get_annotations_plus_kegg() |> 
    dplyr::select(mass, retention_time, kegg_entry, pathway_code, pathway_name)
  
  return(results |> dplyr::left_join(annotations_incl_pathways, dplyr::join_by(mass, retention_time), multiple = 'all', relationship = 'many-to-many')) 
}

#' Processing for ICC results, because some annotations refer to the same identity. If that's the case, we only keep 
#' the one with highest ICC
#' @param df_results 
#'
#' @return processed results
#' @export
#'
#' @examples
keep_most_stable_annotation <- function(df_results) {
  df_results |> 
    dplyr::group_by(feature_belongs_to_same_compound) |> 
    dplyr::arrange(dplyr::desc(median)) |> 
    dplyr::slice(1) |> 
    dplyr::ungroup() 
}

#' Some basic cleaning and reformatting of the data
#'
#' @param df_results 
#' @param df_annotations 
#'
#' @return processed results
#' @export
#'
#' @examples
preprocess_base_layer <- function(df_results, df_annotations) {
  df_results |> 
    dplyr::left_join(df_annotations) |> 
    dplyr::mutate(identified = !is.na(annotation), 
           base = dplyr::if_else(identified, base, compound),
           feature_belongs_to_same_compound = dplyr::if_else(identified, feature_belongs_to_same_compound, base),
           exogenous = dplyr::recode(exogenous, "T" = "1", "F" = "0"))
}

#' String the processing functions together and call the keep_most_stable_annotation function in such a way 
#' (i.e. group_by followed by group_modify) that the most stable annotation of a compound depends on what ICC is of interest. 
#'
#' @param df_results 
#' @param df_annotations 
#'
#' @return processed results
#' @export
#'
#' @examples
preprocess_results <- function(df_results, df_annotations) {
  df_results |> 
    preprocess_base_layer(df_annotations) |> 
    dplyr::filter(stringr::str_detect(variable, "icc")) |> 
    dplyr::group_by(variable) |>
    dplyr::group_modify(~keep_most_stable_annotation(.x)) |> 
    dplyr::ungroup()
}
