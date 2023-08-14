source('code/utils.R')

annotations <- get_annotations()
annotations_incl_pathways <- get_annotations_plus_kegg()

results <- readr::read_csv(here::here('results/unadjusted/three_level_model.csv'), show_col_types = FALSE) |> 
  preprocess_results(annotations)

filtered_features <- results |> dplyr::filter(variable == 'icc_log_level_2') |> 
  dplyr::select(compound, annotation, identified) |> dplyr::distinct()
filtered_features_identified <- filtered_features |> dplyr::filter(identified)

results_adjusted_smooth <- readr::read_csv(here::here('results/adjusted/three_level_model_all_w_fixef_smooth.csv'), show_col_types = FALSE) |> 
  preprocess_base_layer(annotations) |> 
  dplyr::filter(compound %in% filtered_features_identified$compound)  # use processed features from the unadjusted model (= fair comparison)

results_adjusted_smooth_traf <- readr::read_csv(here::here('results/adjusted/three_level_model_all_w_fixef_smooth_traf.csv'), show_col_types = FALSE) |> 
  preprocess_base_layer(annotations) |> 
  dplyr::filter(compound %in% filtered_features_identified$compound)  # use processed features from the unadjusted model

results_two_level <- readr::read_csv(here::here('results/unadjusted/two_level_model.csv'), show_col_types = FALSE) |> 
  preprocess_base_layer(annotations) |> 
  dplyr::filter(compound %in% filtered_features$compound) # use processed features from the unadjusted three-level model

results_two_level_strat_by_center <- readr::read_csv(here::here('results/unadjusted/two_level_model_strat_by_centre.csv'), show_col_types = FALSE) |> 
  preprocess_base_layer(annotations) |> 
  dplyr::filter(variable == 'icc_log', compound %in% filtered_features_identified$compound) # to make an eventual comparison with three-level model, use same features

#### create numerical summaries ####
get_stats <- \(mean, median) tibble::tibble(median_mean = stats::median(mean), 
                                    mean_mean = mean(mean), 
                                    median_median = stats::median(median),
                                    mean_median = mean(median),
                                    iqr_median = stats::IQR(median),
                                    q25_median = stats::quantile(median, 0.25),
                                    q75_median = stats::quantile(median, 0.75),
                                    n = dplyr::n())

summary <- results |> 
  dplyr::filter(stringr::str_detect(variable, 'icc')) |> 
  dplyr::group_by(variable, identified) |> 
  dplyr::summarize(get_stats(mean, median))
summary

summary_classes <- results_w_classes(results) |> 
  dplyr::filter(stringr::str_detect(variable, 'icc'), identified, !is.na(class)) |> 
  dplyr::group_by(variable, class) |> 
  dplyr::summarize(get_stats(mean, median))
summary_classes

summary_pathways <- results_w_pathways(results) |> 
  dplyr::filter(stringr::str_detect(variable, 'icc'), identified, !is.na(pathway_name)) |> 
  dplyr::group_by(variable, pathway_name) |> 
  dplyr::summarize(get_stats(mean, median))
summary_pathways

# classes highest
summary_classes |> 
  dplyr::filter(variable == 'icc_log_level_2', !is.na(class), n >= 7) |> dplyr::arrange(dplyr::desc(median_median)) |> print(n = 25)

summary_pathways |> 
  dplyr::filter(variable == 'icc_log_level_2', !is.na(pathway_name), n >= 4) |> dplyr::arrange(dplyr::desc(median_median)) |> print(n = 25)

# compounds with kegg entry
results_w_pathways(results) |> 
  dplyr::filter(identified, variable == 'icc_log_level_2', !is.na(kegg_entry)) |> 
  dplyr::select(feature_belongs_to_same_compound) |> dplyr::distinct() |> dplyr::tally()

# compounds with pathway
results_w_pathways(results) |> 
  dplyr::filter(identified, variable == 'icc_log_level_2', !is.na(pathway_name)) |> 
  dplyr::select(feature_belongs_to_same_compound) |> dplyr::distinct() |> dplyr::tally()

# unique compounds
summary |> dplyr::select(identified, n)

# top 10 compounds
results |> dplyr::filter(identified, variable == 'icc_log_level_2') |> 
  dplyr::arrange(dplyr::desc(median)) |> 
  dplyr::select(annotation, curated_name, median, '2.5%', '97.5%', perc_present) |> 
  print(n = 10)

# nr of annotations
annotations |> dim() |> purrr::pluck(1)

# nr of times annotation referred to same identity
annotations |> 
  dplyr::group_by(feature_belongs_to_same_compound) |> 
  dplyr::tally() |> dplyr::group_by(n) |> dplyr::tally() |> 
  dplyr::rename("number of times identity was annotated" = n, "number of times this occured" = nn)

# nr classes
# number of compounds per class
results_w_classes(results) |> 
  dplyr::filter(identified, variable == 'icc_log_level_2', !is.na(class)) |> 
  dplyr::group_by(class) |> dplyr::tally() |> dplyr::arrange(dplyr::desc(n)) |> print(n = 25)

# number of exogenous compounds
annotations |> dplyr::select(feature_belongs_to_same_compound, exogenous) |> 
  dplyr::distinct() |> dplyr::group_by(exogenous) |> dplyr::tally() |> dplyr::arrange(dplyr::desc(n)) 

# exogenous average
results |> dplyr::filter(variable == "icc_log_level_2", identified) |> 
  dplyr::group_by(exogenous) |> 
  dplyr::summarize(get_stats(mean, median)) 

## exogenous, per class
results_w_classes(results) |> 
  dplyr::filter(variable == 'icc_log_level_2', identified, exogenous == 1) |> 
  dplyr::select(curated_name, class) 

summary_adjusted_smooth <- results_adjusted_smooth |> 
  dplyr::filter(stringr::str_detect(variable, 'icc')) |> 
  dplyr::group_by(variable) |> 
  dplyr::summarize(get_stats(mean, median))
summary_adjusted_smooth

summary_adjusted_smooth_traf <- results_adjusted_smooth_traf |> 
  dplyr::filter(stringr::str_detect(variable, 'icc')) |> 
  dplyr::group_by(variable) |> 
  dplyr::summarize(get_stats(mean, median))
summary_adjusted_smooth_traf

summary_by_centre <- results_two_level_strat_by_center |> 
  dplyr::filter(stringr::str_detect(variable, 'icc_log')) |> 
  dplyr::group_by(variable, centre) |> 
  dplyr::summarize(get_stats(mean, median)) 

summary_two_level <- results_two_level |>
  dplyr::filter(stringr::str_detect(variable, 'icc')) |> 
  dplyr::group_by(variable, identified) |> 
  dplyr::summarize(get_stats(mean, median))

#### graphical summaries ####
labels <- list(x_icc = "ICC of compound (log scale)") 

color_palette <- c(unname(grDevices::palette.colors()['skyblue']), unname(grDevices::palette.colors()['bluishgreen']))
color_palette_soft <- c('#b2df8a', '#a6cee3')
okabe_2 <- c(unname(grDevices::palette.colors()['gray']), unname(grDevices::palette.colors()['blue']))
okabe_4 <- c(okabe_2, unname(grDevices::palette.colors()['bluishgreen']), unname(grDevices::palette.colors()['orange']))

ggplot2::theme_set(
  ggplot2::theme_light() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold", size = 9),
      strip.text = ggplot2::element_text(face = "bold", 
                                size = 8, 
                                hjust = 0),
      legend.position = 'none')
)

font_sizes = ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8), 
                   axis.text.x = ggplot2::element_text(size = 7), 
                   axis.title.x = ggplot2::element_text(size = 9), 
                   axis.title.y = ggplot2::element_text(size = 9))

p_single_icc <- function(results, summary, icc_string, color_palette, by_annotation) {
  p <- results |> 
    dplyr::filter(variable == icc_string) |> 
    ggplot2::ggplot(ggplot2::aes(x = median, y = ggplot2::after_stat(density)))
  
  if (by_annotation) {
    
    nr_identified <- summary[summary$variable == icc_string, 'n']
    
    p <- p + 
      ggplot2::geom_histogram(bins = 20, position = "identity", alpha = 0.8, color = 'white', ggplot2::aes(fill = identified)) +
      ggplot2::facet_wrap(~factor(identified, levels = c('TRUE', 'FALSE')), nrow = 2, 
                         labeller = ggplot2::as_labeller(c('TRUE' = glue::glue('Identified compounds (n = {nr_identified$n[2]})'), 
                                                  'FALSE' = glue::glue('Unidentified features (n = {nr_identified$n[1]})'))))
  } else {
    p <- p + 
      ggplot2::geom_histogram(bins = 20, position = "identity", alpha = 0.8, color = 'white', ggplot2::aes(fill = variable))
  }
  p + 
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1),
                       minor_breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::geom_vline(data = summary |> dplyr::filter(variable == icc_string), 
               mapping = ggplot2::aes(xintercept = median_median), linetype = 'longdash', 
               linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = color_palette) +
    ggplot2::geom_text(data = summary |> dplyr::filter(variable == icc_string), 
              mapping = ggplot2::aes(x = median_median, label = round(median_median, digits = 2), 
                            y = 0.1, hjust = 1.2), size = 3) +
    ggplot2::labs(x = labels$x_icc)
}

p_icc_by_class <- function(results, summary, threshold, icc_string, color_palette) {
  class_n <- summary |> dplyr::filter(variable == icc_string) |> 
    dplyr::select(variable, class, n)

  results |> dplyr::left_join(class_n) |> 
    dplyr::filter(variable == icc_string, identified, !is.na(class), n >= threshold) |>
    ggplot2::ggplot(ggplot2::aes(x = median, y = tolower(forcats::fct_rev(forcats::fct_infreq(class))))) + 
    ggplot2::geom_boxplot(color = 'black', fill = color_palette[2], size = 0.4, outlier.size = 0.4) +
    ggplot2::geom_jitter(alpha = 0.1, size = 0.8) +
    ggplot2::scale_x_continuous(limits = c(0, 1), # introduces warnings on removing
                       breaks = seq(0, 1, 0.1),
                       minor_breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    ggplot2::labs(x = labels$x_icc, y = "Class")
}

p_icc_by_pathway <- function(results, summary, threshold, icc_string, color_palette) {
  class_n <- summary |> dplyr::filter(variable == icc_string) |> 
    dplyr::select(variable, pathway_name, n)
  
  results |> dplyr::left_join(class_n) |> 
    dplyr::filter(variable == icc_string, identified, 
           !is.na(pathway_name), n >= threshold, pathway_name != 'Metabolic pathways') |>
    ggplot2::ggplot(ggplot2::aes(x = median, y = tolower(forcats::fct_rev(forcats::fct_infreq(pathway_name)))))+ 
    ggplot2::geom_boxplot(color = 'black', fill = color_palette[2], size = 0.4, outlier.size = 0.4) +
    ggplot2::geom_jitter(alpha = 0.1, size = 0.8) +
    ggplot2::scale_x_continuous(limits = c(0, 1), # introduces warnings on removing
                       breaks = seq(0, 1, 0.1),
                       minor_breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    ggplot2::labs(x = labels$x_icc, y = "Pathway")
}

p_by_annotation_status <- p_single_icc(results, summary, "icc_log_level_2", color_palette_soft, TRUE)
# ggsave(filename = here('figures/icc_three_level_model.tiff'),
#        plot = p_by_annotation_status, compression = "lzw")

p_classes <- p_icc_by_class(results_w_classes(results), summary_classes, 3, "icc_log_level_2", color_palette_soft)
# ggsave(filename = here('figures/icc_three_level_model_classes.tiff'),
#        plot = p_classes, compression = "lzw")


p_pathways <- p_icc_by_pathway(results_w_pathways(results), summary_pathways, 4,'icc_log_level_2', color_palette_soft)
# ggsave(filename = here('figures/icc_three_level_model_pathways.tiff'),
#        plot = p_pathways, compression = "lzw")

# figure 1
# maximum width 19.05cm (7.5 inch), tiff LZW compression
p_classes_adj <- p_classes + 
  font_sizes +
  ggplot2::theme(axis.title.y = ggplot2::element_blank(), 
         axis.title.x = ggplot2::element_blank(),
         axis.text.x = ggplot2::element_blank(),
         axis.ticks.x = ggplot2::element_blank()) +
  ggplot2::ggtitle('Classes')

p_pathways_adj <- p_pathways + 
  font_sizes +
  ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
  ggplot2::ggtitle('Pathways')

p_by_annotation_status_adj <- p_by_annotation_status + font_sizes
layout <- c(patchwork::area(1, 1, 24, 4), # 5 / 8
            patchwork::area(1, 5, 15, 7), 
            patchwork::area(16, 5, 24, 7))

fig_1 <- (p_by_annotation_status_adj + p_classes_adj + p_pathways_adj) + patchwork::plot_layout(design = layout)
ggplot2::ggsave(here::here('figures/fig_1.tiff'),
       plot = fig_1, compression = 'lzw', width = 19, height = 10, units = 'cm')

#### SUPPLEMENTAL FIGURES ####
p_log_vs_data_scale <- results_two_level |> 
  dplyr::filter(stringr::str_detect(variable, 'icc')) |> 
  ggplot2::ggplot(ggplot2::aes(x = mean, y = ggplot2::after_stat(density))) +
  ggplot2::geom_histogram(ggplot2::aes(fill = identified), bins = 20, position = "identity", alpha = 0.6, color = 'white') +
  ggplot2::facet_wrap(ggplot2::vars(variable), nrow = 2, labeller = ggplot2::as_labeller(c(icc_data  = 'On data scale', icc_log = 'On log scale'))) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1),
                     minor_breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  ggplot2::coord_cartesian(xlim = c(0, 1)) +
  ggplot2::scale_fill_manual(values = okabe_2) +
  ggplot2::geom_vline(data = summary_two_level, mapping = ggplot2::aes(xintercept = median_median, linetype = identified)) +
  ggplot2::geom_text(data = summary_two_level, mapping = ggplot2::aes(x = median_median, label = round(median_median, digits = 2), y = 0, hjust = 1.2)) + 
  ggplot2::labs(title = "Model: log(y1) | cens(cen1) ~ 1 + (1|subjectid)",
       x = labels$x_icc) +
  ggplot2::theme(legend.position = 'right')

ggplot2::ggsave(filename = here::here('figures/icc_two_level_model_log_vs_data_scale.tiff'),
       plot = p_log_vs_data_scale + font_sizes + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8)), 
       compression = "lzw", width = 14.5, height = 11, unit = 'cm')

var_dic <- c(`icc_log_level_2` = "frac(sigma[3] + sigma[2], sigma[3] + sigma[2] + sigma[1])",
             `icc_log_level_3` = "frac(sigma[3], sigma[3] + sigma[2] + sigma[1])", 
             `icc_log_minus_centre` = "frac(sigma[2], sigma[2] + sigma[1])", 
             `icc_log_minus_centre_nominator` = "frac(sigma[2], sigma[3] + sigma[2] + sigma[1])")

p_iccs_grid <- function(results, summary, label_dic, title_string, by_annotation) {
  p <- results |> 
    dplyr::filter(stringr::str_detect(variable, 'icc_log')) |> 
    ggplot2::ggplot(ggplot2::aes(x = median, y = ggplot2::after_stat(density))) + 
    ggplot2::facet_wrap(ggplot2::vars(variable), 
               labeller = ggplot2::labeller(variable = ggplot2::as_labeller(label_dic, ggplot2::label_parsed)),
               ncol = 2) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1),
                       minor_breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(title = title_string,
         x = labels$x_icc) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", 
                                    size = 14, 
                                    hjust = 0))
  
  if(by_annotation) {
    p <- p + ggplot2::geom_histogram(ggplot2::aes(fill = identified), alpha = 0.5, bins = 20, position = "identity", color = 'white') +
      ggplot2::scale_fill_manual(values = okabe_2) +
      ggplot2::geom_vline(data = summary, mapping = ggplot2::aes(xintercept = median_median, linetype = identified)) +
      ggplot2::theme(legend.position = 'top')
  }
  else {
    p <- p + ggplot2::geom_histogram(ggplot2::aes(fill = variable), alpha = 0.7, bins = 20, position = "identity", color = 'white') +
      ggplot2::scale_fill_manual(values = okabe_4) +
      ggplot2::geom_vline(data = summary, mapping = ggplot2::aes(xintercept = median_median))
  }
  
  p + ggplot2::geom_text(data = summary, mapping = ggplot2::aes(x = median_median, label = round(median_median, digits = 2), y = 0, hjust = 1.2)) 
}

p_all_iccs_m1 <- p_iccs_grid(results, summary, var_dic,
            "Model: log(y1) | cens(cen1) ~ 1 + (1|centre/subjectid) + (1|centre)", TRUE)
ggplot2::ggsave(filename = here::here('figures/all_iccs_main_model.tiff'),
       plot = p_all_iccs_m1 + font_sizes, compression = 'lzw', height = 15, width = 19, unit = 'cm')

p_all_iccs_m2 <- p_iccs_grid(results_adjusted_smooth, summary_adjusted_smooth, var_dic,
                             "Model: log(y1) | cens(cen1) ~ 1 + s(age) + bmi + sex + (1|centre/subjectid) + (1|centre)", FALSE)
ggplot2::ggsave(filename = here::here('figures/all_iccs_adjusted_model.tiff'),
       plot = p_all_iccs_m2 + font_sizes, height = 14, width = 19, unit = 'cm',
       compression = "lzw")

p_classes_adjusted <- p_icc_by_class(results_w_classes(results_adjusted_smooth), summary_classes, 3, "icc_log_level_2", color_palette_soft) +
  font_sizes +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25), limit = c(0, 1))  +
  ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
  ggplot2::ggtitle('Classes')

p_pathways_adjusted <- p_icc_by_pathway(results_w_pathways(results_adjusted_smooth), summary_pathways, 4, 'icc_log_level_2', color_palette_soft) +
  font_sizes +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.25), limit = c(0, 1))  +
  ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
  ggplot2::labs(title = 'Pathways', y = "")

class_pathway_adj <- p_classes_adjusted + p_pathways_adjusted
ggplot2::ggsave(filename = here::here('figures/adjusted_model_pathways_classes.tiff'),
       plot = class_pathway_adj + font_sizes, height = 5, width = 19, unit = 'cm',
       compression = "lzw")

missings <- results |>
  dplyr::filter(stringr::str_detect(variable, 'icc_log_level_2')) |> 
  ggplot2::ggplot(ggplot2::aes(x = perc_present, colour = identified)) +
  ggplot2::stat_ecdf(linewidth = 1.2) +
  ggplot2::scale_colour_manual(values = okabe_4) +
  ggplot2::labs(x = 'Proportion of samples where compound was present', y = 'Proportion of compounds') +
  ggplot2::facet_wrap(~identified, nrow = 1, 
             labeller = ggplot2::as_labeller(c('TRUE' = 'Identified compounds', 'FALSE' = 'Unidentified features')))

ggplot2::ggsave(filename = here::here('figures/missings.tiff'),
       plot = missings + font_sizes + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8)), 
       compression = "lzw", width = 14.5, height = 7, unit = 'cm')

p_icc_by_centre <- results_two_level_strat_by_center |> 
  dplyr::filter(stringr::str_detect(variable, 'icc_log')) |> 
  ggplot2::ggplot(ggplot2::aes(x = median, y = ggplot2::after_stat(density))) + 
  ggplot2::geom_histogram(ggplot2::aes(fill = centre), alpha = 0.7, bins = 20, position = "identity", color = 'white') +
  ggplot2::scale_fill_manual(values = okabe_4) +
  ggplot2::geom_vline(data = summary_by_centre, mapping = ggplot2::aes(xintercept = median_median)) + 
  ggplot2::geom_text(data = summary_by_centre, mapping = ggplot2::aes(x = median_median, label = round(median_median, digits = 2), y = 0, hjust = 1.2)) +
  ggplot2::facet_wrap(dplyr::vars(centre)) +
  ggplot2::scale_y_continuous(limits = c(0, 4)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1),
                     minor_breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  ggplot2::coord_cartesian(xlim = c(0, 1)) +
  ggplot2::labs(x = labels$x_icc) +
  ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", 
                                  size = 14, 
                                  hjust = 0))

ggplot2::ggsave(filename = here::here('figures/iccs_by_centre.tiff'),
       plot = p_icc_by_centre + font_sizes, height = 15, width = 19, unit = 'cm',
       compression = "lzw")

# variance attribution plot
# earlier preprocessing only leaves us with iccs
# first we have to get the annotations
icc_var_plot <- 'icc_log_level_2'

one_annotation_per_compound <- results |> 
  dplyr::filter(identified, variable == icc_var_plot) |> 
  dplyr::pull(annotation)

sds <- c('bw_sd_log_level_3', 'bw_sd_log_level_2', 'wi_sd_log')
results_sds <- readr::read_csv(here::here('results/unadjusted/three_level_model.csv'), show_col_types = FALSE)  |> 
  dplyr::left_join(annotations) |> 
  dplyr::select(compound, annotation, variable, median, curated_name) |> 
  dplyr::filter(variable %in% sds, annotation %in% one_annotation_per_compound) |>
  dplyr::mutate(variable = forcats::fct_relevel(variable, sds))

# ranking such that order annotations appear in is high to low icc and we can split plot in two columns
order <- results |> 
  dplyr::filter(identified, variable == icc_var_plot) |> 
  dplyr::arrange(dplyr::desc(median)) |> 
  dplyr::mutate(icc_rank = dplyr::row_number()) |> 
  dplyr::select(compound, icc_rank)

f_bar_chart <- function(data) {
  data |> 
  ggplot2::ggplot(ggplot2::aes(x  = median, y = stats::reorder(curated_name, -icc_rank), fill = variable)) + 
    ggplot2::geom_bar(position = "fill", stat = "identity") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.5)) +
    ggplot2::scale_fill_manual(name = "Estimated variance: ",
                      values = c("#56B4E9", "#E69F00", "#009E73"), 
                      labels = c("Between-center", "Between-subjects", "Within-subject")) +
    ggplot2::labs(x = "Proportion \n of variance", y = "")
}

p_variance_components_1 <- 
  results_sds |> dplyr::left_join(order) |> dplyr::filter(icc_rank <=  ceiling(length(one_annotation_per_compound)/2)) |> f_bar_chart() 

p_variance_components_2 <- 
  results_sds |> dplyr::left_join(order) |> dplyr::filter(icc_rank > ceiling(length(one_annotation_per_compound)/2)) |> f_bar_chart() 

p_variance_components <- p_variance_components_1 + p_variance_components_2 + 
  patchwork::plot_layout(ncol = 2, guides = 'collect') &
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7),
        legend.position = 'top', 
        legend.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_text(size = 9))

ggplot2::ggsave(here::here('figures/variance_components_identified.tiff'),
       compression = 'lzw',
       dpi = 600, 
       plot = p_variance_components, 
       width = 19, height = 24.5, unit = "cm")

### COMPARISON TARGETED STUDIES ####
# generate template for comparison
readr::read_csv(here::here('results/unadjusted/three_level_model.csv'), show_col_types = FALSE) |> 
  preprocess_results(annotations) |> 
  dplyr::filter(identified, variable == 'icc_log_level_2') |> 
  dplyr::select(annotation, curated_name, median, '2.5%', '97.5%') |> 
  dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, 2))) |> 
  readr::write_csv(here::here('results', 'template_comparison_studies.csv'))

# after filling it in manually
comparison <- readxl::read_excel(here::here('results', 'comparison_targeted_studies.xlsx'))  |> 
  dplyr::filter(!is.na(yin_icc))

comparison_reformat <- comparison |> 
  tidyr::pivot_longer(cols = dplyr::starts_with(c('oosterwegel', 'yin', 'floegel', 'breier')), 
                      names_to = c("study", "variable"),
                      names_sep = "_",
                      values_to = "median") |> 
  tidyr::pivot_wider(names_from = 'variable', values_from = median) |> 
  dplyr::mutate(study = dplyr::recode(study, 
                                      oosterwegel = "Oosterwegel (2023)", 
                                      floegel = "Floegel (2011)", 
                                      yin = "Yin (2022)",
                                      breier = "Breier (2014)"),
                study = forcats::fct_relevel(study,  "Floegel (2011)", "Breier (2014)", "Yin (2022)", "Oosterwegel (2023)"))

comparison_reformat_summary <- comparison_reformat |> 
  dplyr::group_by(study) |> 
  dplyr::summarise(median = stats::median(icc, na.rm = TRUE))

p_comparison <- comparison_reformat |> 
  ggplot2::ggplot(ggplot2::aes(y = study, x = icc)) +
  ggplot2::geom_point() +
  ggplot2::geom_linerange(ggplot2::aes(xmin = `2.5%`, xmax = `97.5%`)) +
  ggplot2::facet_wrap(~curated_name, ncol = 2) +
  ggplot2::geom_vline(data = comparison_reformat_summary, mapping = ggplot2::aes(xintercept = median, linetype = study, color = study),  
                      size = 1.3,
                      alpha = 0.5) +
  ggplot2::scale_x_continuous(limits = c(0,1)) +
  ggplot2::scale_linetype_manual(values = c("dashed", "longdash", "dotdash", "solid")) +
  ggplot2::scale_color_manual(values = okabe_4[c(2, 5, 4, 3)]) +
  ggplot2::labs(x = labels$x_icc, y = "", color = "", linetype = "") +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 4, ncol = 1, byrow = TRUE)) +
  ggplot2::theme(legend.position = 'right',
                 strip.text = ggplot2::element_text(face = "bold",
                                                    margin = ggplot2::margin(0.2, -0.1, 0.1, 0.1, "cm"),
                                                    size = 10,
                                                    hjust = 0))

ggplot2::ggsave(here::here('figures/comparison_targeted.tiff'),
                compression = "lzw",
                plot = p_comparison + font_sizes, width = 19, height = 15, units = 'cm')

# distribution missingness
results |>
  dplyr::filter(variable == 'icc_log_level_2') |> 
  dplyr::group_by(identified) |> 
  dplyr::summarize(median_perc_present= median(perc_present), iqr = stats::IQR(perc_present))

# compounds with missings
results |>
  dplyr::filter(variable == 'icc_log_level_2') |> 
  dplyr::group_by(identified) |> dplyr::summarize(compound_with_missings = sum(perc_present != 1) / sum(!is.na(perc_present)))

### brms diagnostic functions ###
insufficient_diagnostics <- function(results, statistic, diag_measure, sufficient_threshold, operator) {
  nr_of_compounds <- results |> dplyr::filter(variable == statistic) |> dplyr::tally()
  
  insufficient <- results |> 
    dplyr::filter(variable == statistic, operator(.data[[diag_measure]], sufficient_threshold)) 
  
  return(list(perc = dplyr::tally(insufficient) / nr_of_compounds * 100,
              tibble = insufficient))
}

100 - insufficient_diagnostics(results, 'icc_log_level_2', 'ess_tail', 400, `<`)$perc
100 - insufficient_diagnostics(results, 'icc_log_level_2', 'ess_bulk', 400, `<`)$perc
100 - insufficient_diagnostics(results, 'icc_log_level_2', 'rhat', 1.05, `>`)$perc

100 - insufficient_diagnostics(results_adjusted_smooth, 'icc_log_level_2', 'ess_tail', 400, `<`)$perc
100 - insufficient_diagnostics(results_adjusted_smooth, 'icc_log_level_2', 'ess_bulk', 400, `<`)$perc
100 - insufficient_diagnostics(results_adjusted_smooth, 'icc_log_level_2', 'rhat', 1.05, `>`)$perc
