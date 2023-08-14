source('code/run_brms_models.R')

# takes 25 minutes on my machine
run_experiment_wrapper(experiment = 'two-level-model-stratified-by-center', nr_threads_arg = 12, strat_arg = 'centre', annotated_only_arg = TRUE, pre_comp_arg = FALSE)

# takes 2.5 hours on my machine
run_experiment_wrapper(experiment = 'two-level-model', nr_threads_arg = 18, strat_arg = NULL, annotated_only_arg = FALSE, pre_comp_arg = FALSE)

# takes 17.5 - 18 hours on my machine
run_experiment_wrapper(experiment = 'three-level-model', nr_threads_arg = 16, strat_arg = NULL, annotated_only_arg = FALSE, pre_comp_arg = FALSE)

# takes 2.25 hours on my machine, can't be reproduced with open data
run_experiment_wrapper(experiment = 'three-level-model-fixef', nr_threads_arg = 6, strat_arg = NULL, annotated_only_arg = TRUE, pre_comp_arg = FALSE)

# takes 2.25 hours on my machine, can't be reproduced with open data
run_experiment_wrapper(experiment = 'three-level-model-fixef-traf', nr_threads_arg = 6, strat_arg = NULL, annotated_only_arg = TRUE, pre_comp_arg = FALSE)
