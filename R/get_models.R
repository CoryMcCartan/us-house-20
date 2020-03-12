library(rstan)
library(readr)
library(rstanarm)
library(tidybayes)

get_prior_m = function(path, results, recompile=F) {
    compiled_model = paste0(path, ".rdata")
    if (recompile && file.exists(compiled_model))
        file.remove(compiled_model)

    if (file.exists(compiled_model)) {
        prior_model = read_rds(compiled_model)
    } else {
        prior_model = stan_glm(logit_pop ~ pres:(appr+earn+gdp+unemp) +
                                   midterm*pres + lag_pop,
                               data=results, prior=normal(), chains=1, warmup=500)
        write_rds(prior_model, compiled_model, compress="gz")
    }

    prior_model
}

get_intent_m = function(path, recompile=F) {
    compiled_model = paste0(path, ".rdata")
    if (recompile && file.exists(compiled_model))
        file.remove(compiled_model)

    if (file.exists(compiled_model)) {
        intent_model = read_rds(compiled_model)
    } else {
        intent_model =  stan_model(paste0(path, ".stan"))
        write_rds(intent_model, compiled_model, compress="gz")
    }

    intent_model
}

get_results_m = function(path, results, recompile=F) {
    compiled_model = paste0(path, ".rdata")
    if (recompile && file.exists(compiled_model))
        file.remove(compiled_model)

    if (file.exists(compiled_model)) {
        results_model = read_rds(compiled_model)
    } else {
        results_model = stan_glm(seats ~ before + logit_pop + midterm*pres +
                                     pres*house,
                                 data=results, prior=cauchy(), chains=2,
                                 warmup=500, iter=1500)
        write_rds(results_model, compiled_model, compress="gz")
    }

    results_model
}
