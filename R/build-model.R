# HOUSE MODEL DEVELOPMENT
# CORY McCARTAN

library(tidyverse)
library(rstan)
library(rstanarm)
library(brms)
library(loo)
library(tidybayes)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# DATA ##############################c
polls = read_csv("data/current_polls.csv")
results = read_csv("data/combined.csv") %>%
    mutate(lag_pop = lag(logit_pop))
econ = read_csv("data/economy.csv") %>%
    drop_na %>% tail(1)
pres_appr = read_csv("../house-data/trump_polls.csv") %>%
    transmute(approval = approve/100,
              date = mdy(enddate)) %>%
    tail(20) %>%
    summarize(appr = mean(qlogis(approval)))

results_20 = tibble(year=2020, before=235, pres=-1, house=1, midterm=0,
                    retire_lgt=log(9/28), unemp=econ$unemp, gdp=econ$gdp,
                    infl=econ$infl, earn=econ$earn, appr=pres_appr$appr,
                    lag_pop=tail(results, 1)$logit_pop)

####################
#   SEATS WON      #
####################

# TESTING ########
form = seats - before ~ logit_pop + midterm*(pres + I(before-218))
form2 = seats - before ~ logit_pop + midterm*pres + pres:appr +
    midterm*I(before-218) + midterm*retire_lgt

loocv = function(i, form) {
    m = stan_glm(form, data=results[-i,], prior=cauchy(), chains=1,
                  warmup=500, refresh=0)
    test_data = select(results[i,], -logit_pop) %>%
        crossing(logit_pop=rnorm(200, .$logit_mean, .$logit_sd))
    test_pred = as.numeric(posterior_predict(m, newdata=test_data, draws=50))
    test_act = with(results[i,], seats - before)
    mean(abs(test_pred - test_act))
}

cv_err = map_dbl(c(1, 4:15), ~ loocv(., form))
cv_err2 = map_dbl(c(1, 4:15), ~ loocv(., form2))
cv_err3 = map_dbl(c(1, 4:15), ~ loocv(., form3))

formL = seats ~ before + logit_pop + midterm*pres
m = stan_glm(formL, data=results[-i,], prior=cauchy(), chains=1,
             warmup=500, refresh=0)

cvL = map_dbl(c(1, 4:15), ~ loocv(., formL))

qplot(fitted(m), seats, label=year, geom='text', data=results[-i,]) +
    geom_abline(slope=1)


i=15
test_data = select(results[i,], -logit_pop) %>%
    crossing(logit_pop=rnorm(200, .$logit_mean, .$logit_sd))
test_pred = as.numeric(posterior_predict(m, newdata=test_data, draws=50))
qplot(test_pred) + geom_vline(xintercept=results$seats[i])

plot(fitted(m), resid(m))
plot(resid(m))
qqnorm(resid(m))
quantile(test_pred-194, c(0.025, 0.05, 0.5, 0.95, 0.975))


# FINAL MODEL

results_model = stan_glm(seats ~ before + logit_pop + midterm*pres + pres*house,
                         data=results, prior=cauchy(), chains=2,
                         warmup=500, iter=1500)

####################
#   VOTER INTENT   #
####################
prior_model = stan_glm(logit_pop ~ pres:(appr+earn+gdp+unemp) + midterm*pres + lag_pop,
                       data=results, prior=normal(), chains=1, warmup=500)
#post_final = posterior_predict(prior_model, newdata=results[23,])[,1] # 2018
post_final = posterior_predict(prior_model, newdata=results_20)[,1]
c(mean(post_final),  sd(post_final))

model_d = select(polls, -date) %>%
    compose_data(.n_name=n_prefix("N"), W = .$n_weeks[1],
                 prior_final_mean = mean(post_final),
                 prior_final_sd = 1.5*sd(post_final)/4, # arb. inflation for model misspec.
                 prior_natl_poll_error = 4*0.011/2, # pew number adjusted for other errors
                 prior_rv_bias = 0.011, # RV and A bias from 538
                 prior_lv_bias = 0.0,
                 prior_a_bias = 0.02,
                 lv_rv_ratio = 5)

intent_model = stan_model("stan/intent-model.stan")
intent_fit = sampling(intent_model, data=model_d, chains=2, iter=1500, warmup=500,
                      pars=c("sigma", "sigma_u", "sd_firm", "bias_rv", "bias_lv",
                             "house_effects", "natl_dem", "log_lik"),
                      control=list(adapt_delta=0.95, max_treedepth=15))
print(intent_fit, pars=c("natl_dem", "log_lik", "house_effects"), include=F)
plot(intent_fit, pars=c("natl_dem", "log_lik", "house_effects", "lp__"), include=F)
pairs(intent_fit, pars=c("bias_rv", "bias_lv"), include=T)

draws = recover_types(intent_fit) %>% spread_draws(natl_dem[week])
ggplot(draws, aes(week, natl_dem)) +
    geom_point(aes(week, dem, shape=type_lv, color=type_lv), data=polls) +
    stat_lineribbon(alpha=0.3, fill="#666666") +
    theme_minimal()

firms = recover_types(intent_fit, polls) %>%
    spread_draws(house_effects[firm])
firms %>%
    group_by(firm) %>%
    summarize(house_effect = median(house_effects)) %>%
    arrange(desc(house_effect)) %>%
    View

# FINAL INFERENCES

logit_pop = draws %>%
    filter(week==n_weeks) %>%
    pull(natl_dem) %>%
    qlogis %>%
    sample(100)

logit_pop = rnorm(100, 0.1183, 0.1)

inf_data = crossing(results_20, logit_pop=logit_pop)
inf_pred = as.numeric(posterior_predict(results_model, newdata=inf_data, draws=500))
qplot(inf_pred) + geom_vline(xintercept=218)
mean(inf_pred >= 218)

test_data = crossing(results_20, logit_pop=qlogis(seq(0.4, 0.6, 0.01)))
hyp_pred = posterior_predict(results_model, newdata=test_data, draws=500)
qplot(seq(0.4, 0.6, 0.01), colMeans(hyp_pred>=218))
