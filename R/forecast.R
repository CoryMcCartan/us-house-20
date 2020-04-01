#!/usr/bin/env Rscript

#####################################
#     U.S. HOUSE MODEL 2020         #
#     CORY McCARTAN                 #
#     (c) 2020                      #
#####################################

library(optparse)

###
### Parse options
###

option_list = list(
    make_option("--dry", action="store_true", default=F,
                help="Dry run, results not saved."),
    make_option("--date", type="character", default=as.character(Sys.Date()),
                help="The date to estimate from."),
    make_option("--iter", type="integer", default=2000,
                help="Number of MCMC iterations for voter intent estimation,
                      including 500 warmup iterations."),
    make_option("--chains", type="integer", default=2,
                help="Number of MCMC chains for voter intent estimation."),
    make_option("--recompile", action="store_true", default=F,
                help="Force recompile of STAN model."),
    make_option("--model_dir", type="character", default="stan",
                help="The directory in which the models are stored"),
    make_option("--output_file", type="character", default="docs/estimate.json",
                help="The file to save estimates to."),
    make_option("--history_file", type="character", default="docs/history.csv",
                help="The file to save model history to.")
)
opt = parse_args(OptionParser(option_list=option_list,
                              description="Forecast the 2020 U.S. House election."))

suppressMessages(library(tibble))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(tidyr))
suppressMessages(library(readr))
suppressMessages(library(fredr))
suppressMessages(library(lubridate))
suppressMessages(library(glue))
suppressMessages(library(jsonlite))


start_date = ymd("2020-01-01")
from_date = as.Date(opt$date)
intent_model_path = file.path(opt$model_dir, "intent-model")
results_model_path = file.path(opt$model_dir, "results-model")
prior_model_path = file.path(opt$model_dir, "prior-model")

current_seats = 235
election_day = as.Date("2020-11-03")


###
### Download new polling data
###
cat("Downloading data.\n")

source("R/polls.R")
appr_url = "https://projects.fivethirtyeight.com/trump-approval-data/approval_polllist.csv"

# generic ballot
old_polls = suppressMessages(read_csv("docs/polls.csv", col_types="Dcd"))
polls = load_polls(start_date, election_day, write=T) %>%
    filter(date <= from_date) %>%
    drop_na
polls %>%
    select(date, firm, dem) %>%
    write_csv("docs/polls.csv")
if (from_date == Sys.Date() &&
        isTRUE(all.equal(old_polls, select(polls, date, firm, dem)))) {
    cat("No new polls.\n")
    system("osascript -e 'display notification \"No new polls.\" with title \"House Model\"'")
    system("osascript -e beep"); system("osascript -e beep")
    Sys.sleep(10)
}
# pres. approval
pres_appr = suppressMessages(read_csv(appr_url)) %>%
    transmute(approval = approve/100,
              date = mdy(enddate)) %>%
    tail(20) %>%
    summarize(appr = mean(qlogis(approval))) %>%
    pull
# economy
fredr_set_key(Sys.getenv("FRED_KEY"))
econ_params = list(
    series = c("UNRATE", "A939RX0Q048SBEA", "CPIAUCSL", "AHETPI"),
    units = c("lin", "pc1", "pc1", "pc1")
)
econ = pmap_dfr(econ_params, ~ fredr(.x, observation_start=start_date-365,
                              frequency="q", units=.y)) %>%
    pivot_wider(names_from=series_id, values_from=value) %>%
    filter(date <= from_date) %>%
    transmute(unemp = as.numeric(UNRATE) / 100,
              gdp = as.numeric(A939RX0Q048SBEA) / 100,
              infl = as.numeric(CPIAUCSL) / 100,
              earn = as.numeric(AHETPI) / 100) %>%
    fill(everything()) %>%
    tail(1)

# everything
results = suppressMessages(read_csv("data/combined.csv")) %>%
    mutate(lag_pop = lag(logit_pop))

results_20 = tibble(year=2020, before=current_seats, pres=-1, house=1, midterm=0,
                    retire_lgt=log(9/28), unemp=econ$unemp, gdp=econ$gdp,
                    infl=econ$infl, earn=econ$earn, appr=pres_appr,
                    lag_pop=tail(results, 1)$logit_pop)

###
### Estimate voter intent
###

cat("Estimating voter intent.\n")

suppressMessages(library(rstan))
suppressMessages(library(rstanarm))
options(mc.cores=4)
source("R/get_models.R")

prior_model = get_prior_m(prior_model_path, results, opt$recompile)
post_final = posterior_predict(prior_model, newdata=results_20)[,1]

intent_d = select(polls, -date) %>%
    compose_data(.n_name=n_prefix("N"), W = .$n_weeks[1] + 1,
                 prior_final_mean = mean(post_final),
                 prior_final_sd = 2*sd(post_final)/4, # arb. inflation for model misspec.
                 prior_natl_poll_error = 4*0.011/2, # pew number adjusted for other errors
                 prior_rv_bias = 0.011, # RV and A bias from 538
                 prior_lv_bias = 0.0,
                 prior_a_bias = 0.02,
                 lv_rv_ratio = 5)

intent_model = get_intent_m(intent_model_path, opt$recompile)
intent_fit = sampling(intent_model, data=intent_d, chains=opt$chains, iter=opt$iter, warmup=500,
                      pars=c("sigma", "sigma_u", "sd_firm", "bias_rv", "bias_lv",
                             "bias_a", "house_effects", "natl_dem"),
                      control=list(adapt_delta=0.99, max_treedepth=15))

draws = recover_types(intent_fit) %>%
    spread_draws(natl_dem[week])

firm_ct = polls %>%
    group_by(firm) %>%
    summarize(n=n())
firms = recover_types(intent_fit, polls) %>%
    spread_draws(house_effects[firm]) %>%
    group_by(firm) %>%
    summarize(effect = median(house_effects)) %>%
    left_join(firm_ct, by="firm")

###
### Estimate distribution of seats
###

cat("Producing seat estimates.\n\n")

results_model = get_results_m(results_model_path, results, opt$recompile)

logit_pop = draws %>%
    filter(week==intent_d$W) %>%
    pull(natl_dem) %>%
    qlogis

inf_data = crossing(results_20, logit_pop=logit_pop)
seats = posterior_predict(results_model, newdata=inf_data, draws=10) %>%
    as.numeric %>%
    round

###
### Output predictions
###

entry = tibble(
    date = from_date,
    s_exp = median(seats),
    s_min = min(seats),
    s_q05 = quantile(seats, 0.05),
    s_q25 = quantile(seats, 0.25),
    s_q75 = quantile(seats, 0.75),
    s_q95 = quantile(seats, 0.95),
    s_max = max(seats),
    i_exp = plogis(median(logit_pop)),
    i_q05 = plogis(quantile(logit_pop, 0.05)),
    i_q25 = plogis(quantile(logit_pop, 0.25)),
    i_q75 = plogis(quantile(logit_pop, 0.75)),
    i_q95 = plogis(quantile(logit_pop, 0.95)),
    prob = mean(seats >= 218),
    prob_gain = mean(seats > current_seats),
    prob_pop = mean(logit_pop > 0)
)
entry = mutate_if(entry, is.numeric, ~ round(., 4))

intent_tbl = draws %>%
    group_by(week) %>%
    summarize(i_exp = median(natl_dem),
              i_q05 = quantile(natl_dem, 0.05),
              i_q25 = quantile(natl_dem, 0.25),
              i_q75 = quantile(natl_dem, 0.75),
              i_q95 = quantile(natl_dem, 0.95)) %>%
    mutate(date = election_day - 7*intent_d$W + 7*week)
output = append(as.list(entry), list(
    time = Sys.time(),
    n_polls = nrow(polls),
    gain = entry$s_exp - current_seats,
    intent = intent_tbl,
    firm_effects = firms,
    hist = map_int(entry$s_min:entry$s_max, ~ sum(seats == .))
))

with(output, cat(glue("
    ===========================================
     2020 U.S. House Forecast
     {as.character(Sys.Date(), format='%B %d, %Y')}
    -------------------------------------------
     Forecast from: {as.character(from_date, format='%B %d, %Y')}
     {round((election_day - from_date)/7)} weeks until the election.
     {nrow(polls)} polls.

     Dem. share of two-party vote:   {round(100*i_exp)}%
     Median seat estimate:           {round(s_exp)}
     Estimated seat range:           {round(s_q05)} - {round(s_q95)}
     Median seat gain:               {ifelse(gain>=0, '+', '-')}{abs(round(gain))}
     Probability of taking control:  {round(100*prob)}%
    ===========================================
    ")))
cat("\n\n")
system("osascript -e 'display notification \"Model run complete.\" with title \"House Model\"'")

if (opt$dry) quit("no")

if (file.exists(opt$history_file)) {
    history = bind_rows(
        suppressMessages(read_csv(opt$history_file)),
        entry
    ) %>%
        group_by(date) %>%
        slice(n()) %>%
        ungroup %>%
        arrange(date)
} else {
    history = entry
}
write_csv(history, opt$history_file, na="")

# only save full output if current run
if (from_date != Sys.Date()) quit("no")
write_json(output, opt$output_file, auto_unbox=T, digits=7)
