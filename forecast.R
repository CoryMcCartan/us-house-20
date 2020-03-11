#!/usr/bin/env Rscript

#####################################
#     U.S. HOUSE MODEL              #
#     CORY McCARTAN                 #
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
    make_option("--iter", type="integer", default=1000,
                help="Number of MCMC iterations for voter intent estimation."),
    make_option("--recompile", action="store_true", default=F,
            help="Force recompile of STAN model."),
    make_option("--model_dir", type="character", default="model/stan",
                help="The directory in which the models are stored"),
    make_option("--samples_file", type="character", default="model/results-model-samples.rds",
                help="The path to the posterior samples from the results model."),
    make_option("--output_file", type="character", default="docs/data/output.json",
                help="The file to save estimates to."),
    make_option("--history_file", type="character", default="docs/data/history.csv",
                help="The file where past estimates are saved.")
)
opt = parse_args(OptionParser(option_list=option_list, 
                              description="Forecast 2018 U.S. House election."))

from.date = as.Date(opt$date)
intent.model.path = file.path(opt$model_dir, "intent-model")
results.model.path = file.path(opt$model_dir, "results-model")


suppressMessages(library(magrittr))
suppressMessages(library(dplyr))
suppressMessages(library(glue))
suppressMessages(library(jsonlite))
suppressMessages(library(pollstR))
suppressMessages(library(lubridate))
suppressMessages(library(rstan))

current.seats = 194
election.day = ymd("2018-11-06")


###
### Download new polling data 
###
cat("Downloading polling data.\n")

# generic ballot polls
polls = suppressMessages(pollster_charts_polls("2018-national-house-race"))[["content"]] %>%
    filter(partisanship == "Nonpartisan") %>%
    transmute(dem = Democrat / 100,
              gop = Republican / 100,
              other = Other / 100,
              undecided = Undecided / 100,
              n_resp = observations,
              n_side = round((dem + gop) * n_resp),
              n_dem = round(dem * n_resp),
              pollster = survey_house,
              firm.id = as.numeric(as.factor(pollster)),
              date = ymd(end_date),
              type = recode(sample_subpopulation, `Registered Voters`="RV",
                            `Likely Votesrs`="LV")) %>%
    filter(date <= from.date) %>%
    as.data.frame

#polls = read.csv("generic_polllist.csv") %>% 
#    transmute(dem = dem / 100,
#              gop = rep / 100,
#              n_resp = as.integer(samplesize),
#              n_side = round((dem + gop) * n_resp),
#              n_dem = round(dem * n_resp),
#              pollster = pollster,
#              firm.id = as.numeric(as.factor(pollster)),
#              date = mdy(enddate),
#              type = population) %>%
#    filter(date <= from.date) %>% 
#    as.data.frame

# drop missing data
polls = polls[complete.cases(polls),]

# week IDs
start.day = min(polls$date)
n.weeks = ceiling(as.numeric(election.day - start.day) / 7)
polls$week = floor(as.numeric(election.day - polls$date) / 7)

write.csv(polls, "model/data/current_polls.csv", row.names=F) # save for reference

# presidential approval rating
appr = (suppressMessages(pollster_charts_trendlines("trump-job-approval"))[["content"]] %>%
    filter(label == "Approve", date <= from.date) %>%
    tail(1))$value
logit.appr = log((appr/100) / (1 - appr/100))



###
### Estimate voter intent
###

cat("Estimating voter intent.\n")

compiled.model = paste0(intent.model.path, ".rds")
if (opt$recompile && file.exists(compiled.model))
    file.remove(compiled.model)

model.data = list(W = 1 + max(polls$week),
                  P = max(polls$firm.id),
                  N = nrow(polls),
                  w = 1 + max(polls$week) - polls$week,
                  p = polls$firm.id,
                  n_resp = polls$n_resp,
                  n_side = polls$n_side,
                  n_dem = polls$n_dem)
poll.w = max(model.data$w)

# run intent model
if (file.exists(compiled.model)) {
    model.obj = readRDS(compiled.model)
} else {
    model.obj = stan_model(file=paste0(intent.model.path, ".stan"),
                           model_name="intent")
    saveRDS(model.obj, compiled.model)
}
intent.model = sampling(model.obj, data=model.data, 
                        iter=opt$iter+1000, warmup=1000, chains=1,
                        show_messages=F, refresh=-1,
                        control=list(adapt_delta=0.99, max_treedepth=15))

#int.mod.2 = optimizing(model.obj, data=model.data, as_vector=F)

# extract samples and estimates
samples = rstan::extract(intent.model, pars=c("dem_margin", "logit_dem"))
#logit.est = int.mod.2$par$logit_dem[model.data$W]

# get mean and sd of expected voter intent
logit.est = mean(samples$logit_dem[,model.data$W])
logit.sd = sd(samples$logit_dem[,model.data$W])
 


###
### Estimate distribution of seats
###

cat("Producing seat estimates.\n\n")

if (file.exists(opt$samples_file)) {
    est = readRDS(opt$samples_file)
} else {
    if (opt$recompile && file.exists(paste(results.model.path, ".rds")))
        file.remove(paste0(results.model.path, ".rds"))

    form = ~ midterm*pres + pres:appr + pres:earn + pres:unemp + midterm*I(before-218)

    model.data = read.csv("model/data/combined.csv") %>% filter(weeks_until == 0)
    data.list = as.list(model.data)
    data.list$X = model.matrix(form, data=model.data)
    data.list$K = ncol(data.list$X)
    data.list$year = with(data.list, (year - min(year))/2 + 1)
    data.list$N = length(data.list$year)
    # use estimated s.d. intent from this year, but increased slightly
    data.list$sd_intent = rep(1.5*sd(samples$logit_dem[,poll.w]), data.list$N)

    
    # run results model
    results.model = stan(file=paste0(results.model.path, ".stan"), model_name="results",
                        data=data.list, iter=11000, warmup=1000, chains=1,
                        control=list(adapt_delta=0.999, max_treedepth=10))
    
    est = rstan::extract(results.model)
    saveRDS(est, opt$samples_file)
}

election.d = data.frame(logit_intent=logit.est, sd_intent=logit.sd, 
                        midterm=1, pres=-1, house=-1, appr=logit.appr, 
                        unemp=0.043, earn=0.0240, before=current.seats)

post_pred = function(est, new.d) {
    form = ~ midterm*pres + pres:appr + pres:earn + pres:unemp + midterm*I(before-218)
    X = model.matrix(form, data=new.d)
    N = length(est$error)
    logit_true = sample(samples$logit_dem[,model.data$W], N, replace=T)
    #rnorm(N, new.d$logit_intent, new.d$sd_intent)
    mu = est$beta_intent*logit_true + c(est$betas %*% t(X))
    rnorm(N, mu, est$error) + new.d$before
}

seats = post_pred(est, election.d)

s.expected = median(seats)
s.min = quantile(seats, 0.05)
s.max = quantile(seats, 0.95)
gain = s.expected - current.seats
prob = mean(seats >= 218)
prob.gain = mean(seats >= current.seats) # dem. gain
prob.popv = pnorm(0, mean=logit.est, sd=logit.sd, lower.tail=F) # win pop. vote
seats.dist = hist(seats, br=0:436, plot=F)$density
# voter intent
mrg = as.data.frame(t(apply(samples$dem_margin, 2, quantile, 
                            probs=c(0.05, 0.50, 0.95))))
names(mrg) = c("low", "median", "high")
mrg$week = election.day - 7*(nrow(mrg):1 - 1)

###
### Output predictions
###

cat(glue("
    ===========================================
     2018 U.S. House Forecast
     {as.character(Sys.Date(), format='%B %d, %Y')}
    -------------------------------------------
     Forecast from: {as.character(from.date, format='%B %d, %Y')}
     {round((election.day - from.date)/7)} weeks until the election.
     {nrow(polls)} polls.

     Dem. margin in generic ballot:  {round(100*tail(mrg, 1)$median, 1)}%
     Median seat estimate:           {round(s.expected)}
     Estimated seat range:           {round(s.min)} - {round(s.max)}
     Median seat gain:               {ifelse(gain>=0, '+', '-')}{abs(round(gain))}
     Probability of taking control:  {round(100*prob, 1)}%
    ===========================================
    "))
cat("\n")


if (opt$dry) quit("no")

# history
history = data.frame(date=as.Date(character()), prob=numeric(), gain=numeric())
if (file.exists(opt$history_file)) {
    history = read.csv(opt$history_file, 
                       colClasses=c(date="Date", prob="numeric", gain="numeric"))
}
history = rbind(history, data.frame(date=from.date, prob=prob, gain=gain))
write.csv(history, opt$history_file, row.names=F)

if (from.date != Sys.Date()) quit("no") # only save full output if current run

# snapshot
output.data = list(
    date = from.date,
    time = Sys.time(),
    prob = prob,
    gain = gain,
    seats = s.expected,
    seats_min = s.min,
    seats_max = s.max,
    seats_dist = seats.dist,
    intent = mrg
)
write(toJSON(output.data, auto_unbox=T), opt$output_file)
