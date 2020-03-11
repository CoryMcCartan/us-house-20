# HOUSE MODEL DEVELOPMENT
# CORY McCARTAN

library(loo)
library(shinystan)
library(rstan)
library(magrittr)
library(gtools)
library(dplyr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

####################
#   VOTER INTENT   #
####################
polls = read.csv("model/data/current_polls.csv", colClasses=c(date="Date"))
model.data = list(W = 1 + max(polls$week),
                  P = max(polls$firm.id),
                  N = nrow(polls),
                  w = 1 + max(polls$week) - polls$week,
                  d = as.numeric((polls$date - min(polls$date)) / 7),
                  dp = 1:(1 + max(polls$week)),
                  p = polls$firm.id,
                  ldem = with(polls, logit(dem / (dem + gop))),
                  n_resp = polls$n_resp,
                  n_side = polls$n_side,
                  n_dem = polls$n_dem)
poll.w = max(model.data$w)


intent.model = stan(file="model/stan/gp-intent-model.stan", model_name="intent",
                    data=model.data, iter=1000, warmup=500, chains=1)#,
                    #control=list(adapt_delta=0.99, max_treedepth=15))

intent.model = stan_model(file="model/stan/gp-intent-model.stan")
model.results = sampling(intent.model, data=model.data, iter=100, chains=1)

# run intent model
intent.model = stan(file="model/stan/intent-model.stan", model_name="intent",
                    data=model.data, iter=5000, warmup=1000, chains=1,
                    control=list(adapt_delta=0.99, max_treedepth=15))
# extract samples and estimates
print(intent.model, pars=c("sd_poll", "sd_walk", "sd_pollster", "nu", "rho", "mu_pollster"))
samples = rstan::extract(intent.model, pars=c("dem_margin", "logit_dem"))

# get mean and sd of expected voter intent
logit.est = mean(samples$logit_dem[,model.data$W])
logit.sd = sd(samples$logit_dem[,model.data$W])

# plot dem. support estimates
mrg = as.data.frame(t(apply(samples$dem_margin, 2, quantile, 
                            probs=c(0.05, 0.5, 0.95))))
names(mrg) = c("low", "median", "high")
mrg$week = 1:nrow(mrg)

ggplot(mrg, aes(x=week)) + geom_line(aes(y=median)) + 
    geom_ribbon(aes(ymin=low, ymax=high), alpha=0.5) + 
    coord_cartesian(ylim=c(-0.05, 0.20)) + 
    geom_hline(yintercept=0, linetype="dashed") +
    geom_vline(xintercept=poll.w)


####################
#  RESULTS MODEL   #
####################
form = ~ midterm*pres + pres:appr + pres:earn + pres:unemp + midterm*I(before-218)

# set up data list input
model.data = read.csv("data/combined.csv") %>% filter(weeks_until == 0)
data.list = as.list(model.data)
data.list$X = model.matrix(form, data=model.data)
data.list$K = ncol(data.list$X)
data.list$year = with(data.list, (year - min(year))/2 + 1)
data.list$N = length(data.list$year)
data.list$Y = max(data.list$year)
# use estimated s.d. intent from this year, but increased slightly
data.list$sd_intent = rep(1.5*sd(samples$logit_dem[,poll.w]), data.list$N)

# run results model
results.model = stan(file="stan/results-model.stan", model_name="results",
                     data=data.list, iter=11000, warmup=1000, chains=1,
                     control=list(adapt_delta=0.999, max_treedepth=10))
# output estimates
print(results.model, pars=c("beta_intent", "betas", "error"))
plot(results.model, pars=c("beta_intent", "betas", "error"))
# extract and save samples
est = rstan::extract(results.model)
saveRDS(est, "results-model-samples.rds")
samples_pred = est$seats_pred
seats_pred = apply(samples_pred, 2, median)
resid = model.data$seats - seats_pred

# plot actual vs predicted seats for each year
ggplot(model.data) + 
    geom_point(aes(x=factor(year), y=seats_pred)) + 
    geom_point(aes(x=factor(year), y=seats), color='red')

ggplot(model.data) + 
    geom_point(aes(x=factor(year), y=resid))

# calculate in-sample predicted probabilities of control
res = model.data %>% 
    transmute(year=year, 
              control=ifelse(seats >= 218, 1, 0),
              seats=seats,
              deficit=seats - 218,
              pr_control=0)

for (i in 1:nrow(res)) {
    res$pr_control[i] = mean(samples_pred[,i] >= 218)
}

brier = mean((res$pr_control - res$control)^2)
brier_1 = mean((1/3 - res$control)^2) # always guess 0.333
brier_skill = 1 - brier / brier_1

# test 2018 predicitons
election.d = data.frame(logit_intent=logit.est, sd_intent=logit.sd, 
                        midterm=1, pres=-1, house=-1, appr=log(0.37/(1-0.37)), 
                        unemp=0.043, earn=0.0240, before=194)

post_pred = function(est, new.d) {
    X = model.matrix(form, data=new.d)
    N = length(est$error)
    logit_true = rnorm(N, new.d$logit_intent, new.d$sd_intent)
    mu = est$beta_intent*logit_true + c(est$betas %*% t(X))
    rnorm(N, mu, est$error) + new.d$before
}
                        
# histogram of expected seats 
pr = post_pred(est, election.d)
print(mean(pr >= 218))
ggplot() + 
    geom_histogram(aes(pr, fill=(pr >= 218)), binwidth=3, center=0.5) + 
    scale_fill_manual(values=c("FALSE"="red", "TRUE"="blue"), 
                      labels=c("FALSE"="GOP", "TRUE"="DEM")) +
    geom_vline(xintercept=election.d$before, linetype="dashed") + 
    geom_vline(xintercept=203) + 
    geom_vline(xintercept=193) + 
    xlim(0, 435) + labs(x="Seats", y="", fill="Control")


# calculate how probability of control changes w.r.t. some input variable
probs = data.frame(weeks=1:57, low=0, med=0, high=0, prob=0)
for (i in 1:nrow(probs)) {
    logit.est = mean(samples$logit_dem[,i])
    logit.sd = sd(samples$logit_dem[,i])
    election.d = data.frame(logit_intent=logit.est, sd_intent=logit.sd, 
                            midterm=1, pres=-1, house=-1, appr=log(0.37/(1-0.37)), 
                            unemp=0.043, earn=0.0240, before=194)
    pr = post_pred(est, election.d)
    
    q = quantile(pr, c(0.1, 0.5, 0.9))
    probs$low[i] = q[1]
    probs$med[i] = q[2]
    probs$high[i] = q[3]
    probs$prob[i] = mean(pr >= 218)
}

ggplot(probs, aes(x=weeks)) + geom_line(aes(y=med)) + 
    geom_ribbon(aes(ymin=low, ymax=high), alpha=0.5) +
    geom_hline(yintercept=218, linetype="dashed")

ggplot(probs, aes(x=weeks, y=prob)) + geom_line()