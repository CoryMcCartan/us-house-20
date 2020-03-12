/************************************************
 * U.S. HOUSE ELECTION POLLING MODEL            *
 * CORY McCARTAN                                *
 * (c) 2020                                     *
 ************************************************/

data {
    int N; // number of polls
    int W; // number of weeks
    int N_firm; // number of polling firms

    int<lower=1, upper=W> week[N]; // the week of a given poll
    int<lower=1, upper=N_firm> firm[N]; // the polling firm for a given poll
    vector<lower=0, upper=1>[N] type_rv; // RV indicator
    vector<lower=0, upper=1>[N] type_lv; // LV indicator
    vector<lower=0, upper=1>[N] type_a; // adult indicator
    vector[N] dem;
    vector[N] var_poll;

    real prior_final_mean;
    real<lower=0> prior_final_sd;
    real prior_rv_bias;
    real prior_lv_bias;
    real prior_a_bias;
    real<lower=0> prior_natl_poll_error;
    real<lower=0> lv_rv_ratio;
}

transformed data {
    vector[N] logit_dem = logit(dem);
    vector[N] sd_poll = sqrt(var_poll);
}

parameters {
    real<lower=0> sigma; // week-to-week variance

    vector[N] u; // poll-specific error
    real<lower=0> sigma_u;  // add'l polling error (beyond sampling)
    real<lower=0> sd_firm; // hyperparameter for pollster errors
    real natl_error;
    real bias_rv; // registered voter poll bias
    real bias_a; // adult voter poll bias
    real bias_lv; // likely voter poll bias

    // reparametrization
    vector[W] delta_dem; // steps of random walk
    vector[N_firm] delta_firm; // polling firm standardized errors
}

transformed parameters {
    vector[W] mu; // national voter intention, logit
    vector[N] val_poll; // support for dem. in specific poll

    mu[W] = prior_final_mean + 4*prior_final_sd*delta_dem[W];
    for (i in 1:(W-1)) {
        int w = W - i;
        mu[w] = mu[w+1] + sigma*delta_dem[w];
    }

    val_poll = mu[week] + 4*sd_firm*delta_firm[firm] + sigma_u*u + natl_error
        + 4*bias_rv*type_rv + 4*bias_lv*type_lv + 4*bias_a*type_a;
}

model {
    logit_dem ~ normal(val_poll, sd_poll);

    u ~ std_normal();
    // reparametrizations
    delta_dem ~ std_normal();
    delta_firm ~ std_normal();

    sigma ~ gamma(2, 2/0.02);
    sigma_u ~ exponential(1/0.02);
    natl_error ~ normal(0, prior_natl_poll_error);
    sd_firm ~ exponential(1/0.05);
    bias_rv ~ normal(prior_rv_bias, 0.03);
    bias_a ~ normal(prior_a_bias, 0.03);
    bias_lv ~ normal(prior_lv_bias, 0.03/lv_rv_ratio);
}

generated quantities {
    vector[N] log_lik;
    vector[W] natl_dem = inv_logit(mu);
    vector[N_firm] house_effects = sd_firm * delta_firm;

    for (i in 1:N) {
        log_lik[i] = normal_lpdf(logit_dem[i] | val_poll[i], sd_poll[i]);
    }
}
