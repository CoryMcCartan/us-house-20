data {
    int W; // number of weeks 
    int N; // number of polls
    int P; // number of polling firms
    
    int<lower=1> w[N]; // the week of a given poll
    int<lower=1> n_resp[N]; // the size a given poll
    int<lower=1> n_side[N]; // respondents who pick either DEM or GOP
    int<lower=1> n_dem[N]; // respondents who pick DEM
    int<lower=1> p[N]; // the polling firm for a given poll
}

parameters {
    vector[W] delta_dem; // steps of random walk
    real<lower=0> sd_walk; // week-to-week variance
    real<lower=1> nu; // df for random walk
    real<lower=0,upper=1> rho; // strength of mean reversion
    
    vector[P] RP_pollster;  
    real<lower=0> sd_pollster; // hyperparameter for pollster errors
    real mu_pollster; // hyperparameter for global polling error
    
    real<lower=0,upper=1> prop_undecided; 
    
    vector[N] RP_poll;  
    real<lower=0> sd_poll; // hyperparameter for polling errors
}

transformed parameters {
    vector[W] logit_dem; // democratic support by week (main param)
    vector[N] logit_poll; // support for dem. in specific poll
    
    vector[P] pollster_error;  
    vector[N] poll_error; 
    
    pollster_error = mu_pollster + RP_pollster * sd_pollster;
    poll_error = RP_poll * sd_poll;
    
    logit_dem[1] = 10 * delta_dem[1]; // equiv. to t(nu, 0, 10) prior on logit_dem[1]
    for (i in 2:W)
        logit_dem[i] = rho*logit_dem[i-1] + sd_walk*delta_dem[i];
    
    for (i in 1:N)
        logit_poll[i] = logit_dem[w[i]] + pollster_error[p[i]] + poll_error[i];
}

model {
    n_side ~ binomial(n_resp, 1 - prop_undecided);
    n_dem ~ binomial_logit(n_side, logit_poll);
    
    delta_dem ~ student_t(nu, 0, 1);
    
    RP_pollster ~ student_t(2, 0, 1);
    mu_pollster ~ normal(0, 0.02);
    RP_poll ~ student_t(2, 0, 1);
    prop_undecided ~ beta(2, 2);
    
    nu ~ gamma(2, 0.1);
    
    sd_walk ~ cauchy(0, 10);
    rho ~ beta(2, 1);
    sd_pollster ~ student_t(4, 0, 1);
    sd_poll ~ student_t(4, 0, 1);
}

generated quantities {
    vector[N] log_lik;
    vector[W] dem_margin; 
    
    for (i in 1:N) {
        log_lik[i] = binomial_lpmf(n_side[i] | n_resp[i], 1 - prop_undecided) +
                        binomial_logit_lpmf(n_dem[i] | n_side[i], logit_poll[i]);
    }
    
    dem_margin = 2*inv_logit(logit_dem) - 1;
}
