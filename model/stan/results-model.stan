data {
    int N; // number of observations
    int K; // number of addl covariates
    
    vector[N] logit_intent; // logit of voter intent
    vector<lower=0>[N] sd_intent; // measurement error in intent
    vector[N] seats; // response variable
    vector[N] before;  // seats won in previous election
    matrix[N,K] X; // addl. covariates
}

transformed data {
    // QR decomposition
    matrix[N,K] Q = qr_Q(X)[, 1:K] * N;
    matrix[K,K] R = qr_R(X)[1:K, ] / N;
    matrix[K,K] R_inv = inverse(R);
}

parameters {
    vector[N] logit_true; // true voter intent, before meas. error
    // model coefficients 
    real beta_intent;
    vector[K] thetas; // (transformed for QR decomp)
    
    real<lower=0> error; // model error
}

transformed parameters {
    vector[K] betas = R_inv * thetas;
    vector[N] mu = beta_intent*logit_true + Q*thetas;
}

model {
    // measurement error in logit-intentions
    logit_true ~ normal(0, 1);
    logit_intent ~ normal(logit_true, sd_intent);
    
    seats - before ~ normal(mu, error);
    
    beta_intent ~ student_t(3, 0, 10);
    betas ~ student_t(3, 0, 10);
    
    error ~ cauchy(0, 20);
}

generated quantities {
    vector[N] log_lik;
    vector[N] seats_pred;
    
    for (i in 1:N) {
        log_lik[i] = normal_lpdf(seats[i]-before[i] | mu[i], error);
        seats_pred[i] = before[i] + normal_rng(mu[i], error);
    }
}
