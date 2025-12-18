data {
    int<lower=1> K; //number of years
    array[K] int n_k; //number of samples in year
    array[K] int y_k; //number of positive samples in year
}
parameters {
    real mu; 
    real<lower = 0> sigma; 
    real<lower = 0, upper = 1> rho;
    array[K] real eps;
}

transformed parameters {
    array[K] real mu_k;

    mu_k[1] = mu + sigma*eps[1];
    for(k in 2:K){
        mu_k[k] = mu*(1-rho)+rho*mu_k[k-1]+sigma*eps[k];
    }
}

model {
    mu ~ normal(0, 1);
    rho ~ beta(2, 2);
    eps ~ normal(0, 1);

    y_k ~ binomial_logit(n_k, mu_k);
}

generated quantities{
    array[K] real<lower = 0, upper = 1> p_k;

    for (k in 1:K) {
        p_k[k] = inv_logit(mu_k[k]);
    }
}