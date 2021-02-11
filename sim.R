## simulate data and see if stan model gets it right

sim <- function() {

    set.seed(10)
    n <- 200
    b <- c(1,1.5)
    ind <- as.numeric(runif(n) > 0.5) + 1
    print(ind)
    k <- 0.7
    muu <- 1
    sdu <- 1
    sdy <- 0.1
    u <- rnorm(n, mean = muu, sd = sdu)
    x <- rnorm(n, mean = 1.5, sd = 1)
    err <- rnorm(n, mean = 0, sd = sdy)
    y <- b[ind]*exp(x)^k + exp(u)
    yobs <- log(y) + err

    plot(x, yobs)
    points(x[ind==2], yobs[ind == 2], pch = 16)
    abline(0.08, 0.66, lty = "dashed")
    abline(0, 0.7)

    datstan <- list(n = n, yobs = yobs, x = exp(x), ind = ind)

    modstan <- '
       data {
          int n;
          vector[n] yobs;
          vector[n] x;
          int ind[n];
       }
       parameters {
          real muu;
          real<lower = 0> sigu;
          vector[n] etau;
          vector[2] b;
          real k;
          real<lower = 0> sigy;
       }
       transformed parameters {
          vector[n] u;
          vector[n] ymean;
          u = muu + sigu*etau;
          for (i in 1:n) {
             ymean[i] = exp(b[ind[i]])*x[i]^k + exp(u[i]);
          }
       }
       model {
          muu ~ normal(0,3);
          sigu ~ cauchy(0,3);
          etau ~ normal(0,1);
          sigy ~ cauchy(0,3);
          b ~ normal(0,3);
          k ~ normal(1,1);
          yobs ~ student_t(4,log(ymean), sigy);
       }
    '
    require(rstan)
    rstan_options(auto_write = TRUE)
    nchains <- 3
    options(mc.cores = nchains)
    fit <- stan(model_code = modstan,
                data = datstan, iter = 1800, chains = nchains,
                warmup = 600, thin= 2)
    return(fit)

}

fitout <- sim()
