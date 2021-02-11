## simulate data and see if stan model gets it right

sim <- function(runmod, varout = NULL) {

    set.seed(13)
    n <- 50
    b <- c(1,1.5)
    ind <- as.numeric(runif(n) > 0.5) + 1
    k <- 0.7
    muu <- 1
    sdu <- 2
    sdy <- 0.1
    u <- rnorm(n, mean = muu, sd = sdu)
    x <- rnorm(n, mean = 0, sd = 1)
    err <- rnorm(n, mean = 0, sd = sdy)
    y <- b[1]*exp(x)^k + exp(u)
    yobs <- log(y) + err

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2))
    plot(x, yobs)
    abline(log(b[1]), k)


    beta <- exp(u)/(b[1]*exp(x)^k)
    incvec <- rep(T, times = n)

#    incvec <- x > 0.6

    points(x[incvec], yobs[incvec],pch = 16)

    plot(x, exp(u)/(b[1]*exp(x)^k))

    datstan <- list(n = sum(incvec), yobs = yobs[incvec],
                    x = x[incvec], ind = ind[incvec])
    print(str(datstan))


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
          real b;
          real k;
          real<lower = 0> sigy;
       }
       transformed parameters {
          vector[n] u;
          vector[n] ymean;
          u = muu + sigu*etau;
          for (i in 1:n) {
             ymean[i] = log_sum_exp(k*x[i] + b,u[i]);
          }
       }
       model {
          muu ~ normal(0,3);
          sigu ~ cauchy(0,3);
          etau ~ normal(0,1);
          sigy ~ normal(0.1, 0.02);
          b ~ normal(0,3);
          k ~ normal(1,1);
          yobs ~ normal(ymean, sigy);
       }
    '
    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 3
        options(mc.cores = nchains)
        fit <- stan(model_code = modstan,
                    data = datstan, iter = 4000, chains = nchains,
                    warmup = 1000, thin= 3,
                    control = list(adapt_delta = 0.98, max_treedepth = 14))
        return(fit)
    }
    else {
        umod <- apply(varout$u, 2, mean)
        bmod <- mean(varout$b)
        kmod <- mean(varout$k)
        predout <- exp(bmod)*exp(x[incvec])^kmod + exp(umod)

        dev.new()
        par(mar = c(4,4,1,1), mfrow = c(2,2))
        plot(x[incvec], u[incvec]-umod)
        abline(h=0)
        plot(x[incvec], log(b[1]) + k*x[incvec] -
            (bmod + kmod*x[incvec]))
        abline(h=0)
#        u2 <- apply(varout2$u, 2, mean)
#        b2 <- apply(varout2$b, 2, mean)
#        k2 <- mean(varout2$k)
#        predout2 <- exp(b2[ind])*exp(x)^k2 + exp(u2)


        plot(x[incvec], yobs[incvec] - log(predout))
#        points(yobs, log(predout2) - yobs, pch = 16)
        abline(h = 0)
        plot(log(predout), yobs[incvec])
        abline(0,1)
        stop()
        rmsout <- function(x, y) sqrt(sum((x-y)^2)/length(x))
        print(rmsout(log(predout), yobs))
#        print(rmsout(log(predout2), yobs))
        stop()


    }

}

fitout <- sim(runmod = T)
print(fitout, par = c("b", "k", "muu", "sigu"))
varout <- extract(fitout, pars = c("b", "k", "u"))
sim(runmod = F, varout = varout)
