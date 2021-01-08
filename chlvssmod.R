# 1/23/2019: model TSS rather than nvss and vssn
# 1/31/2019: model TSS with changing coef with chl
## 12.18.2020: Model seston
## add diagnostics
chlvssmod <- function(df1, matout = NULL,varout = NULL,
                        runmod = T, xvalid = F) {

    df1$vss <- df1$VSS
    df1$chl <- df1$Chl

    incvec <- ! is.na(df1$vss) & ! is.na(df1$chl)
    df1 <- df1[incvec,]
    df1$lake <- factor(df1$lake)

    df1$lakenum <- as.numeric(df1$lake)

    varlist<- c("vss", "chl")
    mn.val <- apply(df1[, varlist],2,function(x) exp(mean(log(x))))
    print(mn.val)

    for (i in varlist) df1[,i] <- df1[,i]/mn.val[i]

    modstan <- '
        data {
            int n;
            int nlake;
            int lakenum[n];
            vector[n] chl;
            vector[n] vss;
        }
        parameters {
            real muu;
            real<lower = 0> sigu;
            vector[n] etau;

            real mub;
            real k;

            real<lower = 0> sigvss;

        }
        transformed parameters {
            vector[n] u;
            vector[n] vss_mn;

            u = muu + etau*sigu;
            for (i in 1:n) {
               vss_mn[i] = exp(mub)*chl[i]^k + exp(u[i]);
            }
        }
        model {
            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);
            mub ~ normal(0,3);
            k ~ normal(1,1);
            sigvss ~ cauchy(0,3);
            vss ~ lognormal(log(vss_mn), sigvss);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2)/length(x))

    extractvars <- c("k", "mub")

    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 1
        options(mc.cores = nchains)

        datstan <- list(n = nrow(df1),
                        nlake = max(df1$lakenum),lakenum = df1$lakenum,
                        vss = df1$vss,
                        chl = df1$chl)
        print(str(datstan))

        fit <- stan(model_code = modstan,
                    data = datstan, iter = 600, chains = nchains,
                    warmup = 300, thin= 2,
                    control = list(adapt_delta = 0.98, max_treedepth = 14))
        varout <- extract(fit, pars = extractvars)
        return(fit)
    }

    dev.new()

    mud <- apply(varout$mud, 2, mean)
    k <- apply(varout$k, 2, mean)
    u <- apply(varout$u, 2, mean)
    selvec <- log(df1$ntu) < -3

    psed <- exp(mud[1])*df1$nvss^k[2] + exp(mud[2])*exp(u)^k[3]
    psed.vss <- exp(mud[2])*exp(u)^k[3]
    plot(log(psed), psed.vss/psed)
    stop()
    mud.ntu <- apply(varntu.01$mud, 2, mean)
    k.ntu <- apply(varntu.01$k, 2, mean)
    u.ntu <- apply(varntu.01$u, 2, mean)

    psed2 <- exp(mud.ntu[1])*exp(u.ntu)^k.ntu[2]
    plot(log(psed)[!selvec], log(psed2))
    abline(0,1)
    stop()
    plot(log(df1$chl), u)
    require(mgcv)
    mod <- gam(u ~ s(log(chl), k = 4), data = df1)
    iord <- order(df1$chl)
    predout <- predict(mod)
    lines(log(df1$chl)[iord], predout[iord])

    stop()

    plot(log(df1$chl), log(df1$tp-df1$dtp))
    k <- apply(varout$k, 2, mean)
    mud <- apply(varout$mud, 2, mean)
    abline(mud[3], k[4])
    abline(mud.ntu[2], k.ntu[3], col = "red")
    stop()

    tppred <- gettp(df1, varout, withu = TRUE)
    print(summary(tppred))

    print(summary(matout))
    dev.new()
    plot(matout[,1], matout[,2])
    points(log(tppred), log(df1$tp), pch = 16, col = "red", cex = 0.6)
    abline(0,1)


    print(rmsout(log(tppred), log(df1$tp)))
    print(rmsout(matout[,1], matout[,2]))

    return()

}
chlvssmod(mo.all)
