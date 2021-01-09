# 1/23/2019: model TSS rather than nvss and vssn
# 1/31/2019: model TSS with changing coef with chl
## 12.18.2020: Model seston
## add diagnostics
chlvssmod <- function(df1,varout = NULL,
                        runmod = T) {

    df1$vss <- df1$VSS
    df1$chl <- df1$Chl

    incvec <- ! is.na(df1$vss) & ! is.na(df1$chl)
    df1 <- df1[incvec,]
    df1$lake <- factor(df1$lake)
    df1$lakenum <- as.numeric(df1$lake)

    ## select data starting with ~300 samples per year
    incvec <- df1$year >= 1989
    df1 <- df1[incvec,]

    ## select at least 20 samples per lake
    print(nrow(df1))
    nsamp <- table(df1$lake)
    idsav <- names(nsamp)[nsamp > 20]
    incvec <- rep(F, times = nrow(df1))
    for (i in idsav) incvec <- incvec | i == df1$lake
    df1 <- df1[incvec,]
    print(nrow(df1))

    ## reset lake ids
    df1$lake <- factor(df1$lake)
    df1$lakenum <- as.numeric(df1$lake)

    ## scale variables
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
            real<lower = 0> sigu[2];  //two levels for VSSnp
            vector[nlake] etau1;
            vector[n] etau2;

            real mub;
            real k;

            real<lower = 0> sigvss;

        }
        transformed parameters {
            vector[nlake] u1;
            vector[n] u2;
            vector[n] vss_mn;

            u1 = muu + etau1*sigu[1];
            u2 = u1[lakenum] + etau2*sigu[2];

            for (i in 1:n) {
               vss_mn[i] = exp(mub)*chl[i]^k + exp(u2[i]);
            }
        }
        model {
            muu ~ normal(0,3);
            etau1 ~ normal(0,1);
            etau2 ~ normal(0,1);
            sigu ~ cauchy(0,3);
            mub ~ normal(0,3);
            k ~ normal(1,1);
            sigvss ~ cauchy(0,3);
            vss ~ lognormal(log(vss_mn), sigvss);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2)/length(x))

    extractvars <- c("k", "mub", "u1")

    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 3
        options(mc.cores = nchains)

        datstan <- list(n = nrow(df1),
                        nlake = max(df1$lakenum),lakenum = df1$lakenum,
                        vss = df1$vss,
                        chl = df1$chl)
        print(str(datstan))

        fit <- stan(model_code = modstan,
                    data = datstan, iter = 800, chains = nchains,
                    warmup = 300, thin= 2,
                    control = list(adapt_delta = 0.98, max_treedepth = 14))
        return(fit)

    }

    k <- mean(varout$k)
    mub <- mean(varout$mub)
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2))
    plot(log(df1$chl), log(df1$vss), col = "grey")
    abline(mub, k)

    df2 <- merge(df1, res.dat, by.x = "lake", by.y = "MU..")
    df2 <- unique.data.frame(df2[, c("lake", "lakenum", "flush.rate", "mean.depth")])
    df2 <- df2[order(df2$lakenum),]
    u1 <- apply(varout$u1, 2, mean)
    plot(log(df2$mean.depth), u1)

    stop()

    return()

}
#fitout <- chlvssmod(mo.all, runmod = T)
chlvssmod(mo.all, varout, runmod= F)
