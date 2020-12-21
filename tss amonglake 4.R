# 1/23/2019: model TSS rather than nvss and vssn
# 1/31/2019: model TSS with changing coef with chl
## 12.18.2020: Model seston
tss.explore <- function(df1, varout = NULL, runmod = T) {

    df1$vss <- as.numeric(as.character(df1$vss))
    df1$dtp <- as.numeric(as.character(df1$dtp))
    df1$nvss <- as.numeric(as.character(df1$nvss))
    df1$tp <- as.numeric(as.character(df1$tp))
    df1$chl <- as.numeric(as.character(df1$chl))

    date0 <- strptime(paste(df1$month, df1$day, "2004", sep = "-"),
                      format = "%m-%d-%Y")
    df1$yday <- date0$yday

    print(nrow(df1))
    df1 <- merge(df1, res.dat[, c("MU..", "flush.rate")], by.x = "lake",
                 by.y = "MU..")
    print(nrow(df1))
    print(summary(df1$flush.rate))

    nperiod <- 4
#    cutp <- seq(0, 365, length = nperiod)
#    df1$yday.q <- cut(df1$yday, cutp,include.lowest = T)
    df1$yday.q <- 1
    incvec <- df1$month >= 3 & df1$month <= 5
    df1$yday.q[incvec]<- 2
    incvec <- df1$month >= 6 & df1$month <= 8
    df1$yday.q[incvec]<- 3
    incvec <- df1$month >= 9 & df1$month <= 11
    df1$yday.q[incvec]<- 4
    df1$yday.q <- factor(df1$yday.q)

    incvec <- ! is.na(df1$vss) & ! is.na(df1$chl) & ! is.na(df1$tp) &
        ! is.na(df1$nvss)
    df1 <- df1[incvec,]
    df1$lake <- factor(df1$lake)

    df1$lake <- factor(df1$lake)
    print(table(df1$lake))
    df1$lakenum <- as.numeric(df1$lake)

    df1$seasnum <- as.numeric(df1$yday.q)

    varlist<- c("nvss", "chl", "tp", "vss")
    mn.val <- apply(df1[, varlist],2,function(x) exp(mean(log(x))))
    print(mn.val)
    for (i in varlist) df1[,i] <- df1[,i]/mn.val[i]
    df1$dtp <- df1$dtp/mn.val["tp"]

    plot(log(df1$chl), log(df1$tp - df1$dtp))
    abline(-1.24, 0.89)

    datstan <- list(n = nrow(df1),
                    nlake = max(df1$lakenum),lakenum = df1$lakenum,
                    nseas = nperiod, seasnum = df1$seasnum,
                    nvss = df1$nvss,
                    tp = df1$tp,
                    dtp = df1$dtp,
                    vss = df1$vss,
                    chl = df1$chl)
    print(str(datstan))

    modstan <- '
        data {
            int n;
            int nlake;
            int lakenum[n];
            vector[n] tp;
            vector[n] dtp;
            vector[n] chl;
            vector[n] nvss;
            vector[n] vss;
        }
        parameters {
            real muu;
            real<lower = 0> sigu;
            vector[n] etau;

            real mub;
//            real<lower = 0> sigb;
//            vector[nseas] etab;

            real k[2];

            vector[3] mud;
            vector<lower = 0>[3] sigd;
            matrix[nlake,3] etad;
//            vector[nlake] etad2;
//            vector[nlake] etad3;
//            vector[nseas] etad4;

            real<lower = 0> sigvss;
            real<lower = 0> sigtp;

        }
        transformed parameters {
            vector[n] u;
            vector[n] tp_mn;
            vector[n] vss_mn;
            matrix[nlake, 3] d;

            u = muu + etau*sigu;

            for (i in 1:3) d[,i] = mud[i] + etad[,i]*sigd[i];

//            b = mub + etab*sigb;

//            d1 = mud[1] + sigd[1]*etad1;
//            d2 = mud[2] + sigd[2]*etad2;
//            for (i in 1:nlake) {
//                for (j in 1:nseas) {
//                    d3[i,j] = mud[3] + sigd[3]*etad3[i] + sigd[4]*etad4[j];
//                }
//            }
//            d3 = mud[3] + sigd[3]*etad3

            for (i in 1:n) {
               vss_mn[i] = exp(mub)*chl[i]^k[1] + exp(u[i]);

               tp_mn[i] = exp(d[lakenum[i],1])*nvss[i] +
                           exp(d[lakenum[i],2])*exp(u[i]) +
                           exp(d[lakenum[i],3])*chl[i]^k[2] + dtp[i];
            }
        }
        model {
            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);

            mub ~ normal(0,3);

            mud ~ normal(0,3);
            sigd ~ cauchy(0,3);
            for (i in 1:3) etad[,i] ~ normal(0,1);

            k[1] ~ normal(0.807,0.012);
            k[2] ~ normal(1,1);

            sigtp ~ cauchy(0,3);
            sigvss ~ cauchy(0,3);

            vss ~ lognormal(log(vss_mn), sigvss);
            tp ~ lognormal(log(tp_mn), sigtp);
        }
    '


    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)

        nchains <- 1
        options(mc.cores = nchains)


        fit <- stan(model_code = modstan,
                    data = datstan, iter = 1000, chains = nchains,
                    warmup = 500, thin= 1,
                    control = list(adapt_delta = 0.98, max_treedepth = 14))
        return(fit)
    }

#    varout <- extract(fit, pars = c("b", "d1", "d2", "d3", "u", "k", "sigd",
#                                    "mud"))

    df1$u <- apply(varout$u, 2, mean)
    b <- apply(varout$b, 2, mean)
    d1 <- apply(varout$d1, 2, mean)
    d2 <- apply(varout$d2, 2, mean)
    d3 <- apply(varout$d3, c(2,3), mean)
    k <- apply(varout$k, 2, mean)

    df1$dtp <- as.numeric(as.character(df1$dtp))
    df1$nvss <- as.numeric(as.character(df1$nvss))

    tsspred <- exp(b[df1$seasnum])*df1$chl.sc^k[1] + exp(df1$u)
    tppred <- rep(NA, times = nrow(df1))
    for (i in 1:nrow(df1)) {
        tppred[i] <- exp(d1[df1$lakenum[i]]) + exp(d2[df1$lakenum[i]])*exp(df1$u[i]) +
            exp(d3[df1$lakenum[i], df1$seasnum[i]])*df1$chl.sc[i]^k[2]
    }


    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0))
    plot(log(tsspred), log(df1$tss.sc))
    abline(0,1)
    plot(log(tppred), log(df1$tp.sc))
    abline(0,1)
    rmsout <- function(x,y) sqrt(sum((x-y)^2)/length(x))
    print(rmsout(log(tsspred), log(df1$tss.sc)))
    print(rmsout(log(tppred), log(df1$tp.sc)))
    return(varout)
    stop()

    dev.new()
    meanload <- rep(NA, times = 15)
    par(mar = c(4,4,1,1), mfrow = c(4,4), mgp = c(2.3,1,0))
    dp.all <- numeric(0)
    dchl.all <- numeric(0)
    for (i in 1:15){
        incvec <- df1$lakenum == i
#        dtp.pred <- df1$tp.sc[incvec] - tppred[incvec] +
#            exp(df1$d1)[incvec]
        dtp.pred <- df1$dtp[incvec]
        chl.loc <- df1$chl[incvec]
        chl.loc.sc <- (chl.loc - min(chl.loc))/diff(range(chl.loc))*
            diff(range(dtp.pred)) + min(dtp.pred)

        dftemp <- df1[incvec,]
        dftemp <- dftemp[order(dftemp$yday),]
        dt <- dftemp$yday[-1] - dftemp$yday[-nrow(dftemp)]
        t.mn <- 0.5*(dftemp$yday[-1] + dftemp$yday[-nrow(dftemp)])
        dp <- (dftemp$dtp[-1] - dftemp$dtp[-nrow(dftemp)])/dt
        p.mn <- 0.5*(dftemp$dtp[-1] + dftemp$dtp[-nrow(dftemp)])
        dchl <- (dftemp$chl[-1] - dftemp$chl[-nrow(dftemp)])/dt
        chl.mn <- 0.5*(dftemp$chl[-1] - dftemp$chl[-nrow(dftemp)])
        selvec <- dchl > 0

        d3 <- dftemp$d3[1]
        ret <- 1/dftemp$flush.rate[1]*365
        coef <- exp(mn.val["tp"])*exp(d3)/(exp(mn.val["chl"])^k[2])*(chl.mn[selvec])^(k[2]-1)

        pin <- rep(NA, times = length(t.mn))
        pin[selvec] <- ret*dp[selvec] + p.mn[selvec]+
            coef*(ret*dchl[selvec] + chl.mn[selvec])
        pin[!selvec] <- ret*dp[!selvec] + p.mn[!selvec]


        plot(t.mn, pin)
        points(t.mn[!selvec], pin[!selvec], pch = 16, col = "brown")
        abline(h=0)
        meanload[i] <- mean(pin, na.rm = T)

    }

    dftemp <- unique.data.frame(df1[, c("lakenum", "lake", "logit_crop")])
    dftemp <- dftemp[order(dftemp$lakenum),]
    print(dftemp)
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(2,2))
    plot(plogis(dftemp$logit_crop), meanload)
    incvec <- meanload > 500 & plogis(dftemp$logit_crop) < 0.2
    points(plogis(dftemp$logit_crop)[incvec], meanload[incvec], pch = 16)
    print(dftemp[incvec,])
    stop()

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0))
    plot(df1$yday, log(tsspred)-log(df1$tss.sc), pch = 21, col = "grey",
         bg = "white")
    abline(h=0)
    require(mgcv)
    resp <- log(tsspred)-log(df1$tss.sc)
    mod <- gam(resp ~ s(df1$yday,k = 6))
    predout <- predict(mod, se.fit = T)
    mn <- predout$fit
    up <- predout$fit + 2*predout$se.fit
    dn <- predout$fit - 2*predout$se.fit
    iord <- order(df1$yday)
    lines(df1$yday[iord], mn[iord])
    lines(df1$yday[iord], up[iord], lty = "dashed")
    lines(df1$yday[iord], dn[iord], lty = "dashed")

    plot(df1$yday, log(tppred)-log(df1$tp.sc), pch = 21, col = "grey",
         bg = "white")
    abline(h=0)
    resp <- log(tppred)-log(df1$tp.sc)
    mod <- gam(resp ~ s(df1$yday, k = 6))
    predout <- predict(mod, se.fit = T)
    mn <- predout$fit
    up <- predout$fit + 2*predout$se.fit
    dn <- predout$fit - 2*predout$se.fit
    iord <- order(df1$yday)
    lines(df1$yday[iord], mn[iord])
    lines(df1$yday[iord], up[iord], lty = "dashed")
    lines(df1$yday[iord], dn[iord], lty = "dashed")

    rms <- function(x,y) {
        val <- sqrt((sum((x-y)^2))/length(x))
        return(val)
    }
    print(rms(log(tsspred), log(df1$tss.sc)))
    print(rms(log(tppred), log(df1$tp.sc)))

    stop()

    boxplot(split(log(df1$nvss.sc), df1$lakenum))
    stop()
    a <- mean(varout$a)
    k <- apply(varout$k, 2, mean)
    u <- apply(varout$u, 2, mean)
    plot(log(exp(u) + exp(a)*df1$chl.sc^k[3]), log(df1$nvss.sc))
    abline(0,1)
    stop()
    mod <- lm(log(nvss.sc) ~ u)
    print(summary(mod))
    abline(mod, lty = "dashed")
    stop()
#    incvec <- df1$lakenum == 8
#    plot(log(df1$chl.sc)[incvec], log(df1$tss.sc)[incvec])
#    stop()

    df1$u <- apply(varout$u, 2, mean)
    b <- apply(varout$b, 2, mean)

    dftemp <- data.frame(seasnum = 1:nperiod, b = b)
    k <- apply(varout$k, 2, mean)

    print(nrow(df1))
    df1 <- merge(df1, dftemp, by = c( "seasnum"))
    print(nrow(df1))

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2))
    df1$pred.tss <- exp(df1$b)*df1$chl.sc^k[1] + exp(df1$u)
    plot(df1$yday,log(df1$tss.sc)- log(df1$pred.tss), col = "grey")
    err <- log(df1$tss.sc)- log(df1$pred.tss)
    abline(h = 0)
    require(mgcv)
    mod <- gam(err ~ s(yday), data = df1)
    iord <- order(df1$yday)
    predout <- predict(mod, se.fit = T)
    mn <- predout$fit
    up <- predout$fit + 2*predout$se.fit
    dn <- predout$fit - 2*predout$se.fit
    lines(df1$yday[iord], mn[iord])
    lines(df1$yday[iord], up[iord], lty = "dashed")
    lines(df1$yday[iord], dn[iord], lty= "dashed")

    d <- apply(varout$d, c(2,3), mean)
    dfd <- data.frame(lakenum = 1:nrow(d), d)
    names(dfd) <- c("lakenum", "d1", "d2old", "d3")
    df1 <- merge(df1, dfd, by = "lakenum")
    d2 <- apply(varout$d2, c(2,3), mean)
    print(d2)
    print(as.vector(d2))

    nlake <- max(df1$lakenum)
    seasnum <- rep(1:nperiod, times = nlake)
    print(seasnum)
    lakenum <- rep(1:nlake, times = rep(nperiod, times = nlake))
    print(lakenum)
    dftemp <- data.frame(seasnum = seasnum, lakenum = lakenum,
                        d2 = as.vector(d2))
    print(nrow(df1))
    df1 <- merge(df1, dftemp, by = c("seasnum", "lakenum"))
    print(nrow(df1))

    df1$pred <- exp(df1$d1) + exp(df1$d2)*exp(df1$u) +
        exp(df1$d3)*df1$chl.sc^k[2]

    dev.new()
    plot(log(df1$pred), log(df1$pred)-log(df1$tp.sc))
    abline(h=0)
    stop()

    plot(df1$yday, log(df1$tp.sc) - log(df1$pred), pch = 21,
         col = "grey", bg ="white")
    y <- log(df1$tp.sc) - log(df1$pred)
    require(mgcv)
    mod <- gam(y ~ s(yday), data = df1)
    iord <- order(df1$yday)
    predout <- predict(mod, se.fit = T)
    mn <- predout$fit
    up <- predout$fit + 2*predout$se.fit
    dn <- predout$fit - 2*predout$se.fit
    lines(df1$yday[iord], mn[iord])
    lines(df1$yday[iord], up[iord], lty = "dashed")
    lines(df1$yday[iord], dn[iord], lty = "dashed")
    abline(h = 0)
    stop()

    d <- apply(varout$d, 2, mean)
    dftemp <- data.frame(lakenum = 1:length(d), d = d)
    print(nrow(df1))
    df1 <- merge(df1, dftemp, by = "lakenum")
    print(nrow(df1))
    stop()




    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(3,3))
    lake.u <- sort(unique(df1$lakenum))
    for (i in lake.u[3]) {

        incvec <- i == df1$lakenum
        plot(log(df1$chl.sc)[incvec], log(df1$tss.sc)[incvec])
        abline(b, k)
        delt <- log(df1$tss.sc)[incvec] - (b + k*log(df1$chl.sc)[incvec])
        hist(delt)
        abline(v = exp(muu[2]))

        plot(df1$yday[incvec], log(df1$pred.tss)[incvec]-log(df1$tss.sc)[incvec]
             )
        abline(h=0)
        x <- df1$yday[incvec]
        y <- log(df1$pred.tss)[incvec] - log(df1$tss.sc)[incvec]
        mod <- gam(y ~ s(x, k = 6))
        iord <- order(x)
        lines(x[iord], predict(mod)[iord])
    }


    stop()
    k <- apply(varout$k, 2, mean)
    df1$u <- apply(varout$u, 2, mean)
}
fitout <- tss.explore(mo.all, runmod = T)
#tss.explore(moi3.all, varout.tp3, runmod = F)
## varout.tp.1 : b: time, all d: lake
## varout.tp.2 : b: time, d3 time
## varout.tp.3 : b: time, d3 time and lake.
