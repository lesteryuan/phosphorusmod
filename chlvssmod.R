# 1/23/2019: model TSS rather than nvss and vssn
# 1/31/2019: model TSS with changing coef with chl
## 12.18.2020: Model seston
## add diagnostics
tss.explore <- function(df1, matout = NULL,varout = NULL,
                        runmod = T, xvalid = F) {

    df1$vss <- as.numeric(as.character(df1$vss))
    df1$ntu <- as.numeric(as.character(df1$ntu))
    df1$tss <- as.numeric(as.character(df1$tss))

    df1$dtp <- as.numeric(as.character(df1$dtp))
    df1$nvss <- as.numeric(as.character(df1$nvss))
    df1$tp <- as.numeric(as.character(df1$tp))
    df1$chl <- as.numeric(as.character(df1$chl))
    df1$srp <- as.numeric(as.character(df1$srp))

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

    ## model for dtp to chl relationship
    dev.new()
    plot(log(df1$chl), log(df1$dtp - df1$srp), axes= F,xlab="Chl",
         ylab = "DTP - SRP")
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    mod <- lm(log(df1$dtp - df1$srp) ~ log(df1$chl))
    abline(mod)
    print(summary(mod))

    dftemp <- unique.data.frame(df1[, c("lake", "lakenum", "logit_crop",
                                        "flush.rate")])
    dftemp <- dftemp[order(dftemp$lakenum),]
    dftemp$d <- apply(varout.vss$d, 2, mean)
    dev.new()
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0))
    plot(dftemp$flush.rate, exp(dftemp$d), xlab = "Flush rate",
         ylab = "VSS P-content", bty = "l")


    incvec <- dat.merge.all.small$us.l3code == 40
    adj <- 3.375*dat.merge.all.small$chl^0.373
    points(log(dat.merge.all.small$chl)[incvec],
           log(dat.merge.all.small$ptl.result - adj)[incvec],
           pch = 16, col = "red")


    varlist<- c("tss", "chl", "tp", "ntu")
    mn.val <- apply(df1[, varlist],2,function(x) exp(mean(log(x))))
    print(mn.val)
    save(mn.val, file = "mn.val.mo.rda")

    for (i in varlist) df1[,i] <- df1[,i]/mn.val[i]
    df1$dtp <- df1$dtp/mn.val["tp"]
    df1$srp <- df1$srp/mn.val["tp"]
    df1$vss <- df1$vss/mn.val["tss"]
    df1$nvss <- df1$nvss/mn.val["tss"]

    ## drop one big outlier
    incvec <- log(df1$chl) < -2 & log(df1$tp - df1$dtp) < -3
#    points(log(df1$chl)[incvec], log(df1$tp - df1$dtp)[incvec], pch = 16)
    df1 <- df1[!incvec, ]

    k <- apply(varout.vss$k, 2, mean)
    mud <- apply(varout.vss$mud, 2, mean)
    df1$u <- apply(varout.vss$u, 2, mean)
    vss.p <- exp(mud[2])*exp(df1$u)^k[3]
    nvss.p <- exp(mud[1])*df1$nvss^k[2]
    dev.new()
    plot(log(df1$chl), log(df1$tss))

    stop()

    dev.new()
    plot(log(df1$chl), log(df1$tp - df1$dtp), col = "grey",
         ylim = range(c(log(df1$tp - df1$dtp - vss.p), log(df1$tp-df1$dtp)), na.rm = T))
    incvec <- df1$flush.rate < 0.3
    points(log(df1$chl), log(df1$tp - df1$dtp - vss.p), pch = 16,
           cex = 0.6)
    abline(-1.644, 0.95)

    stop()

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

            real k[4];

            vector[3] mud;
//            real<lower = 0> sigd[2];
//            vector[nlake] etad1;
//            vector[nlake] etad2;
//            vector[nlake] etad3;

//            vector<lower = 0>[3] sigd;
//            matrix[nlake,3] etad;

            real<lower = 0> sigvss;
            real<lower = 0> sigtp;

        }
        transformed parameters {
            vector[n] u;
            vector[n] tp_mn;
            vector[n] vss_mn;
//            vector[nlake] d1;
 //           vector[nlake] d2;
 //           vector[nlake] d3;
//            matrix[nlake, 3] d;

            u = muu + etau*sigu;

//            d1 = mud[1] + etad1*sigd[1];
//            d2 = mud[2] + etad2*sigd[2];
//            d3 = mud[3] + etad3*sigd[3];
//            for (i in 1:3) d[,i] = mud[i] + etad[,i]*sigd[i];

            for (i in 1:n) {
               vss_mn[i] = exp(mub)*chl[i]^k[1] + exp(u[i]);

               tp_mn[i] = exp(mud[1])*nvss[i]^k[2] +
                          exp(mud[2])*exp(u[i])^k[3] +
                          exp(mud[3])*chl[i]^k[4] + dtp[i];
            }
        }
        model {
            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);

            mub ~ normal(0,3);

            mud ~ normal(0,3);
//            sigd ~ cauchy(0,3);
//            etad1 ~ normal(0,1);
//            etad2 ~ normal(0,1);
//            etad3 ~ normal(0,1);
//            for (i in 1:3) etad[,i] ~ normal(0,1);

            k[1] ~ normal(0.807,0.012);
            k[2] ~ normal(1,1);
            k[3] ~ normal(1,1);
            k[4] ~ normal(1,1);
//            k[3] ~ normal(0.877,0.033);

            sigtp ~ cauchy(0,3);
            sigvss ~ cauchy(0,3);

            vss ~ lognormal(log(vss_mn), sigvss);
            tp ~ lognormal(log(tp_mn), sigtp);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2)/length(x))
    gettp <- function(df, varout, withu) {
        k <- apply(varout$k, 2, mean)
        mud <- apply(varout$mud, 2, mean)
 #       d1 <- apply(varout$d1, 2, mean)
 #       d2 <- apply(varout$d2, 2, mean)
#        d3 <- apply(varout$d3, 2, mean)
        mub <- mean(varout$mub)

        if (! withu) {
            u <- df$vss - exp(mub)*df$chl^k[1]
            incvec <- u < 0
            u[incvec] <- 0
           }
        else {
            u <- exp(apply(varout$u, 2, mean))
        }

        dev.new()
        plot(log(df$chl), df$nvss/(df$nvss+u))

        tppred <- exp(mud[1])*df$nvss^k[2] +
                  exp(mud[2])*u^k[3] +
                  exp(mud[3])*df$chl^k[4] + df$dtp

        return(tppred)
    }
    extractvars <- c("mud", "k", "mub", "u")

    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 3
        options(mc.cores = nchains)

        if (! xvalid) {
            datstan <- list(n = nrow(df1),
                            nlake = max(df1$lakenum),lakenum = df1$lakenum,
                            nseas = nperiod, seasnum = df1$seasnum,
                            nvss = df1$nvss,
                            tp = df1$tp,
                            dtp = df1$dtp,
                            vss = df1$vss,
                            chl = df1$chl)
            print(str(datstan))

            fit <- stan(model_code = modstan,
                        data = datstan, iter = 1000, chains = nchains,
                        warmup = 500, thin= 2,
                        control = list(adapt_delta = 0.98, max_treedepth = 14))
            varout <- extract(fit, pars = extractvars)
            tp.pred <- gettp(df1, varout, withu = T)
            dev.new()
            plot(log(tp.pred), log(df1$tp))
            abline(0,1)
            print(rmsout(log(tp.pred), log(df1$tp)))

            return(varout)
        }
        else {

            set.seed(1)
            require(loo)
            nfold <- 5
            ik <- kfold_split_random(nfold, nrow(df1))

            for (jj in 1:nfold) {

                dftemp2 <- df1[ik == jj,]
                dftemp <- df1[ik != jj,]
                datstan <- list(n = nrow(dftemp),
                                nlake = max(dftemp$lakenum),lakenum = dftemp$lakenum,
                                nseas = nperiod, seasnum = dftemp$seasnum,
                                nvss = dftemp$nvss,
                                tp = dftemp$tp,
                                dtp = dftemp$dtp,
                                vss = dftemp$vss,
                                chl = dftemp$chl)
                print(str(datstan))

                fit <- stan(model_code = modstan,
                            data = datstan, iter = 1000, chains = nchains,
                            warmup = 500, thin= 2,
                            control = list(adapt_delta = 0.98, max_treedepth = 14))

                varout <- extract(fit, pars = extractvars)
                tp.pred <- gettp(dftemp2, varout, withu = FALSE)

                matval <- data.frame(pred = log(tp.pred), obs = log(dftemp2$tp))

                if (jj == 1) {
                    matall <- matval

                }
                else {
                    matall <- rbind(matall, matval)
                }
            }

            return(matall)
        }
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
#varout.00 <- tss.explore(mo.all, runmod = T, xvalid= F)
