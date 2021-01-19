# 1/23/2019: model TSS rather than nvss and vssn
# 1/31/2019: model TSS with changing coef with chl
## 12.18.2020: Model seston
## add diagnostics
tss.explore <- function(df1, matout = NULL,varout = NULL,
                        runmod = T, xvalid = F) {

    df1$vss <- as.numeric(as.character(df1$vss))
    df1$tss <- as.numeric(as.character(df1$tss))
    df1$dtp <- as.numeric(as.character(df1$dtp))
    df1$nvss <- as.numeric(as.character(df1$nvss))
    df1$tp <- as.numeric(as.character(df1$tp))
    df1$chl <- as.numeric(as.character(df1$chl))
    df1$ntu <- as.numeric(as.character(df1$ntu))
    df1$srp <- as.numeric(as.character(df1$srp))

    date0 <- strptime(paste(df1$month, df1$day, "2004", sep = "-"),
                      format = "%m-%d-%Y")
    df1$yday <- date0$yday

    print(nrow(df1))
    df1 <- merge(df1, res.dat[, c("MU..", "flush.rate", "mean.depth")], by.x = "lake",
                 by.y = "MU..")
    print(nrow(df1))
    print(summary(df1$flush.rate))

    quart <- FALSE
    if (quart) {
        nperiod <- 4
        df1$yday.q <- 1
        incvec <- df1$month >= 3 & df1$month <= 5
        df1$yday.q[incvec]<- 2
        incvec <- df1$month >= 6 & df1$month <= 8
        df1$yday.q[incvec]<- 3
        incvec <- df1$month >= 9 & df1$month <= 11
        df1$yday.q[incvec]<- 4
    }
    else {
        df1$yday.q <- ceiling(df1$month*0.5)
    }

    print(table(df1$month))
    df1$yday.q <- factor(df1$yday.q)
    print(table(df1$yday.q))

    incvec <- ! is.na(df1$vss) & ! is.na(df1$chl) & ! is.na(df1$tp) &
        ! is.na(df1$nvss)
    df1 <- df1[incvec,]
    df1$lake <- factor(df1$lake)

    df1$lake <- factor(df1$lake)
    print(table(df1$lake))
    df1$lakenum <- as.numeric(df1$lake)

    df1$seasnum <- as.numeric(df1$yday.q)

    ## model for dtp to chl relationship
    incvec <- df1$dtp > df1$srp
    dftemp <- df1[incvec,]
    pmean <- tapply(log(dftemp$dtp), dftemp$lake, mean)
    cmean <- tapply(log(dftemp$chl), dftemp$lake, mean)
    dmean <- tapply(log(dftemp$mean.depth), dftemp$lake, mean)
    mod1 <- lm(pmean ~ cmean)
    mod2 <- lm(pmean ~ dmean)
    print(summary(mod2))
#    plot(dmean, pmean)
#    abline(mod2)

    predout <- predict(mod2, se.fit = T)
    print(predout)


    load("cutp.depth.rda")
#    abline(v = cutp.depth)
    cutm <- 0.5*(cutp.depth[-1] + cutp.depth[-length(cutp.depth)])
    cc <- coef(mod2)
    ## create shape of priors for dtp with informative
    ## priors in the range of MO I3 lakes
    prior1 <- cc[1] + cc[2]*cutm
    imin <- which(cutm < min(dmean))
#    prior1[imin] <- prior1[max(imin)+1]
    imax <- which(cutm > max(dmean))
#    prior1[imax] <- prior1[min(imax)-1]
#    plot(cutm, prior1)
#    points(dmean, pmean, pch = 16)
#    save(prior1, file = "prior1.rda")

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

       modstan <- '
        data {
            int n;
            int nlake;
            int lakenum[n];
            int nseas;
            int seasnum[n];
            vector[n] tp;
            vector[n] dtp;
            vector[n] chl;
            vector[n] tss;

        }
        parameters {
            real muu;
            real<lower = 0> sigu;
            vector[n] etau;

            real mub;

            real k[3];

            vector[2] mud;
            real<lower = 0> sigd[2];
            vector[nlake] etad1;
//            vector[nlake] etad2;
            vector[nseas] etad2a;
//            vector[nseas] etad2b;

            real<lower = 0> sigtss;
            real<lower = 0> sigtp;

        }
        transformed parameters {
            vector[n] u;
            vector[n] tp_mn;
            vector[n] tss_mn;
            vector[nlake] d1;
            vector[nseas] d2;

            u = muu + etau*sigu;

            d1 = mud[1] + etad1*sigd[1];
            d2 = mud[2] + etad2a*sigd[2];
//            for (i in 1:nseas) {
//                d2[,i] = mud[2] + etad2*sigd[2] + etad2a[i]*sigd[3];
//            }

            for (i in 1:n) {
               tss_mn[i] = exp(mub)*chl[i]^k[1] + exp(u[i]);

               tp_mn[i] = exp(d1[lakenum[i]])*exp(u[i])^k[2] +
                          exp(d2[seasnum[i]])*chl[i]^k[3] + dtp[i];
            }
        }
        model {
            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);

            mub ~ normal(0,3);

            mud ~ normal(0,3);
            sigd ~ cauchy(0,3);
            etad1 ~ normal(0,1);
//            etad2 ~ normal(0,1);
            etad2a ~ normal(0,1);
//            etad2b ~ normal(0,1);

            k[1] ~ normal(0.832,0.013);  // from VSS model
 //           k[1] ~ normal(1,1);
            k[2] ~ normal(1,1);
            k[3] ~ normal(1,1);

            sigtp ~ cauchy(0,3);
            sigtss ~ normal(0.1, 0.02);

            tss ~ student_t(4,log(tss_mn), sigtss);
            tp ~ student_t(4,log(tp_mn), sigtp);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2)/length(x))
    gettp <- function(df, varout.loc, withu) {

        k <- apply(varout.loc$k, 2, mean)
        mud <- apply(varout.loc$mud, 2, mean)
        d1 <- apply(varout.loc$d1, 2, mean)
        d2 <- apply(varout.loc$d2, 2, mean)
        mub <- mean(varout.loc$mub)

        if (! withu) {
            u <- df$tss - exp(mub)*df$chl^k[1]
            incvec <- u < 0
            u[incvec] <- 0
           }
        else {
            u <- exp(apply(varout.loc$u, 2, mean))
        }

        tppred <- rep(NA, times = nrow(df))
        for (i in 1:nrow(df)) {
            tppred[i] <- exp(d1[df$lakenum[i]])*u[i]^k[2] +
                exp(d2[df$seasnum[i]])*df$chl[i]^k[3] + df$dtp[i]
        }

        return(tppred)
    }
    extractvars <- c("mud", "k", "mub", "u", "d1", "d2", "sigd", "etad2a", "sigtp")

    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 6
        options(mc.cores = nchains)

        if (! xvalid) {
            datstan <- list(n = nrow(df1),
                            nlake = max(df1$lakenum),lakenum = df1$lakenum,
                            nseas = 6, seasnum = df1$seasnum,
                            nvss = df1$nvss,
                            tp = log(df1$tp),
                            dtp = df1$dtp,
                            vss = df1$vss,
                            chl = df1$chl,
                            tss = log(df1$tss))
            print(str(datstan))

            fit <- stan(model_code = modstan,
                        data = datstan, iter = 2000, chains = nchains,
                        warmup = 800, thin= 3,
                        control = list(adapt_delta = 0.98, max_treedepth = 14))

            save(fit, file = "fitmo.rda")
            varout <- extract(fit, pars = extractvars)

            tp.pred <- gettp(df =df1, varout.loc = varout, withu = T)
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
                           nseas = 4, seasnum = dftemp$seasnum,
                                nvss = dftemp$nvss,
                                tp = dftemp$tp,
                                dtp = dftemp$dtp,
                                vss = dftemp$vss,
                                chl = dftemp$chl,
                                tss = dftemp$tss)
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



    mub <- mean(varout$mub)
    u <- apply(varout$u, 2, mean)

    k <- apply(varout$k, 2, mean)
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(2,2))
    d1 <- apply(varout$d1, 2, quantile, prob = c(0.025, 0.5, 0.975))
    nd1 <- ncol(d1)
    plot(1:nd1, d1[2,], ylim = range(d1))
    segments(1:nd1, d1[1,], 1:nd1, d1[3,])
    stop()

    tsspred <- exp(u) + exp(mub)*df1$chl^k[1]

    plot(log(tsspred), log(df1$tss))
    abline(0,1)

    tppred <- gettp(df1, varout, withu = TRUE)
    print(summary(tppred))
    plot(log(tppred), log(df1$tp))
    abline(0,1)
    print(rmsout(log(tppred), log(df1$tp)))
    stop()
    print(summary(matout))

    plot(matout[,1], matout[,2])
    points(log(tppred), log(df1$tp), pch = 16, col = "red", cex = 0.6)
    abline(0,1)



    print(rmsout(matout[,1], matout[,2]))

    return()

}

varout.mo <- tss.explore(moi3.all, runmod = T, xvalid= F)
#mattss.d2l.d1tl <-  tss.explore(moi3.all, runmod = T, xvalid= T)

#tss.explore(moi3.all, matout = mattss.d2l.d1tl, varout = vartss.d2l.d1tl, runmod = F, xvalid = F)

tss.explore(moi3.all, varout = varout.mo, runmod = F, xvalid = F)
