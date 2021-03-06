# 1/23/2019: model TSS rather than nvss and vssn
# 1/31/2019: model TSS with changing coef with chl
## 12.18.2020: Model seston
## add diagnostics
tss.explore <- function(df1, matout = NULL,varout = NULL, varout.n = NULL,
                        runmod = T, xvalid = F) {

    df1$vss <- as.numeric(as.character(df1$vss))
    df1$tss <- as.numeric(as.character(df1$tss))
    df1$dtp <- as.numeric(as.character(df1$dtp))
    df1$nvss <- as.numeric(as.character(df1$nvss))
    df1$tp <- as.numeric(as.character(df1$tp))
    df1$chl <- as.numeric(as.character(df1$chl))
    df1$ntu <- as.numeric(as.character(df1$ntu))
    df1$srp <- as.numeric(as.character(df1$srp))
    df1$don <- as.numeric(as.character(df1$don))

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
        #df1$yday.q <- ceiling(df1$month*0.5)
        df1$yday.q <- df1$month
    }

    print(table(df1$month))
    df1$yday.q <- factor(df1$yday.q)
    print(table(df1$yday.q))

    incvec <- ! is.na(df1$vss) & ! is.na(df1$chl) & ! is.na(df1$tp) &
        ! is.na(df1$nvss)
    df1 <- df1[incvec,]
    df1$lake <- factor(df1$lake)

    df1$lake <- factor(df1$lake)
    df1$lakenum <- as.numeric(df1$lake)
    df1$seasnum <- as.numeric(df1$yday.q)

    ## drop one big outlier
    incvec <- log(df1$tp - df1$dtp) < 0
    df1 <- df1[!incvec, ]

    print(summary(df1$chl))
    print(summary(df1$tp))
    print(summary(df1$tn))
    print(nrow(df1))

    ## simple summary N:P
    print("TN:TP")
    print(summary(df1$tn*1000/df1$tp))

    varlist<- c("tss", "chl", "tp", "ntu", "tn")
    mn.val <- apply(df1[, varlist],2,function(x) exp(mean(log(x))))
    print(mn.val)
    save(mn.val, file = "mn.val.mo.rda")

    for (i in varlist) df1[,i] <- df1[,i]/mn.val[i]
    df1$dtp <- df1$dtp/mn.val["tp"]
    df1$srp <- df1$srp/mn.val["tp"]
    df1$vss <- df1$vss/mn.val["tss"]
    df1$nvss <- df1$nvss/mn.val["tss"]
    df1$dtn <- df1$dtn/mn.val["tn"]
    df1$don <- df1$don/mn.val["tn"]

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
            vector[nseas] etab;
            real<lower = 0> sigb;

            real k[3];

            vector[2] mud;
            real<lower = 0> sigd[2];
            vector[nseas] etad1;
            vector[nseas] etad2;

            real<lower = 0> sigtss;
            real<lower = 0> sigtp;

        }
        transformed parameters {
            vector[n] u;
            vector[n] tp_mn;
            vector[n] tss_mn;
            vector[nseas] d1;
            vector[nseas] d2;
            vector[nseas] b;

            b = mub + etab*sigb;
            u = muu + etau*sigu;

            d1 = mud[1] + etad1*sigd[1];
            d2 = mud[2] + etad2*sigd[2];

            for (i in 1:n) {
               tss_mn[i] = log_sum_exp(b[seasnum[i]] + k[1]*chl[i], u[i]);

               tp_mn[i] = log_sum_exp(mud[1]+k[2]*chl[i],
                                      d2[seasnum[i]]+k[3]*u[i]);

            }
        }
        model {

            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);

            mub ~ normal(0,3);
            etab ~ normal(0,1);
            sigb ~ cauchy(0,3);

            mud ~ normal(0,3);
            sigd ~ cauchy(0,3);
            etad1 ~ normal(0,1);
            etad2 ~ normal(0,1);

       //     k[1] ~ normal(0.85,0.01);  // from VSS model
           k[1] ~ normal(1,1);
            k[2] ~ normal(1,1);
            k[3] ~ normal(1,1);

            sigtp ~ cauchy(0,3);
            sigtss ~ normal(0.1, 0.02);

            tss ~ student_t(4,tss_mn, sigtss);
            tp ~ student_t(4,tp_mn, sigtp);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2, na.rm = T)/
                                 min(sum(! is.na(x)), sum(! is.na(y))))
    gettp <- function(df, varout.loc, withu) {

        k <- apply(varout.loc$k, 2, mean)
        mud <- apply(varout.loc$mud, 2, mean)
        d1 <- apply(varout.loc$d1, 2, mean)
        d2 <- apply(varout.loc$d2, 2, mean)
        mub <- mean(varout.loc$mub)
        b <- apply(varout.loc$b, 2, mean)

        flag <- rep(F, times = nrow(df))

        if (! withu) {
            u <- df$tss - exp(b[df$seasnum])*df$chl^k[1]
            incvec <- u < 0
            u[incvec] <- 0
            flag[incvec] <- T
           }
        else {
            u <- exp(apply(varout.loc$u, 2, mean))
        }

        tppred <- exp(mud[1])*df$chl^k[2] + df$dtp +
            exp(d2[df$lakenum])*u^k[3]

        return(list(tppred = tppred, flag=flag))
    }
    extractvars <- c("k", "u", "d1", "d2", "sigd", "sigtp", "sigtss","mud",
                     "mub", "sigb", "b")

    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 3
        options(mc.cores = nchains)

        if (! xvalid) {
            datstan <- list(n = nrow(df1),
                            nlake = max(df1$lakenum),lakenum = df1$lakenum,
                            nseas = max(df1$seasnum), seasnum = df1$seasnum,
                            nvss = df1$nvss,
                            tp = log(df1$tp - df1$dtp),
                            dtp = df1$dtp,
                            vss = df1$vss,
                            chl = log(df1$chl),
                            tss = log(df1$tss))
            print(str(datstan))

            fit <- stan(model_code = modstan,
                        data = datstan, iter = 1800, chains = nchains,
                        warmup = 600, thin= 2,
                        control = list(adapt_delta = 0.98, max_treedepth = 14))


            save(fit, file = "fitmo.rda")
            varout <- extract(fit, pars = extractvars)

            tp.pred <- gettp(df =df1, varout.loc = varout, withu = T)[[1]]
            dev.new()
            plot(log(tp.pred), log(df1$tp))
            abline(0,1)
            print(rmsout(log(tp.pred), log(df1$tp)))
            return(varout)
        }
        else {

            set.seed(31) # original value is 3
            require(loo)
            nfold <- 5
            ik <- kfold_split_random(nfold, nrow(df1))

            for (jj in 1:nfold) {

                dftemp2 <- df1[ik == jj,]
                dftemp <- df1[ik != jj,]
                datstan <- list(n = nrow(dftemp),
                          nlake = max(dftemp$lakenum),lakenum = dftemp$lakenum,
                           nseas = max(dftemp$seasnum), seasnum = dftemp$seasnum,
                                nvss = dftemp$nvss,
                                tp = log(dftemp$tp),
                                dtp = dftemp$dtp,
                                vss = dftemp$vss,
                                chl = log(dftemp$chl),
                                tss = log(dftemp$tss))
                print(str(datstan))

                fit <- stan(model_code = modstan,
                            data = datstan, iter = 1000, chains = nchains,
                            warmup = 500, thin= 2,
                            control = list(adapt_delta = 0.98, max_treedepth = 14))

                varout <- extract(fit, pars = extractvars)
                tp.pred <- gettp(dftemp2, varout, withu = FALSE)[[1]]

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


    tp.pred <- gettp(df1, varout, withu = T)[[1]]
    cat("Internal validation:", rmsout(log(tp.pred), log(df1$tp)), "\n")

    k <- apply(varout$k, 2, mean)
    b <- apply(varout$b, 2, mean)
    mub <- mean(varout$mub)
    u <- apply(varout$u, 2, mean)
    mud <- apply(varout$mud, 2, mean)
    d2 <- apply(varout$d2, 2, mean)

    psed <- u*k[3] + d2[df1$lakenum]
    prop.sed <- exp(psed - log(df1$tp))
    print(summary(prop.sed))
    stop()

    ## plot dtp vs. chl, dtn vs. chl
    png(width = 6, height = 2.5, pointsize = 6, units = "in",
        res = 600, file = "dissolved.png")
    par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0))
    incvec <- log(df1$dtp*mn.val["tp"]) > 0.5 # one outlier
    plot(log(df1$chl*mn.val["chl"])[incvec], log(df1$dtp*mn.val["tp"])[incvec],
         pch = 21, col = "grey39", bg = "white",
         ylab = expression(P[diss]~(mu*g/L)),
         xlab = expression(Chl~(mu*g/L)), axes = F)
    logtick.exp(0.001, 10, c(1,2), c(F,T))
    abline(log(1.16 - 0.84), 1)
    y <- log(df1$don*mn.val["tn"]*1000)
    incvec <- y > log(0.1)
    plot(log(df1$chl*mn.val["chl"])[incvec], y[incvec],
         pch = 21, col = "grey39", bg = "white",
         ylab = expression(DON~(mu*g/L)),
         xlab = expression(Chl~(mu*g/L)), axes = F)
    logtick.exp(0.001, 10, c(1,2), c(F,T))
    abline(log(9.22 - 6.39), 1)
    dev.off()

    stop()


    prop.sed.low <- tapply(prop.sed, df1$lakenum, quantile, prob = 0.05)
    prop.sed.hi <- tapply(prop.sed, df1$lakenum, quantile, prob = 0.95)
    prop.sed.med <- tapply(prop.sed, df1$lakenum, quantile, prob = 0.5)

    dftemp <- unique(data.frame(df1[, c("lakenum", "flush.rate")]))
    dftemp <- dftemp[order(dftemp$flush.rate),]
    plot(dftemp$flush.rate, prop.sed.med,
         ylim = range(c(prop.sed.low, prop.sed.hi)))
    segments(dftemp$flush.rate, prop.sed.low,
             dftemp$flush.rate, prop.sed.hi)

    stop()

    credint <- c(0.025, 0.5, 0.975)
    b <- apply(varout$b, 2, mean)
    blim <- exp(apply(varout$b, 2, quantile, prob = credint))
    chl0 <- 10
    chl.ss <- 1/blim/mn.val["tss"]*chl0^(1-k[1])*mn.val["chl"]^k[1]

    png(width = 3, height = 2.5, pointsize = 6, units = "in",
        res = 600, file = "chl.ss.png")
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0))
    plot(1:12, chl.ss[2,], ylim= range(chl.ss), type = "n",
         ylab = expression(Chl/SS~(mu*g/m*g)),
         xlab = "Month", axes = F)
    axis(2)
    axis(1)
    box(bty = "l")
    segments(1:12, chl.ss[1,], 1:12, chl.ss[3,], col = "grey39")
    points(1:12, chl.ss[2,], pch = 21, col = "grey39", bg = "white")
    dev.off()

    k <- apply(varout$k, 2, mean)
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(3,4))
    for (i in 1:12) {
        incvec <- df1$seasnum == i
        plot(log(df1$chl), log(df1$tss),type = "n", axes = F,
             xlab = expression(Chl~(mu*g/L)),
             ylab = expression(TSS~(mu*g/L)))
        points(log(df1$chl)[incvec],
               log(df1$tss)[incvec],
               pch = 21, col = "grey39", bg = "white")
        logtick.exp(0.001, 10, c(1,2), c(F,F))
        abline(b[i], k[1])
    }

    print("*** coefficients for chl relationship ***")
    print("*** TP ***")
    print(quantile(exp(varout$mud[,1] - varout$k[,2]*log(mn.val["chl"]) + log(mn.val["tp"])), prob = credint))
    print(quantile(varout$k[,2], prob =credint))

    print("*** TN ***")
    print(quantile(exp(varout.n$mud[,1] - varout.n$k[,1]*log(mn.val["chl"]) + log(mn.val["tn"]) + log(1000)), prob = credint))
    print(quantile(varout.n$k[,1], prob =credint))

    print("*** suspended sediment ***")
    print(quantile(exp(varout$mud[,2] - varout$k[,3]*log(mn.val["tss"]) +
                           log(mn.val["tp"])), prob =credint))
    d2 <- apply(varout$d2, 2, mean)
    print(range(exp(d2 - varout$k[,3]*log(mn.val["tss"]) +
                           log(mn.val["tp"]))))
    print(quantile(varout$k[,3], prob = credint))

    print("*** volatile suspended sediment N ***")
    print(quantile(exp(varout.n$mud[,2] - varout.n$k[,2]*log(mn.val["tss"]) +
                           log(mn.val["tn"]*1000)), prob = credint))
    d2n <- apply(varout.n$d2, 2, mean)
    print(range(exp(d2n - varout.n$k[,2]*log(mn.val["tss"]) +
                        log(mn.val["tn"]*1000))))

    print(quantile(varout.n$k[,2], prob = credint))

    grey.t1 <- adjustcolor("grey20", alpha = 0.5)
    grey.t2 <- adjustcolor("grey70", alpha = 0.5)

    dftemp <- unique.data.frame(df1[, c("lakenum", "flush.rate", "logit_crop")])
    dftemp <- dftemp[order(dftemp$lakenum),]

    um <- mean(apply(varout$u, 2, mean))
    umn <- mean(apply(varout.n$u,2, mean))
    print("*** mean suspended seds ***")
    print(exp(um + log(mn.val["tss"])))
    print(exp(umn + log(mn.val["tss"])))

    predout <- matrix(NA, ncol = 3, nrow = 15)
    for (i in 1:15) {
        y <- mn.val["tp"]*exp(varout$d2[,i])*exp(um)^varout$k[,3]/
            (exp(um)*mn.val["tss"])
        predout[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))

    }
    png(width = 6, height = 2.5, pointsize = 8, units = "in", res = 600,
        file = "lakessp.png")
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l", mfrow = c(1,2))
    plot(log(dftemp$flush.rate), predout[,2], type = "n",axes = F,
         xlab = "Flush rate (1/yr)",
         ylab = expression(P/SS[np]~(mu*g/m*g)), ylim = range(predout))
    logtick.exp(0.001, 10, c(1), c(F,F))
    segments(log(dftemp$flush.rate), predout[,1],
             log(dftemp$flush.rate), predout[,3], col = "grey39")
    points(log(dftemp$flush.rate), predout[,2], col = "grey39", pch = 21,
           bg = "white")
    axis(2)
    box(bty = "l")

    predout2 <- matrix(NA, ncol = 3, nrow = 12)
    for (i in 1:12) {
        y2 <- 1000*mn.val["tn"]*exp(varout.n$d2[,i])*exp(umn)^varout.n$k[,2]/
            (exp(umn)*mn.val["tss"])
        predout2[i,] <- quantile(y2, prob = c(0.025, 0.5, 0.975))
    }

    plot(1:12, predout2[,2], type = "n",axes = F,
         xlab = "Month",
         ylab = expression(N/VSS[np]~(mu*g/m*g)), ylim = range(predout2))
    axis(2)
    axis(1)
    segments(1:12, predout2[,1],
             1:12, predout2[,3], col = "grey39")
    points(1:12, predout2[,2], col = "grey39", pch = 21,
           bg = "white")
    axis(2)
    box(bty = "l")
    dev.off()

    ## plot effect of Chl on N:P
    ## start N:P at 4 mug/L because of limits to TP model
    chl0 <- seq(log(4/mn.val["chl"]), log(max(df1$chl)), length = 40)
    predout <- matrix(NA, ncol = 3, nrow = length(chl0))
    grey.t <- adjustcolor("grey39", alpha = 0.5)
    for (i in 1:length(chl0)) {
        y <- mn.val["tp"]*exp(varout$mud[,1])*exp(chl0[i])^varout$k[,2]
        y2 <- 1000*mn.val["tn"]*exp(varout.n$mud[,1])*exp(chl0[i])^varout.n$k[,1]
        predout[i,] <- quantile(y2/y, prob = c(0.05, 0.5, 0.95))
    }
    print(min(predout[,3]))
    print(max(predout[,1]))

    png(width = 3, height = 2, pointsize = 6, units ="in", res = 600,
        file = "np.png")
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0))
    xraw <- chl0 + log(mn.val["chl"])
    plot(xraw, predout[,2],
         ylim = range(predout), type = "n", xlab = expression(Chl~(mu*g/L)),
         ylab = "N:P", axes= F)
    polygon(c(xraw, rev(xraw)), c(predout[,1], rev(predout[,3])),
            col = grey.t, border = NA)
    lines(xraw, predout[,2])
    axis(2)
    logtick.exp(0.001, 10, c(1), c(F,F))
    dev.off()

    xnew <- seq(log(0.1), log(300), length = 50)
    predouta <- matrix(NA, ncol = 3, nrow = length(xnew))
    predout2a <- matrix(NA, ncol = 3, nrow = length(xnew))
    nsamp <- nrow(varout$mud)
    nsamp2 <- nrow(varout.n$mud)
    for (i in 1:length(xnew)) {
        y <- varout$mud[,1]+
            varout$k[,2]*(xnew[i]-log(mn.val["chl"])) + log(mn.val["tp"])
        predouta[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))
        y2 <- varout.n$mud[,1] +
            varout.n$k[,1]*(xnew[i] - log(mn.val["chl"])) +
                log(mn.val["tn"]) + log(1000)
        predout2a[i,] <- quantile(y2, prob = c(0.025, 0.5, 0.975))
    }

    mub <- mean(varout$mub)
    k <- apply(varout$k, 2, mean)

    png(width = 6, height = 2.5, pointsize = 6, units = "in",
        res = 600, file = "molimits.png")
    par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0))
    plot(log(df1$chl) + log(mn.val["chl"]),
         log(df1$tp - df1$dtp) + log(mn.val["tp"]),
         pch = 21, col = "grey39", bg = "white", axes = F,
         xlab = expression(Chl~(mu*g/L)),
         ylab = expression(P[part]~(mu*g/L)))
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    lines(xnew, predouta[,2])
    polygon(c(xnew, rev(xnew)), c(predouta[,1], rev(predouta[,3])),
            col = grey.t2, border = NA)

    plot(log(df1$chl) + log(mn.val["chl"]),
         log(df1$tn - df1$dtn) + log(mn.val["tn"]) + log(1000),
         pch = 21, col = "grey39", bg = "white", axes = F,
         xlab = expression(Chl~(mu*g/L)),
         ylab = expression(N[part]~(mu*g/L)))
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    lines(xnew, predout2a[,2])
    polygon(c(xnew, rev(xnew)), c(predout2a[,1], rev(predout2a[,3])),
            col = grey.t2,  border = NA)
    dev.off()

    stop()
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


    return()

}

#varout.mo.d10.d2T <- tss.explore(moi3.all, runmod = T, xvalid= F)
#matout.mo.d1T.d2L.4 <-  tss.explore(moi3.all, runmod = T, xvalid= T)

tss.explore(moi3.all, matout = matout.mo.d10.d2L, varout = varout.mo.d10.d2L,
            varout.n = varout.mon.d10.d2T,runmod = F, xvalid = F)
