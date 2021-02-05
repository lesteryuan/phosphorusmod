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
    df1$lakenum <- as.numeric(df1$lake)
    df1$seasnum <- as.numeric(df1$yday.q)

    ## drop one big outlier
    incvec <- log(df1$tp - df1$dtp) < 0
    df1 <- df1[!incvec, ]

    print(summary(df1$chl))
    print(summary(df1$tp))
    print(summary(df1$tn))
    print(nrow(df1))

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
            vector[nseas] etad1;
            vector[nlake] etad2;

            real<lower = 0> sigtss;
            real<lower = 0> sigtp;

        }
        transformed parameters {
            vector[n] u;
            vector[n] tp_mn;
            vector[n] tss_mn;
            vector[nseas] d1;
            vector[nlake] d2;

            u = muu + etau*sigu;

            d1 = mud[1] + etad1*sigd[1];
            d2 = mud[2] + etad2*sigd[2];
//            for (i in 1:nseas) {
//                d2[,i] = mud[2] + etad2*sigd[2] + etad2a[i]*sigd[3];
//            }

            for (i in 1:n) {
               tss_mn[i] = exp(mub)*chl[i]^k[1] + exp(u[i]);

               tp_mn[i] = exp(d1[seasnum[i]])*chl[i]^k[2] + dtp[i] +
                         exp(d2[lakenum[i]])*exp(u[i])^k[3];

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
            etad2 ~ normal(0,1);

            k[1] ~ normal(0.85,0.01);  // from VSS model
      //      k[1] ~ normal(1,1);
            k[2] ~ normal(1,1);
            k[3] ~ normal(1,1);

            sigtp ~ cauchy(0,3);
            sigtss ~ normal(0.1, 0.02);

            tss ~ student_t(4,log(tss_mn), sigtss);
            tp ~ student_t(4,log(tp_mn), sigtp);
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
            tppred[i] <- exp(d1[df$seasnum[i]])*df$chl[i]^k[2] + df$dtp[i] +
                         exp(d2[df$lakenum[i]])*u[i]^k[3]

        }

        return(tppred)
    }
    extractvars <- c("mud", "k", "mub", "u", "d1", "d2", "sigd", "sigtp", "sigtss")

    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 3
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
                        data = datstan, iter = 1800, chains = nchains,
                        warmup = 600, thin= 2,
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

            set.seed(3)
            require(loo)
            nfold <- 5
            ik <- kfold_split_random(nfold, nrow(df1))

            for (jj in 1:nfold) {

                dftemp2 <- df1[ik == jj,]
                dftemp <- df1[ik != jj,]
                datstan <- list(n = nrow(dftemp),
                          nlake = max(dftemp$lakenum),lakenum = dftemp$lakenum,
                           nseas = 6, seasnum = dftemp$seasnum,
                                nvss = dftemp$nvss,
                                tp = log(dftemp$tp),
                                dtp = dftemp$dtp,
                                vss = dftemp$vss,
                                chl = dftemp$chl,
                                tss = log(dftemp$tss))
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

    credint <- c(0.025, 0.5, 0.975)
    print(quantile(exp(varout$mud[,1] - varout$k[,2]*log(mn.val["chl"]) + log(mn.val["tp"])), prob = credint))
    print(quantile(varout$k[,2], prob =credint))

    print("*** suspended sediment ***")
    print(quantile(exp(varout$mud[,2] - varout$k[,3]*log(mn.val["tss"]) +
                           log(mn.val["tp"])), prob =credint))
    d2 <- apply(varout$d2, 2, mean)
    print(range(exp(d2 - varout$k[,3]*log(mn.val["tss"]) +
                           log(mn.val["tp"]))))
    print(quantile(varout$k[,3], prob = credint))

    print("*** volatile suspended sediment N ***")
    print(quantile(exp(varout.n$mud[,2] - varout$k[,2]*log(mn.val["tss"]) +
                           log(mn.val["tn"]*1000)), prob = credint))
    d2n <- apply(varout.n$d2, 2, mean)
    print(sort(exp(d2n)))

    print(range(exp(d2n - varout$k[,2]*log(mn.val["tss"]) +
                        log(mn.val["tn"]*1000))))
    print(quantile(varout.n$k[,2], prob = credint))

    grey.t1 <- adjustcolor("grey20", alpha = 0.5)
    grey.t2 <- adjustcolor("grey70", alpha = 0.5)

    dftemp <- unique.data.frame(df1[, c("lakenum", "flush.rate")])
    dftemp <- dftemp[order(dftemp$lakenum),]

    um <- mean(apply(varout$u, 2, mean))

    umn <- mean(apply(varout.n$u,2, mean))
    print("*** mean suspended seds ***")
    print(exp(um + log(mn.val["tss"])))
    print(exp(umn + log(mn.val["tss"])))

    predout <- matrix(NA, ncol = 3, nrow = 15)
    predout2 <- matrix(NA, ncol = 3, nrow = 15)
    for (i in 1:15) {
        y <- mn.val["tp"]*exp(varout$d2[,i])*exp(um)^varout$k[,3]/
            (exp(um)*mn.val["tss"])
        y2 <- 1000*mn.val["tn"]*exp(varout.n$d2[,i])*exp(umn)^varout.n$k[,2]/
            (exp(umn)*mn.val["tss"])
        predout[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))
        predout2[i,] <- quantile(y2, prob = c(0.025, 0.5, 0.975))
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

    plot(log(dftemp$flush.rate), predout2[,2], type = "n",axes = F,
         xlab = "Flush rate (1/yr)",
         ylab = expression(N/VSS[np]~(mu*g/m*g)), ylim = range(predout2))
    logtick.exp(0.001, 10, c(1), c(F,F))
    segments(log(dftemp$flush.rate), predout2[,1],
             log(dftemp$flush.rate), predout2[,3], col = "grey39")
    points(log(dftemp$flush.rate), predout2[,2], col = "grey39", pch = 21,
           bg = "white")
    axis(2)
    box(bty = "l")

    dev.off()
    ## plot effect of time on d1
    chl0 <- 10
    predout <- matrix(NA, ncol = 3, nrow = 6)
    predout2 <- matrix(NA, ncol = 3, nrow = 6)
    predout3 <- matrix(NA, ncol = 3, nrow = 6)
    for (i in 1:6) {
        y <- mn.val["tp"]*exp(varout$d1[,i])/(mn.val["chl"]^varout$k[,2])*chl0^(varout$k[,2]-1)
        y2 <- 1000*mn.val["tn"]*exp(varout.n$d1[,i])/(mn.val["chl"]^varout.n$k[,1])*chl0^(varout.n$k[,1]-1)

        y3 <- y2/y
        predout[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))
        predout2[i,] <- quantile(y2, prob = c(0.025, 0.5, 0.975))
        predout3[i,] <- quantile(y3, prob = c(0.025, 0.5, 0.975))
    }
    png(width = 6, height = 2, pointsize = 8, units ="in", res = 600,
        file = "timechlp.png")

    par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l", mfrow = c(1,3))
    plot(1:6, predout[,2], ylim = range(predout), type = "n", xlab = "",
         ylab = expression(P/Chl~(mu*g/mu*g)), axes= F)
    axis(2)
    axis(1, at = 1:6, lab = c("Jan/Feb", "Mar/Apr", "May/Jun",
                          "Jul/Aug", "Sep/Oct", "Nov/Dec"))
    box(bty = "l")
    segments(1:6, predout[,1], 1:6, predout[,3], col = "grey39")
    points(1:6, predout[,2], pch = 21, col = "grey39", bg = "white")

    plot(1:6, predout2[,2], ylim = range(predout2), type = "n", xlab = "",
         ylab = expression(N/Chl~(mu*g/mu*g)), axes= F)
    axis(2)
    axis(1, at = 1:6, lab = c("Jan/Feb", "Mar/Apr", "May/Jun",
                          "Jul/Aug", "Sep/Oct", "Nov/Dec"))
    box(bty = "l")
    segments(1:6, predout2[,1], 1:6, predout2[,3], col = "grey39")
    points(1:6, predout2[,2], pch = 21, col = "grey39", bg = "white")

    plot(1:6, predout3[,2], ylim = range(predout3), type = "n", xlab = "",
         ylab = expression(N/P), axes= F)
    axis(2)
    axis(1, at = 1:6, lab = c("Jan/Feb", "Mar/Apr", "May/Jun",
                          "Jul/Aug", "Sep/Oct", "Nov/Dec"))
    box(bty = "l")
    segments(1:6, predout3[,1], 1:6, predout3[,3], col = "grey39")
    points(1:6, predout3[,2], pch = 21, col = "grey39", bg = "white")


    dev.off()
    cat("N:P:", predout2[,2]/predout[,2], "\n")

    xnew <- seq(log(0.1), log(300), length = 50)
    predout <- matrix(NA, ncol = 3, nrow = length(xnew))
    predouta <- matrix(NA, ncol = 3, nrow = length(xnew))
    predout2 <- matrix(NA, ncol = 3, nrow = length(xnew))
    predout2a <- matrix(NA, ncol = 3, nrow = length(xnew))
    nsamp <- nrow(varout$mud)
    nsamp2 <- nrow(varout.n$mud)
    for (i in 1:length(xnew)) {
        y <- rnorm(nsamp, mean = varout$mud[,1], sd = varout$sigd[,1]) +
            varout$k[,2]*(xnew[i]-log(mn.val["chl"])) + log(mn.val["tp"])
        predout[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))
        y <- varout$mud[,1]+
            varout$k[,2]*(xnew[i]-log(mn.val["chl"])) + log(mn.val["tp"])
        predouta[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))
        y2 <- rnorm(nsamp2, mean = varout.n$mud[,1], sd = varout.n$sigd[,1]) +
            varout.n$k[,1]*(xnew[i] - log(mn.val["chl"])) + log(mn.val["tn"])
        predout2[i,] <- quantile(y2, prob = c(0.025, 0.5, 0.975))
        y2 <- varout.n$mud[,1] +
            varout.n$k[,1]*(xnew[i] - log(mn.val["chl"])) + log(mn.val["tn"])
        predout2a[i,] <- quantile(y2, prob = c(0.025, 0.5, 0.975))
    }

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(2,3))
    for (i in 1:6) {
        incvec <- df1$seasnum == i
        plot(log(df1$chl)+ log(mn.val["chl"]),
             log(df1$tp - df1$dtp) + log(mn.val["tp"]),
             type = "n",
             axes = F,
             xlab = expression(Chl~(mu*g/L)),
             ylab = expression(P[part]~(mu*g/L)))
        points(log(df1$chl)[incvec]+ log(mn.val["chl"]),
               log(df1$tp - df1$dtp)[incvec] + log(mn.val["tp"]),
               pch = 21, col = "grey39", bg = "white")
        logtick.exp(0.001, 10, c(1,2), c(F,F))
        lines(xnew, predout[,2])
    }
    mub <- mean(varout$mub)
    k <- apply(varout$k, 2, mean)
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(2,3))
    for (i in 1:6) {
        incvec <- df1$seasnum == i
        plot(log(df1$chl),
             log(df1$vss),
             type = "n",
             axes = F,
             xlab = expression(Chl~(mu*g/L)),
             ylab = expression(P[part]~(mu*g/L)))
        points(log(df1$chl)[incvec],
               log(df1$vss)[incvec],
               pch = 21, col = "grey39", bg = "white")
        logtick.exp(0.001, 10, c(1,2), c(F,F))
        abline(mub, k[1])
    }
    stop()

    png(width = 6, height = 2.5, pointsize = 6, units = "in",
        res = 600, file = "molimits.png")
    par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0))
    plot(log(df1$chl) + log(mn.val["chl"]),
         log(df1$tp - df1$dtp) + log(mn.val["tp"]),
         pch = 21, col = "grey39", bg = "white", axes = F,
         xlab = expression(Chl~(mu*g/L)),
         ylab = expression(P[part]~(mu*g/L)))
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    lines(xnew, predout[,2])
    polygon(c(xnew, rev(xnew)), c(predout[,1], rev(predout[,3])), col = grey.t2,
            border = NA)
    polygon(c(xnew, rev(xnew)), c(predouta[,1], rev(predouta[,3])), col = grey.t1,
            border = NA)

    tout <- approx(xnew, predout[,1], log(10))
    print(exp(tout$y))
    tout <- approx(xnew, predout[,3], log(10))
    print(exp(tout$y))


    plot(log(df1$chl) + log(mn.val["chl"]),
         log(df1$tn - df1$dtn) + log(mn.val["tn"]),
         pch = 21, col = "grey39", bg = "white", axes = F,
         xlab = expression(Chl~(mu*g/L)),
         ylab = expression(N[part]~(mu*g/L)))
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    lines(xnew, predout2[,2])
    polygon(c(xnew, rev(xnew)), c(predout2[,1], rev(predout2[,3])),col = grey.t2,
            border = NA)
    polygon(c(xnew, rev(xnew)), c(predout2a[,1], rev(predout2a[,3])), col = grey.t1,
            border = NA)
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

#varout.mo.d1T.d2L <- tss.explore(moi3.all, runmod = T, xvalid= F)
#matout.mo.d1T.d2T <-  tss.explore(moi3.all, runmod = T, xvalid= T)

tss.explore(moi3.all, matout = matout.mo.d1T.d2L, varout = varout.mo.d1T.d2L,
            varout.n = varout.mon.d1T.d2Lv,
            runmod = F, xvalid = F)
