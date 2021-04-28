# 1/23/2019: model TSS rather than nvss and vssn
# 1/31/2019: model TSS with changing coef with chl
## 12.18.2020: Model seston
## 1.4.2021: Use srp rather than dtp in model, and ntu rather than tss
## 1.5.2021: TN model
## 1.22.2021: Change back to straight particulate matter model
tnmod <- function(df1, matout = NULL,varout = NULL,
                        runmod = T, xvalid = F) {

    print(nrow(df1))
    df1 <- merge(df1, res.dat[, c("MU..", "flush.rate")], by.x = "lake",
                 by.y = "MU..")
    print(nrow(df1))
    print(summary(df1$flush.rate))

    nperiod <- 6
   # df1$yday.q <- ceiling(df1$month*0.5)
    df1$yday.q <- df1$month
    df1$yday.q <- factor(df1$yday.q)
    print(table(df1$yday.q))

    incvec <- ! is.na(df1$don) & ! is.na(df1$chl) & ! is.na(df1$din) &
        ! is.na(df1$tn) & ! is.na(df1$doc) & ! is.na(df1$vss)
    df1 <- df1[incvec,]
    df1$lake <- factor(df1$lake)

    df1$lakenum <- as.numeric(df1$lake)
    df1$seasnum <- as.numeric(df1$yday.q)

    varlist<- c("tn", "chl", "doc", "tss")
    mn.val <- apply(df1[, varlist],2,function(x) exp(mean(log(x))))
    save(mn.val, file = "mn.val.tnmo.rda")


    for (i in varlist) df1[,i] <- df1[,i]/mn.val[i]
    df1$din <- df1$din/mn.val["tn"]
    df1$don <- df1$don/mn.val["tn"]
    df1$dtn <- df1$dtn/mn.val["tn"]
    df1$vss <- df1$vss/mn.val["tss"]

    ## drop pn measurements near zero or negative
    pn <- df1$tn - df1$din - df1$don
    incvec <- pn < 1e-10
    print(sum(incvec))

    df1 <- df1[!incvec,]


       modstan <- '
        data {
            int n;
            int nlake;
            int lakenum[n];
            int nseas;
            int seasnum[n];
            vector[n] tn;
            vector[n] dtn;
            vector[n] chl;
            vector[n] tss;
        }
        parameters {
            real k[3];
            vector[2] mud;
            real<lower = 0> sigd[2];
            vector[nlake] etad1;
            vector[nseas] etad2;

            real<lower = 0> sigtn;

            real muu;
            real<lower =0> sigu;
            vector[n] etau;

            real mub;
            real<lower = 0> sigb;
            vector[nseas] etab;
            real<lower = 0> sigtss;

        }
        transformed parameters {
            vector[n] tn_mn;
            vector[n] u;
            vector[n] tss_mn;

            vector[nseas] b;
            vector[nlake] d1;
            vector[nseas] d2;

            b = mub + sigb*etab;

            d1 = mud[1] + etad1*sigd[1];
            d2 = mud[2] + etad2*sigd[2];

            u = muu + sigu*etau;

            for (i in 1:n) {
                tss_mn[i] = log_sum_exp(b[seasnum[i]] + k[3]*chl[i], u[i]);

                tn_mn[i] = log_sum_exp(mud[1] + k[1]*chl[i],
                        d2[seasnum[i]] + k[2]*u[i]);
            }
        }
        model {
            mud ~ normal(0,3);
            sigd ~ cauchy(0,3);
            etad1 ~ normal(0,1);
            etad2 ~ normal(0,1);

            k[3] ~ normal(0.85,0.01);
            k[1] ~ normal(1,1);
            k[2] ~ normal(1,1);
            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);
            sigtss ~ normal(0.1,0.002);

            mub ~ normal(0,3);
            etab ~ normal(0,1);
            sigb ~ cauchy(0,3);

            sigtn ~ cauchy(0,3);

            tss ~ student_t(4,tss_mn, sigtss);
            tn ~ normal(tn_mn, sigtn);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2, na.rm = T)/
                                     min(sum(! is.na(x)), sum(!is.na(y))))
    gettn <- function(df, varout, withu=T) {
        k <- apply(varout$k,2,mean)
        mud <- apply(varout$mud, 2, mean)
        mub <- mean(varout$mub)
        d1 <- apply(varout$d1, 2, mean)
        d2 <- apply(varout$d2, 2, mean)


        if (withu) {
            u <- exp(apply(varout$u, 2, mean))
        }
        else {
            u <- df$vss - exp(mub)*df$chl^k[3]
            u[u <= 0] <- 0
        }

        tnpred <- exp(mud[1])*df$chl^k[1] +
            exp(d2[df$lakenum])*u^k[2]
        return(tnpred)
    }
    extractvars <- c("mud", "k", "u", "mub", "d1", "d2", "sigd")

    nchains <- 3
    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)

        options(mc.cores = nchains)

        if (! xvalid) {
            datstan <- list(n = nrow(df1),
                            nlake = max(df1$lakenum),lakenum = df1$lakenum,
                            nseas = max(df1$seasnum), seasnum = df1$seasnum,
                            tn = log(df1$tn-df1$dtn),
                            dtn = df1$dtn,
                            chl = log(df1$chl),
                            tss = log(df1$vss))
            print(str(datstan))

            fit <- stan(model_code = modstan,
                        data = datstan, iter = 1800, chains = nchains,
                        warmup = 600, thin= 1,
                        control = list(adapt_delta = 0.98, max_treedepth = 14))
            save(fit, file = "fitout.rda")
            varout <- extract(fit, pars = extractvars)

            tnpred <- gettn(df1, varout)
            dev.new()
            plot(log(tnpred), log(df1$tn-df1$dtn))
            abline(0,1)
            cat("Internal RMS:", rmsout(log(tnpred), log(df1$tn-df1$dtn)), "\n")

            return(varout)
        }
        else {

            set.seed(2)
            require(loo)
            nfold <- 5
            ik <- kfold_split_random(nfold, nrow(df1))

            for (jj in 1:nfold) {

                dftemp2 <- df1[ik == jj,]
                dftemp <- df1[ik != jj,]
                datstan <- list(n = nrow(dftemp),
                                nlake = max(dftemp$lakenum),
                                lakenum = dftemp$lakenum,
                                nseas = 6, seasnum = dftemp$seasnum,
                                tn = log(dftemp$tn),
                                dtn = dftemp$dtn,
                                chl = dftemp$chl,
                                tss = log(dftemp$vss))

                fit <- stan(model_code = modstan,
                            data = datstan, iter = 1000, chains = nchains,
                            warmup = 400, thin= 2,
                            control = list(adapt_delta = 0.98, max_treedepth = 14))

                varout <- extract(fit, pars = extractvars)
                tnpred <- gettn(dftemp2, varout, withu = FALSE)

                matval <- data.frame(pred = log(tnpred), obs = log(dftemp2$tn))

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
    plot(log(dat.merge.all$chl), log(dat.merge.all$ntl.result -
                                         dat.merge.all$no3no2.result),
         col = "grey")
    points(log(df1$chl*mn.val["chl"]), log(1000*(df1$tn-df1$dtn)*mn.val["tn"]),
           pch = 16, cex = 0.5)
    stop()

    dftemp <- df1
    dftemp <- dftemp[dftemp$doc > 0.5,]

    cutp <- quantile(dftemp$chl, prob = seq(0, 1, length = 6))
    cutf <- cut(dftemp$chl, cutp, include.lowest = T)
    dflist <- split(dftemp, cutf)
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(2,3))
    plot(log(dftemp$doc), log(dftemp$don), col = "grey")
    for (i in 1:length(dflist)) {

#        points(log(dflist[[i]]$doc), log(dflist[[i]]$don), pch = 16, col = "blue")
        mod <- lm(log(dflist[[i]]$don) ~ log(dflist[[i]]$doc))
        print(summary(mod))
        abline(mod)
    }
    stop()


    tnpred <- gettn(df1, varout)

    cat("Internal RMS:", rmsout(log(tnpred), log(df1$tn-df1$dtn)), "\n")
    stop()

    print(quantile(exp(varout$mud[,1] - varout$k[,1]*log(mn.val["chl"]) +
                           log(mn.val["tn"]*1000)), prob = c(0.05, 0.5, 0.95)))
    print(quantile(varout$k[,1], prob = c(0.05, 0.5, 0.95)))
    stop()
    k <- apply(varout$k, 2, mean)
    mub <- mean(varout$mub)
    u <- apply(varout$u, 2, mean)

    d1 <- apply(varout$d1, 2, quantile, prob = c(0.05, 0.5, 0.95))
    d2 <- apply(varout$d2, 2, quantile, prob = c(0.05, 0.5, 0.95))

    plot(log(df1$chl), log(exp(d2[df1$lakenum])*exp(u)^k[2]))

    dev.new()
    ucalc <- df1$vss - exp(mub)*df1$chl^k[3]
    plot(ucalc, exp(u))
    abline(0,1)

    tsspred <- exp(mub)*df1$chl^k[3] + exp(u)
    plot(log(df1$chl), log(df1$vss))
    abline(mub, k[3])

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2))
    tnpred <- gettn(df1, varout, withu = T)
    plot(log(df1$chl), log(df1$tn)- log(tnpred))
    abline(h=0)

    cat("Internal RMS:", rmsout(log(tnpred), log(df1$tn)), "\n")
    cat("External RMS:", rmsout(matout[,1], matout[,2]), "\n")
    plot(matout[,1], matout[,2])
    abline(0,1)

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2))

    nd1 <- ncol(d1)
    plot(1:nd1, d1[2,], ylim = range(d1))
    segments(1:nd1, d1[1,], 1:nd1, d1[3,])
    nd2 <- ncol(d2)
    plot(1:nd2, d2[2,], ylim = range(d2))
    segments(1:nd2, d2[1,], 1:nd2, d2[3,])

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(2,3))
    for (i in 1:6) {
        incvec <- df1$seasnum == i
        plot(log(df1$chl), log(df1$tn - df1$dtn), type = "n")
        points(log(df1$chl)[incvec], log(df1$tn - df1$dtn)[incvec], pch = 21)
        abline(d1[2,i], k[1])
    }

    mud <- apply(varout$mud, 2, mean)

    abline(mud[1], k[1])

    plot(log(dat.merge.all$chl/mn.val["chl"]),
         log((dat.merge.all$ntl.result-dat.merge.all$no3no2.result)/
                 (1000*mn.val["tn"])))
    abline(mud[1], k[1])



    stop()
    print(summary(matout))
    dev.new()
    plot(matout[,1], matout[,2])
    points(log(tppred), log(df1$tp), pch = 16, col = "red", cex = 0.6)
    abline(0,1)



    print(rmsout(matout[,1], matout[,2]))

    return()

}
varout.test<- tnmod(moi3.all, runmod = T, xvalid = F)
#matout.mon.d1T.d2Tv <- tnmod(moi3.all, runmod = T, xvalid = T)
#tnmod(moi3.all, matout = matout.mon.d1T.d20,
#      varout = varout.mon.d10.d2L, runmod = F, xvalid = F)
