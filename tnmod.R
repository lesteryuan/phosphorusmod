# 1/23/2019: model TSS rather than nvss and vssn
# 1/31/2019: model TSS with changing coef with chl
## 12.18.2020: Model seston
## 1.4.2021: Use srp rather than dtp in model, and ntu rather than tss
## 1.5.2021: TN model
tnmod <- function(df1, matout = NULL,varout = NULL,
                        runmod = T, xvalid = F) {

    print(nrow(df1))
    df1 <- merge(df1, res.dat[, c("MU..", "flush.rate")], by.x = "lake",
                 by.y = "MU..")
    print(nrow(df1))
    print(summary(df1$flush.rate))

    nperiod <- 4
    df1$yday.q <- 1
    incvec <- df1$month >= 3 & df1$month <= 5
    df1$yday.q[incvec]<- 2
    incvec <- df1$month >= 6 & df1$month <= 8
    df1$yday.q[incvec]<- 3
    incvec <- df1$month >= 9 & df1$month <= 11
    df1$yday.q[incvec]<- 4
    df1$yday.q <- factor(df1$yday.q)

    incvec <- ! is.na(df1$don) & ! is.na(df1$chl) & ! is.na(df1$din) &
        ! is.na(df1$tn) & ! is.na(df1$doc)
    df1 <- df1[incvec,]
    df1$lake <- factor(df1$lake)
    print(table(df1$lake))
    print(df1[1:10, c("tn", "din", "don")])
    df1$lakenum <- as.numeric(df1$lake)
    df1$seasnum <- as.numeric(df1$yday.q)

    varlist<- c("tn", "chl", "doc")
    mn.val <- apply(df1[, varlist],2,function(x) exp(mean(log(x))))
    save(mn.val, file = "mn.val.tnmo.rda")

    for (i in varlist) df1[,i] <- df1[,i]/mn.val[i]
    df1$din <- df1$din/mn.val["tn"]
    df1$don <- df1$don/mn.val["tn"]
    df1$dtn <- df1$dtn/mn.val["tn"]

    ## drop one big outlier
    print(summary(df1$chl))
    print(summary(df1$tn))
    print(summary(df1$din))

    ## drop pn measurements near zero or negative
    pn <- df1$tn - df1$din - df1$don
    incvec <- pn < 1e-10
    print(sum(incvec))
    df1 <- df1[!incvec,]

    plot(log(df1$doc), log(df1$don))

       modstan <- '
        data {
            int n;
            int nlake;
            int lakenum[n];
            vector[n] tn;
            vector[n] don;
            vector[n] din;
            vector[n] chl;
        }
        parameters {
            real k;
            real mud;
//            real<lower = 0> sigd;
//            vector[nlake] etad1;

            real<lower = 0> sigtn;

        }
        transformed parameters {
            vector[n] tn_mn;
//            vector[nlake] d1;
//            d1 = mud[1] + etad1*sigd[1];
            for (i in 1:n) {
               tn_mn[i] = din[i] + don[i] + exp(mud)*chl[i]^k;
            }
        }
        model {
            mud ~ normal(0,3);
//            sigd ~ cauchy(0,3);
//            etad1 ~ normal(0,1);
            k ~ normal(1,1);
            sigtn ~ cauchy(0,3);

            tn ~ lognormal(log(tn_mn), sigtn);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2)/length(x))
    gettn <- function(df, varout) {
        k <- mean(varout$k)
        mud <- mean(varout$mud)
 #       d1 <- apply(varout$d1, 2, mean)

        tnpred <- df$din + df$don +
            exp(mud)*df$chl^k

        return(tnpred)
    }
    extractvars <- c("mud", "k")

    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 2
        options(mc.cores = nchains)

        if (! xvalid) {
            datstan <- list(n = nrow(df1),
                            nlake = max(df1$lakenum),lakenum = df1$lakenum,
                            nseas = nperiod, seasnum = df1$seasnum,
                            tn = df1$tn,
                            din = df1$din,
                            don = df1$don,
                            chl = df1$chl)
            print(str(datstan))

            fit <- stan(model_code = modstan,
                        data = datstan, iter = 600, chains = nchains,
                        warmup = 300, thin= 1,
                        control = list(adapt_delta = 0.98, max_treedepth = 14))
            varout <- extract(fit, pars = extractvars)
            tnpred <- gettp(df1, varout)
            dev.new()
            plot(log(tnpred), log(df1$tn))
            abline(0,1)
            print(rmsout(log(tnpred), log(df1$tn)))

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
                                ntu = dftemp$ntu,
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
    tnpred <- gettn(df1, varout)
    plot(log(df1$tn), log(tnpred))
    abline(0,1)
    cat("Internal RMS:", rmsout(log(tppred), log(df1$tp)), "\n")

    print(summary(matout))
    dev.new()
    plot(matout[,1], matout[,2])
    points(log(tppred), log(df1$tp), pch = 16, col = "red", cex = 0.6)
    abline(0,1)



    print(rmsout(matout[,1], matout[,2]))

    return()

}
tnmod(moi3.all, runmod = T, xvalid = F)
