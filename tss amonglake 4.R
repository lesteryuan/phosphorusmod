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
#    dev.new()
#    plot(log(df1$chl), log(df1$dtp - df1$srp), axes= F,xlab="Chl",
#         ylab = "DTP - SRP")
#    logtick.exp(0.001, 10, c(1,2), c(F,F))
    mod <- lm(log(df1$dtp - df1$srp) ~ log(df1$chl))
#    abline(mod)
    print(summary(mod))

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
            real<lower = 0> sigd[3];
            vector[nlake] etad1;
            vector[nlake] etad2;
            vector[nseas] etad2a;

            real<lower = 0> sigtss;
            real<lower = 0> sigtp;

        }
        transformed parameters {
            vector[n] u;
            vector[n] tp_mn;
            vector[n] tss_mn;
            vector[nlake] d1;
            matrix[nlake, nseas] d2;

            u = muu + etau*sigu;

            d1 = mud[1] + etad1*sigd[1];
            for (i in 1:nseas) {
                d2[,i] = mud[2] + etad2*sigd[2] + etad2a[i]*sigd[3];
            }

            for (i in 1:n) {
               tss_mn[i] = exp(mub)*chl[i]^k[1] + exp(u[i]);

               tp_mn[i] = exp(d1[lakenum[i]])*exp(u[i])^k[2] +
                          exp(d2[lakenum[i], seasnum[i]])*chl[i]^k[3] + dtp[i];
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
            etad2a ~ normal(0,1);

            k[1] ~ normal(0.832,0.013);  // from VSS model
            k[2] ~ normal(1,1);
            k[3] ~ normal(1,1);

            sigtp ~ cauchy(0,3);
            sigtss ~ cauchy(0,3);

            tss ~ lognormal(log(tss_mn), sigtss);
            tp ~ lognormal(log(tp_mn), sigtp);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2)/length(x))
    gettp <- function(df, varout, withu) {
        k <- apply(varout$k, 2, mean)
        mud <- apply(varout$mud, 2, mean)
        d1 <- apply(varout$d1, 2, mean)
        d2 <- apply(varout$d2, c(2,3), mean)
        mub <- mean(varout$mub)

        if (! withu) {
            u <- df$tss - exp(mub)*df$chl^k[1]
            incvec <- u < 0
            u[incvec] <- 0
           }
        else {
            u <- exp(apply(varout$u, 2, mean))
        }

        tppred <- rep(NA, times = nrow(df))
        for (i in 1:nrow(df)) {
            tppred[i] <- exp(d1[df$lakenum[i]])*u[i]^k[2] +
                exp(d2[df$lakenum[i], df$seasnum[i]])*df$chl[i]^k[3] + df$dtp[i]
        }

        return(tppred)
    }
    extractvars <- c("mud", "k", "mub", "u", "d1", "d2", "sigd")

    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)
        nchains <- 3
        options(mc.cores = nchains)

        if (! xvalid) {
            datstan <- list(n = nrow(df1),
                            nlake = max(df1$lakenum),lakenum = df1$lakenum,
                            nseas = 4, seasnum = df1$seasnum,
                            nvss = df1$nvss,
                            tp = df1$tp,
                            dtp = df1$dtp,
                            vss = df1$vss,
                            chl = df1$chl,
                            tss = df1$tss)
            print(str(datstan))

            fit <- stan(model_code = modstan,
                        data = datstan, iter = 1100, chains = nchains,
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

    tppred <- gettp(df1, varout, withu = TRUE)
    print(summary(tppred))
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2))
    plot(log(tppred), log(df1$tp))
    abline(0,1)
    print(rmsout(log(tppred), log(df1$tp)))

    print(summary(matout))

    plot(matout[,1], matout[,2])
    points(log(tppred), log(df1$tp), pch = 16, col = "red", cex = 0.6)
    abline(0,1)



    print(rmsout(matout[,1], matout[,2]))

    return()

}
#vartss.d2tl.d1l <- tss.explore(moi3.all, runmod = T, xvalid= F)
#mattss.d2tl.d1l <-  tss.explore(moi3.all, runmod = T, xvalid= T)

tss.explore(moi3.all, matout = mattss.d2tl.d1l, varout = vartss.d2tl.d1l, runmod = F, xvalid = F)
##tss.explore(moi3.all, matout.chl, varout.chl, runmod = F, xvalid = F)
