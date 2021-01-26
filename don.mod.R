# 1/23/2019: model TSS rather than nvss and vssn
# 1/31/2019: model TSS with changing coef with chl
## 12.18.2020: Model seston
## 1.4.2021: Use srp rather than dtp in model, and ntu rather than tss
## 1.5.2021: TN model
## 1.22.2021: Change back to straight particulate matter model
## 1.25.2021: model for don
donmod <- function(df1,varout = NULL,
                        runmod = T) {

    print(nrow(df1))
    df1 <- merge(df1, res.dat[, c("MU..", "flush.rate")], by.x = "lake",
                 by.y = "MU..")
    print(nrow(df1))
    print(summary(df1$flush.rate))

    nperiod <- 6
    df1$yday.q <- ceiling(df1$month*0.5)
    df1$yday.q <- factor(df1$yday.q)
    print(table(df1$yday.q))

    incvec <- ! is.na(df1$don) & ! is.na(df1$chl) &
        ! is.na(df1$tn) & ! is.na(df1$doc)
    df1 <- df1[incvec,]
    df1$lake <- factor(df1$lake)

    df1$lakenum <- as.numeric(df1$lake)
    df1$seasnum <- as.numeric(df1$yday.q)

    varlist<- c("tn", "chl", "doc")
    mn.val <- apply(df1[, varlist],2,function(x) exp(mean(log(x))))
    save(mn.val, file = "mn.val.tnmo.rda")

    for (i in varlist) df1[,i] <- df1[,i]/mn.val[i]
    df1$don <- df1$don/mn.val["tn"]

    ## drop pn measurements near zero or negative
    pn <- df1$tn - df1$din - df1$don
    incvec <- pn < 1e-10
    print(sum(incvec))

    df1 <- df1[!incvec,]

       modstan <- '
        data {
            int n;
            vector[n] don;
            vector[n] doc;
            vector[n] chl;
        }
        parameters {
      //      real f;
            real g[2];
            real k[2];
      //      real<lower = 0> sigdoc;
            real<lower = 0> sigdon;

       //     real muu;
       //     real<lower = 0> sigu;
       //     vector[n] etau;
        }
        transformed parameters {
       //     vector[n] u;
       //     vector[n] doc_mn;
            vector[n] don_mn;
       //     u = muu + etau*sigu;

            for (i in 1:n) {
       //         doc_mn[i] = exp(f)*chl[i]^k[1] + exp(u[i]);
                don_mn[i] = exp(g[1])*chl[i]^k[1] + exp(g[2])*doc[i]^k[2];
            }
        }
        model {
      //      f ~ normal(0,4);
            g ~ normal(0,4);
            k ~ normal(1,1);
      //      sigdoc ~ normal(0.1, 0.01);
            sigdon ~ normal(0.1, 0.01);
      //      muu ~ normal(0,4);
      //      sigu ~ cauchy(0,4);
      //      etau ~ normal(0,1);

      //      doc ~ student_t(4,log(doc_mn), sigdoc);
            don ~ student_t(4,log(don_mn), sigdon);
        }
    '
    rmsout <- function(x,y) sqrt(sum((x-y)^2, na.rm = T)/
                                     min(sum(! is.na(x)), sum(!is.na(y))))

    extractvars <- c("g", "k")

    nchains <- 3
    if (runmod) {
        require(rstan)
        rstan_options(auto_write = TRUE)

        options(mc.cores = nchains)

        datstan <- list(n = nrow(df1),
                        don = log(df1$don),
                        doc = df1$doc,
                        chl = df1$chl)

        print(str(datstan))

        fit <- stan(model_code = modstan,
                    data = datstan, iter = 1000, chains = nchains,
                    warmup = 400, thin= 1,
                    control = list(adapt_delta = 0.98, max_treedepth = 14))
        save(fit, file = "fitout.rda")
        varout <- extract(fit, pars = extractvars)
        return(varout)
    }
    else {
        dev.new()
        par(mar = c(4,4,1,1), mfrow = c(1,2))
        plot(log(df1$chl), log(df1$don))
        g <- apply(varout$g, 2, mean)

        k <- apply(varout$k, 2, mean)
        abline(g[1], k[1])
        plot(log(df1$doc), log(df1$don))
        abline(g[2], k[2])


        predout <- exp(g[1])*df1$chl^k[1] + exp(g[2])*df1$doc^k[2]

        dev.new()
        plot(log(predout), log(df1$don))
        abline(0,1)
        rmsout <- function(x,y) sqrt(sum((x-y)^2)/length(x))
        print(rmsout(log(predout), log(df1$don)))

        mod <- lm(log(df1$don) ~ log(df1$doc) + log(df1$chl))
        pred2 <- predict(mod)
        print(rmsout(predict(mod), log(df1$don)))
        points(pred2, log(df1$don), pch = 16)
    }


    return()

}
#varout <- donmod(moi3.all, runmod = T)
donmod(moi3.all, varout, runmod = F)
