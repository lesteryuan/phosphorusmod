## 11.6.2019: Production version of TP-chl model
## 12.17.2019Cleaned and commented
## 1.6.2021: Testing adjustment for dtp
## 1.21.2021: Limit only model for TP
ntumodel <- function(df1, varout = NULL, varout.mo = NULL,
                     varout.n = NULL, varout.mo.n = NULL,
                     runmod = T) {
#    source("logtick.exp.R")
#    source("binv.R")

    require(rstan)

    nchains <- 6    # select number of chains

    depthsel <- 3.2              # lake depth
    ecosel <- 65                # Level III ecoregion code
    credint <- 0.80             # Credible interval for criterion calc
    chltarg <- 10               # Chl target

    ## omit records that are missing data
    incvec <- ! is.na(df1$ptl.result) & ! is.na(df1$turb.result) &
         ! is.na(df1$chl) & ! is.na(df1$index.site.depth)
    cat("N omitted due to missing records:", sum(!incvec), "\n")
    df1 <- df1[incvec,]

    ## select index sites
    df1 <- subset(df1, sample.type == "MICX")

    ## drop below detection limit turb and tp
    norig <- nrow(df1)
    incvec <- df1$turb.result > 0.01
    df1 <- df1[incvec,]
    incvec <- df1$ptl.result > 1
    df1 <- df1[incvec,]
    incvec <- df1$chl >0
    df1 <- df1[incvec,]
    cat("N dropped for detection limit:", norig - nrow(df1), "\n")

    ## drop samples with chl > 100
    incvec <- df1$chl < 108
    df1 <- df1[incvec,]
    incvec <- df1$chl > 1
    df1 <- df1[incvec,]

    print(summary(df1$ptl.result))
    print(summary(df1$chl))
    print(summary(df1$ntl.result))
    print(nrow(df1))
    print(length(unique(df1$site.id)))

    ## scale chl and tp
    chlmn <- mean(log(df1$chl))
    tpsc <- exp(mean(log(df1$ptl.result)))
    chlsc <- exp(chlmn)
    df1$chl.sc <- df1$chl/chlsc
    df1$tp.sc <- df1$ptl.result/tpsc

    ## define 30 depth classes based on quantiles
    cutp.depth <- quantile(log(df1$index.site.depth), prob = seq(0, 1,length = 31))
    cutm <- 0.5*(cutp.depth[-1] + cutp.depth[-length(cutp.depth)])
    df1$dclass <- cut(log(df1$index.site.depth), cutp.depth, include.lowest = T)
    names(cutm) <- levels(df1$dclass)
    df1$dclassnum <- as.numeric(df1$dclass)

    ##drop HI, only dealing with conterminous US
    incvec <- df1$st.nla2012 == "HI"
    df1 <- df1[!incvec,]

    # thin data down by 0.25
 #   set.seed(1)
 #   isamp <- sample(nrow(df1))
 #   isamp <- isamp[1:floor(0.25*nrow(df1))]
 #   df1 <- df1[isamp,]
  #  print(nrow(df1))

    ## make seas factor
    xcut <- seq(min(df1$yday.x), by = 30, length = 6)
    xcut[6] <- max(df1$yday.x)
    hist(df1$yday.x)
    abline(v = xcut)
    seasfac <- cut(df1$yday.x, xcut, include.lowest = T)
    df1$seasnum <- as.numeric(seasfac)

    ## reset factor levels for state
    df1$state <- factor(df1$st.nla2012)
    df1$statenum <- as.numeric(df1$state)

    ## reset L3 ecoregion codes
    df1$us.l3code <- factor(df1$us.l3code)
    df1$econum <- as.numeric(df1$us.l3code)

    incvec <- df1$us.l3code == 40
    print(unique(df1$econum[incvec]))

    ## save data out to disk
    tpchldat <- df1
    save(tpchldat, chlsc, cutp.depth, file = "tpchldat.rda")


    datstan <- list(n = nrow(df1),
                    ndepth = max(df1$dclassnum),depthnum = df1$dclassnum,
                    neco = max(df1$econum), econum = df1$econum,
                    ntu = log(df1$turb.result),
                    tp = log(df1$tp.sc),
                    chl = df1$chl.sc,
                    depth = log(df1$index.site.depth),
                    nseas = max(df1$seasnum),
                    seasnum = df1$seasnum,
                    a = c(4, -1))

    print(str(datstan))


    modstan <- '
        data {
            int n;                 // number of samples
            int ndepth;            // number of depth classes
            int depthnum[n];       // depth class assignment
            int neco;              // number of L3 ecoregions
            int econum[n];         // ecoregion assigment
            vector[n] ntu;         // NTU
            vector[n] tp;          // TP
            vector[n] chl;         // Scaled Chl
            vector[n] depth;
            real a[2];
            int nseas;
            int seasnum[n];
        }
        parameters {
            real muk;  // exponents on chl and u in models
            real mud[2];
          //  vector[nseas] d2;
            real<lower = 0> sigd[3];// SD of ecoregion- or depth-specific coef
            vector[neco] etad1;
            vector[n] etad1a;
            vector[nseas] etad2;

            real<lower = 0> sigtp;    // measurement error of tp

        }
        transformed parameters {
            // dropped depth specific d1 because no sig
            // difference among depths. just using mud[1]
            vector[neco] d1;
            vector[n] d1a;
            vector[nseas] d2;

            d1 = mud[1] + sigd[1]*etad1;
            d1a = d1[econum] + sigd[3]*etad1a;
            d2 = mud[2] + sigd[2]*etad2;
        }
        model {
            vector[n] tp_mn;

            muk ~ normal(1,1);    // muk should be somewhere around 1
            mud ~ normal(0,4);
   //         d2 ~ normal(0,4);

            sigd[1] ~ cauchy(0,3);
            sigd[2] ~ cauchy(0,1);
            sigd[3] ~ cauchy(0,3);

            etad1 ~ normal(0,1);
            etad1a ~ normal(0,1);
            etad2 ~ normal(0,1);

            sigtp ~ normal(0.1, 0.002);

//           mud[1] ~ normal(a[1], 0.8);
//           slp ~ normal(a[2], 0.8);

           for (i in 1:n) {
               tp_mn[i] = exp(d1a[i]) +
                          exp(d2[seasnum[i]])*chl[i]^muk;
           }

            tp ~ student_t(4,log(tp_mn),sigtp);
        }
    '


    if (runmod) {

        rstan_options(auto_write = TRUE)
        options(mc.cores = nchains)
        fit <- stan(model_code = modstan,
                    data = datstan, iter = 2400, chains = nchains,
                    warmup = 600, thin = 3)
        return(fit)
    }


    ## post processing
    grey.t <- adjustcolor("grey39", alpha.f = 0.5)

    credint <- c(0.025, 0.5, 0.975)
    mud <- mean(varout$mud)
    muk <- mean(varout$muk)
##    d1a <- apply(varout$d1a, 2, mean)  ## loaded from the file
    d2 <- apply(varout$d2, 2, mean)
    predout <- exp(d1a) + exp(d2[df1$seasnum])*df1$chl.sc^muk
    plot(log(predout), log(df1$tp.sc))
    abline(0,1)

    rmsout <- function(x, y) sqrt(sum((x-y)^2)/length(x))
    print(rmsout(log(predout), log(df1$tp.sc)))
    stop()

    png(width = 6, height = 2.5, pointsize = 6, units = "in", res = 600,
        file = "nla.mo.comp.png")
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l", mfrow = c(1,2))
    plot(log(df1$chl), log(df1$ptl.result),
         xlab = expression(Chl~(mu*g/L)),
         ylab = expression(TP~(mu*g/L)), pch= 21, col = "grey",
         bg = "white", axes = F)

    logtick.exp(0.001, 10, c(1,2), c(F,F))

    x <- seq(min(log(0.5)), max(log(df1$chl)), length = 50)
    predout <- matrix(NA, ncol = 3, nrow = length(x))

    ## calculate estimate based on MO data
    predout.mo <- matrix(NA, ncol = 3, nrow = length(x))
    load("mn.val.mo.rda")

    nsamp <- nrow(varout$mud)

    for (i in 1:length(x)) {
#        y <- varout$mud[,3] - varout$muk[,3]*log(chlsc) +
#            varout$muk[,3]*x[i]
        y <-  rnorm(nsamp, mean = varout$mud[,2], sd = varout$sigd[,2])  +
            varout$muk*(x[i] - log(chlsc)) + log(tpsc)
        y2 <- varout.mo$d1[,4] +
            varout.mo$k[,2]*(x[i] - log(mn.val["chl"])) +
                log(mn.val["tp"])

        predout[i,] <- quantile(y, prob = credint)
        predout.mo[i,] <- quantile(y2, prob = credint)
    }

    polygon(c(x, rev(x)), c(predout.mo[,1], rev(predout.mo[,3])),
            col = grey.t, border = NA)
    lines(x, predout[,1])
    lines(x, predout[,3])

    d1.mo <- apply(varout.mo.n$d1, 2, mean)
    k.mo <- apply(varout.mo.n$k, 2, mean)
    load("mn.val.tnmo.rda")
    load("tnsc.rda")

    predout1 <- matrix(NA, ncol = 3, nrow = length(x))
    predout2 <- matrix(NA, ncol = 3, nrow = length(x))
    ns1 <- nrow(varout.mo.n$mud)

    for (i in 1:length(x)) {
        y <- rnorm(ns1, mean = varout.n$mud[,1], sd =varout.n$sigd[,1]) +
            varout.n$muk*(x[i]-log(chlsc)) + log(tnsc)
        predout1[i,] <- quantile(y, prob = credint)
        y2 <- varout.mo.n$d1[,4] + varout.mo.n$k[,1]*(x[i] - log(mn.val["chl"])) +
            log(mn.val["tn"]) + log(1000)
        predout2[i,] <- quantile(y2, prob = credint)
    }


    plot(log(df1$chl.sc*chlsc), log(df1$ntl.result - df1$no3no2.result),
         axes = F, xlab = expression(Chl~(mu*g/L)),
         ylab = expression(TN - NO[x]~(mu*g/L)),pch = 21, col = "grey",
         bg = "white")
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    polygon(c(x, rev(x)), c(predout2[,1], rev(predout2[,3])),
            col = grey.t, border= NA)

    lines(x, predout1[,1])
    lines(x, predout1[,3])
    dev.off()





}

## runmod variable set to T to run simulation and set to F to
##  run post processing.
fitout <- ntumodel(dat.merge.all, runmod = T)
## post processing
varout.p.limnat <- extract(fitout, pars = c("muk", "mud", "sigd", "d1", "d2", "d1a"))

#ntumodel(dat.merge.all, varout = varout.p.limnat,
#         varout.mo = varout.mo.d1T.d2L,
#         varout.n = varout.n.limnat,
#         varout.mo.n = varout.mon.d1T.d2Lv,
#         runmod = F)






