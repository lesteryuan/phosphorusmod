## 11.6.2019: Production version of TP-chl model
## 12.17.2019Cleaned and commented
## 1.6.2021: Testing adjustment for dtp
ntumodel <- function(df1, varout = NULL, varout.mo = NULL, runmod = T) {
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
    incvec <- df1$chl < 100
    df1 <- df1[incvec,]
    incvec <- df1$chl > 1
    df1 <- df1[incvec,]

    ## scale chl
    chlmn <- mean(log(df1$chl))
    chlsc <- exp(chlmn)
    df1$chl.sc <- df1$chl/chlsc

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

    sdprior1 <- c(rep(1, times = 5), rep(0.3, times = 16),
                  rep(1, times = 9))


    datstan <- list(n = nrow(df1),
                    ndepth = max(df1$dclassnum),depthnum = df1$dclassnum,
                    neco = max(df1$econum), econum = df1$econum,
                    ntu = log(df1$turb.result),
                    tp = log(df1$ptl.result),
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
            real muu_mn;             // grand mean ntu_np
            real<lower= 0> sigu[2];  // standard deviation of ntu_np
            vector[ndepth] eta_u1;
            vector[n] eta_u2;
            real<lower = 0> signtu;  // measurement error of ntu
            real mub;                // mean coef of chl-ntu relationship
            real<lower = 0> sigb;    // standard deviation of b among ecoregions
            vector[neco] etab;
            real<lower = 0> muk[3];  // exponents on chl and u in models
            vector[3] mud;           // mean coef for tp model
            real<lower = 0> sigd[3];// SD of ecoregion- or depth-specific coef

            vector[n] etad1;
            vector[neco] etad2;
            real<lower = 0> sigtp;    // measurement error of tp

            real am[2];

            vector[nseas] etad3;


        }
        transformed parameters {
            vector[ndepth] muu;
            vector[neco] b;
            vector[n] u;
            // dropped depth specific d1 because no sig
            // difference among depths. just using mud[1]
            vector[n] d1s;
            vector[neco]  d2;
            vector[nseas] d3;

            b = mub + etab*sigb;

            muu = muu_mn + eta_u1*sigu[1];
            u = muu[depthnum] + eta_u2*sigu[2];

            d1s = am[1] + am[2]*depth + sigd[1]*etad1;
            d2 = mud[2] + sigd[2]*etad2;
            d3 = mud[3] + sigd[3]*etad3;
        }
        model {
            vector[n] turb_mn;
            vector[n] tp_mn;

            muu_mn ~ normal(0,3);    // weakly informative priors
            sigu ~ cauchy(0,3);
            eta_u1 ~ normal(0,1);
            eta_u2 ~ normal(0,1);
            signtu ~ normal(0.1,0.002);
            mub ~ normal(0,4);
            sigb ~ cauchy(0,3);
            etab ~ normal(0,1);

            muk ~ normal(1,0.3);    // muk should be somewhere around 1
            mud ~ normal(0,4);

            sigd ~ cauchy(0,3);

            etad1 ~ normal(0,1);
            etad2 ~ normal(0,1);
            etad3 ~ normal(0,1);

            sigtp ~ normal(0.123, 0.002);

           am[1] ~ normal(a[1], 0.4);
           am[2] ~ normal(a[2], 0.4);

           for (i in 1:n) {
               turb_mn[i] = exp(b[econum[i]])*chl[i]^muk[1] + exp(u[i]);
               tp_mn[i] = exp(d1s[i]) +
      //                    exp(d2[econum[i]])*exp(u[i])^muk[2] +
                          exp(d3[seasnum[i]])*chl[i]^muk[3];
           }

            // Eqn 25 from document
            ntu ~ student_t(4,log(turb_mn), signtu);
            // Eqn 30 from document
            tp ~ student_t(4,log(tp_mn),sigtp);
        }
    '


    if (runmod) {

        rstan_options(auto_write = TRUE)
        options(mc.cores = nchains)
        fit <- stan(model_code = modstan,
                    data = datstan, iter = 1800, chains = nchains,
                    warmup = 600, thin = 3)
        return(fit)
    }


    ## post processing
    grey.t <- adjustcolor("grey39", alpha.f = 0.5)

    b <- apply(varout$b, 2, mean)
    mud <- apply(varout$mud, 2, mean)
    muk <- apply(varout$muk, 2, mean)
    d2 <- apply(varout$d2, 2, mean)
    uup <- apply(varout$u, 2, quantile, prob = 0.05)
    udn <- apply(varout$u, 2, quantile, prob = 0.95)
    umean <- apply(varout$u, 2, mean)
    iord <- order(uup)
    ucalc <- df1$turb.result - exp(b[df1$econum])*df1$chl.sc^muk[1]
    print(length(ucalc))

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2))
    plot(umean, log(ucalc))
    abline(0,1)
    d1s <- apply(varout$d1s, 2, mean)

    plot(log(df1$index.site.depth), d1s)
    abline(4, -1)


    predout <- exp(d1s) + exp(d2[df1$econum])*exp(umean)^muk[2] +
        exp(mud[3])*df1$chl.sc^muk[3]
    plot(log(predout), log(df1$ptl.result))
    print(sqrt(sum((log(predout) - log(df1$ptl.result))^2)/length(predout)))
    abline(0,1)


    dev.new()
    plot(log(df1$chl), log(df1$ptl.result),
         xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = expression(TP~(mu*g/L)), pch= 21, col = "grey",
         bg = "white", axes = F)
#    points(        log(moi3.all$chl), log(moi3.all$tp - moi3.all$dtp),
#           pch = "+")
    logtick.exp(0.001, 10, c(1,2), c(F,F))

    x <- seq(min(log(df1$chl)), max(log(df1$chl)), length = 40)
    predout <- matrix(NA, ncol = 3, nrow = length(x))

    ## calculate estimate based on MO data
    predout.mo <- matrix(NA, ncol = 3, nrow = length(x))
    load("mn.val.mo.rda")

    nsamp <- nrow(varout$mud)

    for (i in 1:length(x)) {
#        y <- varout$mud[,3] - varout$muk[,3]*log(chlsc) +
#            varout$muk[,3]*x[i]
        y <-  varout$mud[,3]  -
            varout$muk[,3]*log(chlsc) +
                varout$muk[,3]*x[i]
        mnval <- varout.mo$d2[,4]
        y2 <- log(mn.val["tp"]) +
             mnval -
            varout.mo$k[,3]*log(mn.val["chl"]) +
                varout.mo$k[,3]*x[i]
        predout[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))
        predout.mo[i,] <- quantile(y2, prob = c(0.025, 0.5, 0.975))
    }
    print(predout)
    polygon(c(x, rev(x)), c(predout[,1], rev(predout[,3])),
            col = grey.t, border = NA)
    lines(x, predout[,2])
    lines(x, predout.mo[,2], lty = "dashed")
    lines(x, predout.mo[,1], lty = "dotted")
    lines(x, predout.mo[,3], lty = "dotted")
    stop()

    muk <- apply(varout$muk, 2, mean)
    d1 <- apply(varout$d1, 2, mean)
    d2 <- apply(varout$d2, 2, mean)
    d3 <- apply(varout$d3, 2, mean)
    muu <- apply(varout$muu, 2,mean)
    muumn <- mean(varout$muu_mn)
    d3raw <- d3 - muk[3]*log(chlsc)
    d2raw <- d2-muk[2]*muumn
    mud <- apply(varout$mud, 2, mean)

    mub <- mean(varout$mub)
    mubraw <- mub - muk[1]*log(chlsc)

    ## PLOT: relationship between Chl and turbidity (Fig 24)
    png(width = 3, height = 2.5, pointsize = 8, units = "in", res = 600,
        file = "chlturb.png")
    par(mar = c(4,4,1,1), mfrow = c(1,1), mgp = c(2.3,1,0))
    plot(log(df1$chl), log(df1$turb.result), axes = F,
         xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = "Turbidity (NTU)", pch = 21, col = "grey",
         bg="white")
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    abline(mubraw, muk[1])
    dev.off()

    b <- apply(varout$b, 2, mean)
    braw <- b - muk[1]*log(chlsc)

    ## PLOT: predicted values for Pdiss and muu (Fig 25)
    png(width = 6, height = 2.5, pointsize = 8, units = "in", res = 600,
        file = "Pdiss.turb.depth.png")
    par(mar = c(4,4,1,1), mgp = c(2.3, 1, 0), mfrow= c(1,2))
    plot(cutm, exp(muu), xlab = "Depth (m)", axes = F,
         ylab = expression(Turb[np]~(NTU)),
         pch = 21, col = "grey39", bg="white")
    axis(2)
    logtick.exp(0.001, 10, c(1), c(F,F))

    plot(cutm, exp(d1), xlab = "Depth (m)", axes = F,
         ylab = expression(P[diss]~(mu*g/L)),
         pch = 21, col = "grey39", bg="white")
    axis(2)
    logtick.exp(0.001, 10, c(1), c(F,F))
    dev.off()

    ## load in ntu_np mean values into main data
    df1$umean <- umean

    ## merge ecoregion and depth specific coefficients into main data
    dfd <- data.frame(num = 1:length(muu),  muu)
    names(dfd) <- c("dclassnum",  "muu")
    print(nrow(df1))
    df1 <- merge(df1, dfd, by = "dclassnum")
    print(nrow(df1))
    dfd3 <- data.frame(num = 1:length(d3), d3, d3raw, d2, d2raw)
    print(dim(dfd3))
    names(dfd3) <- c("econum", "d3", "d3raw", "d2", "d2raw")
    df1 <- merge(df1, dfd3, by = "econum")

    ## save ecoregion coefficients to file for mapping
    dfd3 <- merge(dfd3, unique.data.frame(df1[, c("econum", "us.l3code")]),
                  by = "econum")
    save(dfd3, file= "dfd3.rda")  # output to different script for mapping

        ## compute mean predicted TP
    df1$predout <- exp(mud[1]) + exp(df1$d2)*exp(df1$umean - muumn)^muk[2] +
        exp(df1$d3)*df1$chl.sc^muk[3]


    png(width = 4, height = 4, pointsize = 10, units = "in", res = 600,
        file = "pred.v.obs.all.png")
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0), bty = "l")
    plot(log(df1$predout), log(df1$ptl.result), xlab = "Predicted TP",
         ylab = "Observed TP",  pch = 21, col = "grey39",
         bg = "white", axes = F)
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    abline(0,1)
    dev.off()

    rmsfnc <- function(x,y) {
        return(sqrt(sum((x-y)^2)/length(x)))
    }
    cat("Whole data RMS:",rmsfnc(log(df1$predout), log(df1$ptl.result)), "\n")

    ## print numerical summaries of parameters
    print("mub")
    print(exp(quantile(varout$mub - varout$muk[,1]*log(chlsc),
                   prob = c(0.05, 0.5, 0.95))))
    print("k distribution")
    print(apply(varout$muk, 2, quantile, prob = c(0.05, 0.5, 0.95)))
    mud <- apply(varout$mud, 2, mean)

    print("mud[2]")
    print(exp(quantile(varout$mud[,2] - varout$muk[,2]*varout$muu_mn,
                       prob = c(0.05, 0.50, 0.95))))
    print("mud[3]")
    print(exp(quantile(varout$mud[,3] - varout$muk[,3]*log(chlsc),
                       prob = c(0.05, 0.5, 0.95))))

    # PLOT: NTU_np vs TP, Chl vs TP (fig 27)
#    png(width = 6, height = 2.5, pointsize = 8, units = "in", res = 600,
#        file = "tpchl.png")

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0))
    plot(df1$umean, log(df1$ptl.result), xlab = expression(Turb[np]~(NTU)),
         ylab = expression(TP~(mu*g/L)), pch = 21, col = "grey",
         bg="white", axes = F)
    logtick.exp(0.000001, 10, c(1,2), c(F,F))
    x <- seq(min(df1$umean), max(df1$umean), length = 40)
    predout <- matrix(NA, ncol = 3, nrow = length(x))
    for (i in 1:length(x)) {
        y <- varout$mud[,2] + varout$muk[,2]*x[i] - varout$muu_mn*varout$muk[,2]
        predout[i,] <- quantile(y, prob = c(0.05, 0.5, 0.95))
    }
    polygon(c(x, rev(x)), c(predout[,1], rev(predout[,3])),
            col = grey.t, border = NA)
    lines(x, predout[,2])

    plot(log(df1$chl), log(df1$ptl.result),
         xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = expression(TP~(mu*g/L)), pch= 21, col = "grey",
         bg = "white", axes = F)
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    x <- seq(min(log(df1$chl)), max(log(df1$chl)), length = 40)
    predout <- matrix(NA, ncol = 3, nrow = length(x))

    ## calculate estimate based on MO data
    predout.mo <- matrix(NA, ncol = 3, nrow = length(x))
    load("mn.val.mo.rda")
    load("varout.vss.rda")
    load("varntu.00.rda")

    d3 <- apply(varout$d3, 2, mean)
#    ip <- which(d3 == min(d3))
    ip <- 40

    for (i in 1:length(x)) {
#        y <- varout$mud[,3] - varout$muk[,3]*log(chlsc) +
#            varout$muk[,3]*x[i]
        y <- varout$d3[,ip] - varout$muk[,3]*log(chlsc) +
            varout$muk[,3]*x[i]
        y2 <- log(mn.val["tp"]) + varntu.00$mud[,2] -
            varntu.00$k[,3]*log(mn.val["chl"]) +
                varntu.00$k[,3]*x[i] + log(2)
        predout[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))
        predout.mo[i,] <- quantile(y2, prob = c(0.025, 0.5, 0.975))
    }
    polygon(c(x, rev(x)), c(predout[,1], rev(predout[,3])),
            col = grey.t, border = NA)
    lines(x, predout[,2])
    lines(x, predout.mo[,2], lty = "dashed")
    lines(x, predout.mo[,1], lty = "dotted")
    lines(x, predout.mo[,3], lty = "dotted")
    stop()

    dev.off()

    ## Plot criterion derivation figure

    ## find ecoregion index number
    ieco <- which(levels(df1$us.l3code) == ecosel)

    ## find correct depth class
    idepth <- 1
    while(exp(cutp.depth[idepth]) < depthsel & (idepth < length(cutp.depth)))
        idepth <- idepth + 1
    print(idepth)

    grey.t <- adjustcolor("grey39", alpha.f = 0.5)

    ## calculate mean parameter values
    muk <- apply(varout$muk, 2, mean)
    d3 <- apply(varout$d3, 2, mean)
    d3 <- d3 - muk[3]*log(chlsc)

    incvec <- df1$econum == ieco

    ## define regularly spaced values along chl gradient
    ## for computing predictions
    x <- seq(min(log(df1$chl)), max(log(df1$chl)), length = 50)
    xsc <- x - log(chlsc)

    ## compute predicted ambient and limiting values
    predout <- matrix(NA, ncol = 3, nrow = length(x))
    predout2 <- matrix(NA, ncol = 3, nrow = length(x))
    for (i in 1:length(x)) {
        y <- varout$d3[,ieco] - varout$muk[,3]*log(chlsc) +
            varout$muk[,3]*x[i]
        predout[i,] <- quantile(y, prob = c(0.5*(1-credint),
                                       0.5, 1- 0.5*(1-credint)))
        y2 <- exp(varout$mud[,1]) +
            exp(varout$d2[,ieco])*
                (exp(varout$muu[, idepth] -
                         varout$muu_mn))^varout$muk[,2] +
                             exp(varout$d3[, ieco])*exp(xsc[i])^varout$muk[,3]
        predout2[i,] <- quantile(log(y2), prob = c(0.5*(1-credint),
                                              0.5, 1-0.5*(1-credint)))
    }

    dev.new()
    par(mgp = c(2.3,1,0), bty = "l", mar = c(4,4,1,1))
    plot(log(df1$chl), log(df1$ptl.result), type = "p",
         axes = F, xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = expression(TP~(mu*g/L)), col = "grey80",
         pch = 21, bg = "white")
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    ## highlight points in the selected ecoregion
    points(log(df1$chl)[incvec], log(df1$ptl.result)[incvec],
           pch = 21, col = "black", bg = "grey")

    cat("Median depth:", median(df1$index.site.depth[incvec]), "\n")

    polygon(c(x, rev(x)), c(predout[,1], rev(predout[,3])),
            col  = grey.t, border = NA)
    lines(x, predout[,2])
    polygon(c(x, rev(x)), c(predout2[,1], rev(predout2[,3])),
            col  = grey.t, border = NA)
    lines(x, predout2[,2], lty = "dashed")

    crit1 <- approx(x, predout[,1], log(chltarg))$y
    crit2 <- approx(x, predout2[,1], log(chltarg))$y

    ylo <- min(log(df1$ptl.result)) - 0.04*diff(range(log(df1$ptl.result)))
    xlo <- min(log(df1$chl)) - 0.04*diff(range(log(df1$chl)))
    segments(xlo, crit1,
             log(chltarg), crit1,  col = "red")
    segments(xlo, crit2,
             log(chltarg), crit2,  col = "red")
    segments(log(chltarg), crit2,
             log(chltarg), ylo, col = "red")
    segments(xlo, crit1,
             xlo, crit2, col = "red", lwd = 3)

    cat("TP criterion:", round(exp(crit2)), "\n")


}

## runmod variable set to T to run simulation and set to F to
##  run post processing.
fitout <- ntumodel(dat.merge.all, runmod = T)
## post processing
#varout.nou <- extract(fitout, pars = c("muk", "mub", "b",  "d2","d1s",
#                              "muu", "muu_mn", "u", "mud", "sigd"))
#umean <- apply(varout.nou$u, 2, mean)

#ntumodel(dat.merge.all, varout = varout.p.nat, varout.mo = varout.mo.d1L.d2T, runmod = F)
