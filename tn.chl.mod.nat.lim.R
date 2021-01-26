## 11.6.2019: Production run for TN model
## 12.17.2019: Cleaned and commented

tn.model <- function(df1, varout = NULL, varout.mo = NULL, runmod = F) {
    require(rstan)
    nchains <- 3    # number of chains

    ecosel <- 65              # pick ecoregion
    chltarg <- 10             # pick chl target
    doc0 <- 5                 # DOC

    load("tpchldat.rda")
    #source("logtick.exp.R")

    ## omit records that are missing data
    incvec <- ! is.na(df1$ntl.result) & ! is.na(df1$no3no2.result) &
        ! is.na(df1$chl) & ! is.na(df1$doc.result)
    cat("N omitted due to missing records:", sum(!incvec), "\n")
    df1 <- df1[incvec,]
    print(nrow(df1))

    ## drop chl = 0
    incvec <- df1$chl > 0
    df1 <- df1[incvec,]

    ## drop HI
    incvec <- df1$st.nla2012 == "HI"
    df1 <- df1[!incvec,]
    print(sum(incvec))

    ## compute TN-DIN and drop values that are <= 0
    df1$tkn <- df1$ntl.result - df1$no3no2.result
    incvec <- df1$tkn <= 0
    cat("TKN <= 0:", sum(incvec), "\n")
    df1 <- df1[!incvec,]

    ## select index sites
    incvec <- df1$sample.type == "MICX"
    df1 <- df1[incvec,]

    ## reset state and L3 ecoregion factors
    df1$state <- factor(df1$st.nla2012)
    df1$statenum <- as.numeric(df1$state)
    df1$us.l3code <- factor(df1$us.l3code)
    df1$econum <- as.numeric(df1$us.l3code)

    ## drop samples with chl > 100
    incvec <- df1$chl < 100
    df1 <- df1[incvec,]
    incvec <- df1$chl > 1
    df1 <- df1[incvec,]

    ## center chl and doc
    ## use same chl scale as TP model
    df1$chl.sc <- df1$chl/chlsc

    docsc <- exp(mean(log(df1$doc.result)))
    df1$doc.sc <- df1$doc.result/docsc

    ## center tn and nox
    tnsc <- exp(mean(log(df1$ntl.result)))
    df1$tn.sc <- df1$ntl.result/tnsc
    df1$nox.sc <- df1$no3no2.result/tnsc

    ## make seas factor
    xcut <- seq(min(df1$yday.x), by = 30, length = 6)
    xcut[6] <- max(df1$yday.x)
    seasfac <- cut(df1$yday.x, xcut, include.lowest = T)
    df1$seasnum <- as.numeric(seasfac)

    tnchldat <- df1[, c("chl", "chl.sc", "tkn", "doc.sc", "doc.result",
                        "statenum", "state", "ntl.result", "tn.sc",
                        "no3no2.result", "nox.sc", "index.lat.dd",
                        "index.lon.dd", "econum", "us.l3code")]
    save(tnchldat, docsc,tnsc, file = "tnchldat.rda")

    datstan <- list(n = nrow(df1),
                    neco = max(df1$econum),econum = df1$econum,
                    tn = log(df1$tn.sc),
                    nox = df1$nox.sc,
                    doc = df1$doc.sc,
                    chl = df1$chl.sc,
                    nseas = max(df1$seasnum),
                    seasnum = df1$seasnum,
                    dp = log(81*docsc/tnsc))

    print(str(datstan))


    modstan <- '
        data {
            int n;                // number of samples
            int neco;            // number of ecoregions
            int econum[n];       // ecoregion of each sample
            vector[n] tn;        // scaled TN
            vector[n] chl;       // scaled chl
            vector[n] nox;       // scaled NOx
            vector[n] doc;       // scaled DOC
            int nseas;
            int seasnum[n];
            real dp;
        }
        parameters {
            real muk;                // mean value of exponent on chl
            real mud[2];              // mean value of model coefficients
            real<lower = 0> sigd[2]; // SD of model coefficients among ecoregions
            vector[nseas] etad1;
            vector[neco] etad2;

            real<lower = 0> sigtn;  // measurement error of tn
            real muu;
            vector[n] etau;
            real<lower = 0> sigu;

        }
        transformed parameters {
            vector[nseas] d1;
            vector[neco] d2;
            vector[n] u;
            d1 = mud[1] + sigd[1]*etad1;
            d2 = mud[2] + sigd[2]*etad2;
            u = muu + etau*sigu;
        }
        model {
            vector[n] tnmean;
            muk ~ normal(1,1);    // small amount of information for k
            mud[1] ~ normal(0,4);
            mud[2] ~ normal(0,4);
            sigd ~ cauchy(0,4);

            etad1 ~ normal(0,1);
            etad2 ~ normal(0,1);

            muu ~ normal(0,3);
            etau ~ normal(0,1);
            sigu ~ cauchy(0,3);

            sigtn ~ normal(0.1, 0.002);
           // Eqn 34
           for (i in 1:n) tnmean[i] = nox[i] +
                                 exp(d1[seasnum[i]])*chl[i]^muk +
                                 exp(d2[econum[i]])*doc[i] + exp(u[i]);
            tn ~ student_t(4,log(tnmean), sigtn);
        }
    '
    if (runmod) {
        rstan_options(auto_write = TRUE)
        options(mc.cores = nchains)
        fit <- stan(model_code = modstan,
                    data = datstan, iter = 900, chains = nchains,
                    warmup = 300, thin = 2)
        return(fit)
    }



    grey.t <- adjustcolor("grey39", alpha.f = 0.5)

    mud<- apply(varout$mud, 2, mean)
    muk <- mean(varout$muk)

    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(2,2))
    plot(log(moi3.all$doc), log(moi3.all$don) + log(1000))
    mod <- lm(log(don) -log(doc) ~ 1, data = moi3.all)
    print(summary(mod))

    print(confint(mod))
    abline(coef(mod) + log(1000), 1, col = "blue")

    abline(mud[2] - log(docsc) + log(tnsc), 1)
    plot(log(df1$doc.sc), log(df1$tn.sc - df1$nox.sc))
    abline(mud[2], 1)

    nsamp <- nrow(varout$mud)
    y <- rnorm(nsamp, mean = varout$mud[,2], sd = varout$sigd[,2]) -
        log(docsc) + log(tnsc) - log(1000)
    print(quantile(y, prob = c(0.025, 0.5, 0.975)))
    stop()
    u <- apply(varout$u, 2, mean)
    plot(log(df1$chl), u)
    stop()


    d1.mo <- apply(varout.mo$d1, 2, mean)
    k.mo <- apply(varout.mo$k, 2, mean)
    load("mn.val.tnmo.rda")

    xnew <- seq(min(log(1)), max = log(100), length = 50)
    predout1 <- matrix(NA, ncol = 3, nrow = length(xnew))
    predout2 <- matrix(NA, ncol = 3, nrow = length(xnew))
    ns1 <- nrow(varout$mud)

    for (i in 1:length(xnew)) {
        y <- rnorm(ns1, mean = varout$mud[,1], sd =varout$sigd[,1]) +
            varout$muk*(xnew[i]-log(chlsc)) + log(tnsc)
        predout1[i,] <- quantile(y, prob = c(0.05, 0.5, 0.95))
        y2 <- varout.mo$d1[,4] + varout.mo$k[,1]*(xnew[i] - log(mn.val["chl"])) +
            log(mn.val["tn"]) + log(1000)
        predout2[i,] <- quantile(y2, prob = c(0.05, 0.5, 0.95))
    }

    dev.new()
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0))
    plot(log(df1$chl.sc*chlsc), log((df1$tn.sc - df1$nox.sc)*tnsc),
         axes = F, xlab = expression(Chl~(mu*g/L)),
         ylab = expression(TN~(mu*g/L)),pch = 21, col = "grey39",
         bg = "white")
    points(log(moi3.all$chl), log(moi3.all$tn - moi3.all$dtn) + log(1000),
           pch = 16, cex = 0.6)
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    polygon(c(xnew, rev(xnew)), c(predout1[,1], rev(predout1[,3])),
            col = grey.t, border= NA)
#    abline(mud[1] + log(tnsc) - muk*log(chlsc), muk)
#    abline(d1.mo[4] + log(mn.val["tn"]) - k.mo[1]*log(mn.val["chl"]) + log(1000),
#           k.mo[1],lty = "dashed")
    lines(xnew, predout2[,1])
    lines(xnew, predout2[,3])
    stop()

    draw <- d
    draw[,1] <- d[,1] - log(docsc) + log(tnsc)
    draw[,2] <- d[,2] - k*log(chlsc) + log(tnsc)

    print("*** ecoregion specific coefs ***")
    print(range(exp(draw[,1])))
    print(range(exp(draw[,2])))

    # merge in coefficient to main file
    dftemp <- data.frame(econum = 1:length(d),draw, d, k)
    names(dftemp) <- c("econum", "draw1", "draw2", "d1", "d2", "k")
    df1 <- merge(df1, dftemp, by = "econum")

    # save coefficients for mapping
    dftemp <- merge(dftemp, unique.data.frame(df1[, c("econum", "us.l3code")]),
                    by = "econum")

    coefs <- dftemp
    save(coefs, file= "coefs.rda")

    # print numerical summaries of parameters
    print("d1")
    print(exp(quantile(varout$mud[,1] - log(docsc) + log(tnsc),
                       prob = c(0.05, 0.5, 0.95))))
    print("d2")
    print(exp(quantile(varout$mud[,2] - varout$muk*log(chlsc) + log(tnsc),
                       prob = c(0.05, 0.5, 0.95))))


    print("k")
    print(quantile(varout$muk, prob = c(0.05, 0.5, 0.95)))

    ## calculate predicted TN
    df1$predout <- log(exp(df1$d1)*df1$doc.sc +
                           exp(df1$d2)*df1$chl.sc^df1$k) + log(tnsc)

    ## PLOT: Full data set pred.v.obs
    png(width = 3, height = 3, pointsize = 7, units = "in", res = 600,
        file = "predobs.png")
    par(mar = c(4,4,1,1), mgp = c(2.3,1,0))
    plot(df1$predout, log(df1$tkn), xlab = expression(Predicted~TN-DIN~(mu*g/L)),
         ylab = expression(Observed~TN-DIN~(mu*g/L)),
         pch = 21, col = "grey39", bg = "white", axes = F)
    logtick.exp(0.001, 10,  c(1,2), c(F,F))
    abline(0,1)
    dev.off()
    rms <- sqrt(sum((df1$predout - log(df1$tkn))^2)/
                    nrow(df1))
    cat("RMS:", rms, "\n")

    # PLOT: Pred.v.obs by state
    pdf(file = "plots.pdf", width = 9, height = 6, pointsize = 10)
    par(mar = c(4,4,3,2), mfrow = c(2,3), mgp = c(2.3,1,0))
    for (i in 1:max(df1$econum)) {
        incvec <- df1$econum == i
##        plot(df1$predout[incvec], log(df1$tkn)[incvec], axes = F, xlab = "Predicted TN-DIN",
##             ylab = "Obs TN-DIN", pch = 21, col = "grey39", bg = "white",
##             main = i)
        plot(log(df1$chl.sc[incvec]), log(df1$tn.sc[incvec]-df1$nox.sc[incvec]),
             main = i, axes = F)
        abline(d[i,2], k)

        logtick.exp(0.001, 10, c(1,2), c(F,F))
#        abline(0,1)
    }
    dev.off()

    # PLOT: Chl v. TN, DOC v. TN  (Fig 31)
#    png(width = 6, height = 2.5, pointsize =8, units = "in", res = 600,
#        file = "tnchl.png")
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(1,2), mgp = c(2.3,1,0))

    plot(log(df1$chl), log(df1$ntl.result),
         col = "grey", pch = 21,
         bg = "white", axes = F, xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = expression(TN-DIN~(mu*g/L)))
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    x <- seq(min(log(df1$chl.sc)), max(log(df1$chl.sc)), length = 50)
    xraw <- x + log(chlsc)
    predout <- matrix(NA, ncol = 3, nrow = length(x))
    predout.mo <-matrix(NA, ncol = 3, nrow = length(x))

    load("mn.val.tnmo.rda")
    load("vartn.0.rda")

    for (i in 1:length(x)) {
        y <- varout$mud[,2] + varout$muk*x[i] + log(tnsc)
        predout[i,] <- quantile(y, prob = c(0.05, 0.5, 0.95))
        y2 <- log(mn.val["tn"]) + vartn.0$mud[,1] -
            vartn.0$k[,1]*log(mn.val["chl"]) + vartn.0$k[,1]*xraw[i] +
                log(1000) + log(2)
        predout.mo[i,] <- quantile(y2, prob = c(0.05, 0.5, 0.95))
    }
    polygon(c(xraw, rev(xraw)), c(predout[,1], rev(predout[,3])),
            col = grey.t, border= NA)
    lines(xraw, predout[,2])
    lines(xraw, predout.mo[,1], lty = "dashed")
    lines(xraw, predout.mo[,3], lty = "dashed")
    stop()

    plot(log(df1$doc.result), log(df1$ntl.result), col = "grey", pch = 21,
         bg = "white", axes = F, xlab = "DOC (mg/L)",
         ylab = expression(TN-DIN~(mu*g/L)))
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    x <- seq(min(log(df1$doc.sc)), max(log(df1$doc.sc)), length = 50)
    xraw <- x + log(docsc)
    for (i in 1:length(x)) {
        y <- varout$mud[,1] + x[i] + log(tnsc)
        predout[i,] <- quantile(y, prob = c(0.05, 0.5, 0.95))
    }
    polygon(c(xraw, rev(xraw)), c(predout[,1], rev(predout[,3])),
            col = grey.t, border= NA)
    lines(xraw, predout[,2])
    dev.off()

    ieco <- which(levels(df1$us.l3code) == ecosel)

#    doc.mn <- tapply(df1$doc.sc, df1$econum, function(x) mean(log(x)))
#    print(exp(doc.mn[ieco]+log(docsc)))
    doc0.sc <- doc0/docsc

    grey.t <- adjustcolor("grey39", alpha.f = 0.5)
                                        # calculate parameters

    incvec <- df1$econum == ieco
    x <- seq(min(log(df1$chl.sc)), max(log(df1$chl.sc)), length = 50)
    xraw <- x + log(chlsc)
    predout <- matrix(NA, ncol = 3, nrow = length(x))
    predout2 <- matrix(NA, ncol = 3, nrow = length(x))

    for (i in 1:length(x)) {
        y <- varout$d[,ieco,2]  + varout$muk*x[i] + log(tnsc)
        predout[i,] <- quantile(y, prob = c(0.10, 0.5, 0.90))
        y2 <- exp(varout$d[,ieco,1])*doc0.sc +
            exp(varout$d[,ieco,2])*exp(x[i])^varout$muk
        predout2[i,] <- quantile(log(y2) + log(tnsc), prob = c(0.10, 0.5, 0.9))
    }

    cat("Median DOC for ecoregion:", median(df1$doc.result[incvec]), "\n")

    ## Figure 32
    dev.new()
    par(mgp = c(2.3,1,0), bty = "l", mar = c(4,4,1,1))
    plot(log(df1$chl), log(df1$tn.sc*tnsc - df1$nox.sc*tnsc), type = "p",
         axes = F, xlab = expression(Chl~italic(a)~(mu*g/L)),
         ylab = expression(TN-DIN~(mu*g/L)), col = "grey80",
         pch = 21, bg = "white")
    logtick.exp(0.001, 10, c(1,2), c(F,F))
    points(log(df1$chl)[incvec], log(df1$tn.sc*tnsc - df1$nox.sc*tnsc)[incvec],
           pch = 21, col = "black", bg = "grey")

    polygon(c(xraw, rev(xraw)), c(predout[,1], rev(predout[,3])),
            col  = grey.t, border = NA)
    lines(xraw, predout[,2])
    polygon(c(xraw, rev(xraw)), c(predout2[,1], rev(predout2[,3])),
            col  = grey.t, border = NA)
    lines(xraw, predout2[,2], lty = "dashed")

    crit1 <- approx(xraw, predout[,1], log(chltarg))$y
    crit2 <- approx(xraw, predout2[,1], log(chltarg))$y

    incvec <- df1$tn.sc > df1$nox.sc
    ylo <- min(log(df1$tn.sc[incvec]*tnsc - df1$nox.sc[incvec]*tnsc)) -
        0.04*diff(range(log(df1$tn.sc[incvec]*tnsc - df1$nox.sc[incvec]*tnsc)))
    xlo <- min(log(df1$chl)) - 0.04*diff(range(log(df1$chl)))
    segments(xlo, crit1,
             log(chltarg), crit1,  col = "red")
    segments(xlo, crit2,
             log(chltarg), crit2,  col = "red")
    segments(log(chltarg), crit2,
             log(chltarg), ylo, col = "red")
    segments(xlo, crit1,
             xlo, crit2, col = "red", lwd = 3)

    cat("TN Criterion:", round(exp(crit2)), "\n")

    stop()

}

## save extracted variables to varout to post-process
#fitout <- tn.model(dat.merge.all, runmod = T)

varout.n <- extract(fitout, pars = c("muk", "mud", "d1", "sigd", "u"))
tn.model(dat.merge.all, varout = varout.n, varout.mo = varout.mon.d1T.d2Lv,
         runmod = F)
