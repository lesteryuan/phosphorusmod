## plot temperature
## 2.5.2021

tempplot <- function(varout.test = NULL, df1) {

    ## location of radiation estimate
    lat0 <- 39.61
    lon0 <- -93.38
    H <- 5 # assume mixed layer is about 3m for now

    ## get time of each radiation prediction
    datestr <- paste(radiation$Year,"-", radiation$Month, "-", radiation$Day,
                     " ", radiation$Hour, ":", radiation$Minute, sep = "")
    radiation$date0 <- as.POSIXct(strptime(datestr, format = "%Y-%m-%d %H:%M",
                                 tz = "America/Chicago"))
    attr(radiation$date0, "tzone") <- "GMT"

    ## get sun angle at eaech time
    require(oce)
    ang0 <- sunAngle(radiation$date0, longitude= lon0, latitude = lat0)

    ## calculate GHI from DNI and DHI
    altitude <- ang0$altitude
    altitude[altitude <0] <- 0
    r <- sin(altitude*pi/180)
    radiation$GHI <- r*radiation$DNI #+ radiation$DHI  # trying just with dni
    print(summary(radiation$GHI))

    ## get posix data for moi3
    datestr <- paste("2004-", df1$month, "-", df1$day, " 12:00", sep = "")
    df1$date0 <- as.POSIXct(strptime(datestr, format =  "%Y-%m-%d %H:%M",
                                     tz = "America/Chicago"))
    attr(df1$date0, "tzone") <- "GMT"
    print(df1$date0[1:10])

    ## one week before sample
    starttime <- df1$date0 - 7*24*3600
    print(starttime[1:10])

    ## calculate average radiation the week before each sample
    df1$rad.av <- rep(NA, times = nrow(df1))
    for (i in 1:nrow(df1)) {
        incvec <- radiation$date0 >= starttime[i] & radiation$date0 <= df1$date0[i]
        incvec[is.na(incvec)] <- F
        df1$rad.av[i] <- mean(radiation$GHI[incvec], na.rm = T)
    }

    ## calculate attenuation, ka, from secchi
    ## try ka*SD from paddial and thomaz
    ka <- 2.26/df1$secchi
    df1$rad.av.k <- df1$rad.av/(ka*H)*(1-exp(-ka*H))

#    plot(df1$date0, df1$rad.av,
#         ylim = range(c(df1$rad.av, df1$rad.av.k), na.rm = T))
#    points(df1$date0, df1$rad.av.k, pch = 16)

    ## paper gives 0.4 mol quanta/m2-d is approximately 1 W/m2
    f <- 0.4
    df1$chl.c.pred <- 0.003 + 0.0154*exp(0.050*df1$surftemp)*
        exp(-0.059*f*df1$rad.av.k)

    df1$month2 <- ceiling(df1$month*0.5)

    chl.c.month <- tapply(df1$chl.c.pred, df1$month2, mean, na.rm = T)
    dev.new()
    par(mar = c(4,4,1,1), mfrow = c(2,1))
    r0 <- c(min(chl.c.month) - 0.008, max(chl.c.month) + 0.012)
    b <- apply(varout.test$b, 2, mean)
    k <- mean(varout.test$k)
    load("mn.val.vss.rda")
    ns <- length(b)
    predout <- matrix(NA, ncol = 3, nrow = ns)
    for (i in 1:ns) {
        y <- 1/(0.50*1000*exp(varout.test$b[,i])/(mn.val["chl"]^varout.test$k)*mn.val["vss"]*30^(varout.test$k-1))
#        y <- 1/(0.50*1000*exp(varout.test$b[,i])*30^(varout.test$k-1))
        predout[i,] <- quantile(y, prob = c(0.025, 0.5, 0.975))
    }
    plot(as.numeric(names(chl.c.month)), as.vector(chl.c.month),
         ylim = range(c(r0, predout)))
    segments(1, chl.c.month[1] - 0.008, 1, chl.c.month[1] + 0.008)
    segments(4, chl.c.month[4] - 0.012, 4, chl.c.month[4] + 0.012)
    points(1:ns, predout[,2],pch = 16, col = "blue")
    segments(1:ns, predout[,1], 1:ns, predout[,3], col = "blue")

    stop()

    plot(1:ns, predout[,2], ylim = range(predout), type = "n")


}

tempplot(varout.test = varout.v1, df1 = moi3.all)
