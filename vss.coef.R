## look at relationship between vss exponent and land use

vss.coef <- function(varout, df1) {
    d <- apply(varout$d, 2, mean)

    df1$lake <- factor(df1$lake)
    df1$lakenum <- as.numeric(df1$lake)
    dftemp <- unique(data.frame(df1[, c("lakenum", "logit_crop")]))
    dftemp <- dftemp[order(dftemp$lakenum),]
    plot(plogis(dftemp$logit_crop), exp(d))
}

vss.coef(varout.vss, moi3.all)
