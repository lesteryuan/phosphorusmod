## 1.26.2021
## summarize RMS results

sumrms <- function() {
    fname <- c("d1L.d2L", "d1L.d2T", "d1T.d2L", "d1T.d2T")

    rmsout <- rep(NA, times = length(fname))
    names(rmsout) <- fname
    getrms <- function(x,y) sqrt(sum((x-y)^2, na.rm = T)/
                                     min(sum(!is.na(x)), sum(!is.na(y))))

    for (i in fname) {
        load(paste("matout.mo.", i, ".rda", sep = ""))
        df <- get(paste("matout.mo.", i, sep = ""))
        rmsout[i] <- getrms(df[,1], df[,2])
        print(sum(is.na(df$pred)))
    }

    print(rmsout)
}

sumrms()

