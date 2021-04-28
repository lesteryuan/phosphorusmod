## explore size fraction effects on tss-chl relationship
## 3.5.2021

## notes: seasonal difference in size distribution doesn't seem
## very large. May not be worth building a model

## nchl is chl after passing through a "medium filter"
## uchl is chl after passing through a fine filter
tsssize <- function(df1) {

    mn0 <- tapply(df1$nchl/df1$chl, df1$month, mean)
    plot(1:12, mn0)

}

tsssize(moi3.all)
