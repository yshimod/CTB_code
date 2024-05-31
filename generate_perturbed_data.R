library(truncnorm)

add_noise <- function (odata, sdval, sd_idx){
    repnum <- 10
    n <- length(odata$ID)

    for (i in 0:(repnum-1)) {
        original_relcons <- ( odata$cons_sooner * odata$price_sooner ) / odata$income
        relcons <- rtruncnorm(
            length(original_relcons),
            a = 0,
            b = 1,
            mean = original_relcons,
            sd = sdval
        )

        tempdf <- odata
        tempdf$ID <- odata$ID * 10000 + sd_idx * 1000 + i
        tempdf$sd <- sdval
        tempdf$cons_sooner <- odata$income * relcons / odata$price_sooner
        tempdf$cons_later <- odata$income * (1 - relcons) / odata$price_later
        tempdf$noise <- relcons - original_relcons

        if ( i == 0 ) {
            ndata <- tempdf
        }
        else {
            ndata <- merge(ndata, tempdf, all=T)
        }
    }
    return(ndata)
}



originaldata <- read.csv("base_zeronoise.csv")

sd_list = c(0.01, 0.05, 0.10, 0.15, 0.20)

data <- data.frame()
for (sd_idx in 1:length(sd_list)) {
    set.seed(2222+sd_idx)
    data <- rbind(
        data,
        add_noise(originaldata, sd_list[sd_idx], sd_idx)
    )
}

write.csv(data, "base_perturbed.csv", row.names=FALSE)
