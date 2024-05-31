compute_quasi_homothetic_ces_demand_2d <- function(price, income, alpha, beta, neg_rho){
    x <- c(0, 0)

    if (abs(neg_rho) > 1e-9) {
        log_alpha_rho <- (1 / (1 + neg_rho)) * log((alpha[1] / alpha[2]))
        log_p_rho <- (-neg_rho / (1 + neg_rho)) * log((price[1] / price[2]))

        ratio_alpharho_prho <- log_alpha_rho - log_p_rho
        if ( ratio_alpharho_prho > 10 ) {
            ar_o_ar_pr <- 1
            pr_o_ar_pr <- 0
        }
        else if ( ratio_alpharho_prho < -10 ){
            ar_o_ar_pr <- 0
            pr_o_ar_pr <- 1
        }
        else {
            ar_o_ar_pr <- 1 / ( 1 + exp(-ratio_alpharho_prho) )
            pr_o_ar_pr <- 1 / ( 1 + exp(+ratio_alpharho_prho) )
        }

        if (
            (alpha[1]*price[2]) / (alpha[2]*price[1]) <= (beta[1] / (beta[2] + (income / price[2])))^(1+neg_rho)
        ) {
            x <- c(0, income/price[2])
        }
        else if (
            (alpha[1]*price[2]) / (alpha[2]*price[1]) >= ((beta[1] + (income / price[1])) / beta[2])^(1+neg_rho)
        ) {
            x <- c(income/price[1], 0)
        }
        else {
            x[1] <- ar_o_ar_pr * ((price[1] * beta[1] + price[2] * beta[2] + income) / price[1]) - beta[1]
            x[2] <- pr_o_ar_pr * ((price[1] * beta[1] + price[2] * beta[2] + income) / price[2]) - beta[2]
        }
    }
    else {
        if (
            (alpha[1] * price[2]) / (alpha[2] * price[1]) <= beta[1] / (beta[2] + (income / price[2]))
        ) {
            x <- c(0, income / price[2])
        }
        else if ((alpha[1] * price[2]) / (alpha[2] * price[1]) >= (beta[1] + (income / price[1])) / beta[2]) {
            x <- c(income / price[1], 0)
        }
        else {
            x[1] <- (alpha[1] / price[1]) * (price[1] * beta[1] + price[2] * beta[2] + income) - beta[1]
            x[2] <- (alpha[2] / price[2]) * (price[1] * beta[1] + price[2] * beta[2] + income) - beta[2]
        }
    }

    return(x)
}

generate_discount_present_ces_data_2d <- function(prices, incomes, times, lags, delta, beta, neg_rho, omega){
    num_goods <- 2
    num_data <- dim(prices)[1]

    consumptions <- matrix(0, num_data, num_goods)

    for(i in 1:num_data){
        if (times[i] == 0) {
            alpha1 <- 1 / (beta * delta^lags[i] + 1)
        }
        else {
            alpha1 <- 1 / (delta^lags[i] + 1)
        }

        alpha <- c(alpha1, 1 - alpha1)
        consumptions[i,] <- compute_quasi_homothetic_ces_demand_2d(
            price = c(prices[i,]),
            income = incomes[i],
            alpha = alpha,
            beta = omega,
            neg_rho = neg_rho
        )
    }

    return(consumptions)
}

generate_zeronoise_data <- function (delta, rho, beta) {
    subjs <- as.matrix(expand.grid(delta, rho, beta))
    num_sim <- length(subjs[,1])

    colList <- c(
        "ID",
        "delta",
        "rho",
        "beta",
        "sd",
        "time",
        "lag",
        "price_sooner",
        "price_later",
        "income",
        "cons_sooner",
        "cons_later"
    )
    df <- data.frame(matrix(rep(NA, length(colList)), nrow=1))[numeric(0), ]
    colnames(df) <- colList

    for (i in 1:num_sim) {
        sim_consumption <- generate_discount_present_ces_data_2d(
            prices = prices,
            incomes = incomes,
            times = times,
            lags = lags,
            delta = subjs[i, 1],
            beta = subjs[i, 3],
            neg_rho = -subjs[i, 2],
            omega = c(0, 0)
        )

        tempdf <- data.frame(
            ID = rep(i-1, length(sim_consumption[,1])),
            delta = rep(subjs[i, 1], length(sim_consumption[,1])),
            rho = rep(subjs[i, 2], length(sim_consumption[,1])),
            beta = rep(subjs[i, 3], length(sim_consumption[,1])),
            sd = rep(0, length(sim_consumption[,1])),
            time = times,
            lag = lags,
            price_sooner = prices[,1],
            price_later = prices[,2],
            income = incomes,
            cons_sooner = sim_consumption[,1],
            cons_later = sim_consumption[,2]
        )

        df <- merge(df, tempdf, all=T)
    }

    return(df)
}


# Problem set
times <- c(rep(0, 21), rep(1, 21))

lags <- rep(70, 42)

r_plus_1 <- rep(seq(0.6, 2, length=21), 2)
prices <- matrix(c(r_plus_1, rep(1, 42)), 42, 2)

incomes <- rep(20, 42)


# Ground-truth parameters
delta_list <- round(seq(0.9912, 1.0025, length=10) , 4)
rho_list <- -round(exp(seq(-2, 5, length=10))^(-1) - 1, 3)
beta_list <- seq(0.85, 1.12, length=10)


originaldata <- generate_zeronoise_data(delta_list, rho_list, beta_list)

write.csv(originaldata, "base_zeronoise.csv", row.names=FALSE)
