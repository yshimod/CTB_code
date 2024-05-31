utility <- function(idx, price, t0, k, income, rho, delta, beta) {
    discounting <- beta ^ t0 * delta ^ k
    c_sooner <- (idx / 100) * income / price
    c_later <- income - price * c_sooner

    (c_sooner ^ rho + discounting * c_later ^ rho) / rho
}


neglog_probability <- function(idx, price, t0, k, income, rho, delta, beta) {
    lambda <- 1
    numerator <- exp( lambda * utility(idx, price, t0, k, income, rho, delta, beta) )
    denominator <- sum(
        sapply(
            0:100,
            function(jdx) {
                exp( lambda * utility(jdx, price, t0, k, income, rho, delta, beta) )
            }
        )
    )

    -log( numerator / denominator )
}


lnsigma_func <- function(theta) {
    4 * tanh(theta) + 1.5
}


rho_func <- function(theta) {
    1 - exp(-lnsigma_func(theta))
}


objfunc <- function(params, choices, price, t0, k, income) {
    deltai <- params[1]
    betai <- params[2]
    rhoi <- rho_func(params[3])
    sum(
        sapply(
            1:length(choices),
            function(i) {
                neglog_probability(choices[i], price[i], t0[i], k[i], income[i], rhoi, deltai, betai)
            }
        )
    )
}



set.seed(1)

df <- read.csv("base_perturbed.csv")

df$t0 <- ifelse(df$time == 0, 1, 0)

df$cons1 <- df$cons_sooner * df$price_sooner / df$income
df$choice <- round( 100 * df$cons1 )

nprob <- 42
numofindiv <- length(unique(df$ID))

rescolumns <- list(NULL, c(
    "labnumber",  # 1
    "delta",  # 2
    "rho",  # 3
    "beta",  # 4
    "sd",  # 5
    "convergei",  # 6
    "alphai",  # 7
    "betai",  # 8
    "deltai",  # 9
    "se_alpha",  # 10
    "se_beta",  # 11
    "se_delta",  # 12
    "ratei",  # 13
    "lnsigmai",  # 14
    "cparami",  # 15
    "se_rate",  # 16
    "se_lnsigma",  # 17
    "se_cparam",  # 18
    "ll",  # 19
    "ic"  # 20
))
allres <- matrix(NA, nrow = numofindiv, ncol = length(rescolumns[[2]]))
dimnames(allres) <- rescolumns

for (i in 1:numofindiv) {
    allres[i, 1] <- df$ID[(i - 1) * nprob + 1]
    allres[i, 2] <- df$delta[(i - 1) * nprob + 1]
    allres[i, 3] <- df$rho[(i - 1) * nprob + 1]
    allres[i, 4] <- df$beta[(i - 1) * nprob + 1]
    allres[i, 5] <- df$sd[(i - 1) * nprob + 1]

    temprows <- ((i - 1) * nprob + 1):(i * nprob)

    mleres <- try(
        optim(
            par = c(1., 1., -.07),
            fn = objfunc,
            choices = df$choice[temprows],
            price = df$price_sooner[temprows],
            t0 = df$t0[temprows],
            k = df$lag[temprows],
            income = df$income[temprows],
            control = list(maxit = 10000, abstol = 1e-5),
            hessian = TRUE,
            method = "Nelder-Mead"
        )
    )
    if (class(mleres) == "try-error") {
        mleres <- list(
            convergence = 0,
            par = c(NA, NA, NA),
            value = NA,
            counts = list(NA),
            hessian = matrix(NA, nrow = 3, ncol = 3)
        )
    }

    allres[i, 6] <- ifelse(mleres$convergence == 0, 1, 0)

    allres[i, 8] <- mleres$par[2]  # betai
    allres[i, 9] <- mleres$par[1]  # deltai
    allres[i, 15] <- mleres$par[3]  # cparami

    allres[i, 19] <- -mleres$value  # ll
    allres[i, 20] <- mleres$counts[[1]]  # ic

    tempse <- rep(NA, 3)
    try(tempse <- sqrt(abs(diag(solve(mleres$hessian)))))

    allres[i, 11] <- tempse[2]  # se_beta
    allres[i, 12] <- tempse[1]  # se_delta
    allres[i, 18] <- tempse[3]  # se_cparam

    allres[i, 14] <- lnsigma_func(mleres$par[3])  # lnsigmai
    allres[i, 17] <- 4 * ( 1 - tanh( mleres$par[3] ) ^ 2 ) * abs( tempse[3] )  # se_lnsigma

    allres[i, 7] <- rho_func(mleres$par[3])  # alphai
    allres[i, 10] <- 4 * ( 1 - tanh( mleres$par[3] ) ^ 2 ) * exp( - 4 * tanh( mleres$par[3] ) - 1.5 ) * abs( tempse[3] )  # se_alpha

    allres[i, 13] <- max(1e-5, mleres$par[1] ^ (-365)) - 1  # ratei
    allres[i, 16] <- abs( (-365) * max(1e-5, mleres$par[1] ^ (-366)) ) * abs( tempse[1] )  # se_rate
}

save(
    allres,
    file = "base_estimated_mlogit.Rdata"
)
