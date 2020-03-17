library(MASS)

generate_eps <- function(params, N) rnorm(N, 0, params$noise_sd)

simulated_data <- function(x, y, params, class) {
    data <- list(x = x, y = y, params = params)
    class(data) <- class
    data
}

generate_data <- function(params, N) UseMethod("generate_data")

generate_data.univariate_params <- function(params, N) {
    x <- rnorm(N, 0, 1)
    y <- params$beta * x + generate_eps(params, N)
    simulated_data(x, y, params, "univariate_data")
}

generate_data.multivariate_params <- function(params, N) {
    x <- mvrnorm(N, params$mean, diag(params$sd))
    y <- x %*% params$beta + generate_eps(params, N)
    simulated_data(x, y, params, "multivariate_data")
}


calc_coef <- function(data) UseMethod("calc_coef")
calc_coef.univariate_data <- function(d) cov(d$x,d$y) / var(d$x)
calc_coef.multivariate_data <- function(d) {
    x <- d$x
    y <- d$y
    solve(t(x) %*% x, t(x) %*% y)
}


calc_se <- function(data, coef) UseMethod("calc_se")

calc_se.univariate_data <- function(data, coef) {
    n <- length(data$y)
    eps <- data$y - data$x * coef
    e_sd <- mean(eps^2)
    se <- sqrt(e_sd / n*var(data$x))
    se
}
calc_se.multivariate_data <- function(data, coef) {
    n <- length(data$y)
    eps <- data$y - data$x %*% coef
    e_sd <- mean(eps^2)
    se <- sqrt(e_sd / diag(n*var(data$x)))
    se
}

run_regression <- function(data) {
    coef <- calc_coef(data)
    se <- calc_se(data, coef)
    list(coef=coef, se=se)
}

eval_model <- function(params, coef, se) UseMethod("eval_model")
eval_model.univariate_params <- function(params, coef, se, conf=1.96) {
    beta <- params$beta
    up <- coef + se*conf
    down <- coef - se*conf
    (beta > down) & (beta < up)
}
eval_model.multivariate_params <- function(params, coef, se, conf=1.96) {
    beta <- params$beta
    up <- coef + se*conf
    down <- coef - se*conf
    all(beta > down) & all(beta < up)
}

simulate <- function(simulate_params, N) {
    data <- generate_data(simulate_params, N)
    m <- run_regression(data)
    eval_model(simulate_params, m$coef, m$se)
}

avg_simulations <- function(M, N, simulate_params) {
    inside <- sapply(1:M, function(x) {
        simulate(simulate_params, N)
    })
    sum(inside) / M
}

multivariate_params<- function(beta, means, sds, noise_sd) {
    params <- list(beta=beta, means=means, sd=sds, noise_sd=noise_sd)
    class(params) <- "multivariate_params"
    params
}

univariate_params<- function(beta, noise_sd) {
    params <- list(beta=beta, noise_sd=noise_sd)
    class(params) <- "univariate_params"
    params
}

# NOTE:
# This code, the "eval_model" is testing whether or not the WHOLE model is correct,
# which is NOT actually what the SE are supposed to report, but shows the danger
# of thinking of them in this way (multiple testing!)

# avg_simulations(10, 10, multivariate_params(c(1,1,1), c(0,0,0), c(1,1,1), .5))
# avg_simulations(100, 10, multivariate_params(c(1,1,1), c(0,0,0), c(1,1,1), .5))
# avg_simulations(100, 100, univariate_params(1, .5))
