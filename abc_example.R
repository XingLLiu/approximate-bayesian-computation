discrepancy <- function(x, y){
    discrepancy <- sum(abs(x - y)) / length(x)
    return(discrepancy)
}

abc_sample <- function(N, epsilon, y, alpha=1, beta=1){

    samples <- rep(0, N)

    for (i in 1:N) {
        rho <- epsilon + 1
        accept_prob <- rho <= epsilon

        while (runif(n = 1) > accept_prob){
            theta <- rbeta(n = 1, shape1 = alpha, shape2 = beta)
            x <- rbinom(n, prob = theta, size = 1)
            rho <- discrepancy(x, y)
            accept_prob <- rho <= epsilon
        }

        samples[i] <- theta
    }
    
    return(samples)
}


kernel <- function(x, y, epsilon){
    return(discrepancy(x, y) <= epsilon)
}


set.seed(1)
p_0 <- 0.7
sample_size <- 3
n <- 10
y <- rbinom(n, prob = p_0, size = 1)

samples <- abc_sample(N = 500, epsilon = 0, y = y, alpha = 2, beta = 1)

