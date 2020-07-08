library(MASS)

MCMCmetrop <- function(init.state, target, niter, burnin) {
	# Dimension of state space
	dim <- length(init.state)
	# Current state of the Markov chain
	theta <- init.state
	# Generate matrix containing the samples. Initialize first sample with the starting value
    samples <- matrix(0, nrow=niter+1, ncol=dim, dimnames=list(NULL, names(init.state)))
	samples[1, ] <- init.state
    # Generate uniform random numbers in advance, to save computation.
    u <- runif(niter)
    # Proposal is a multivariate standard normal distribution. Generate samples and
    # later on use linearity property of Gaussian distribution
    normal_shift <- mvrnorm(n=niter, mu=integer(dim), Sigma=diag(dim))
    for (i in 1:niter) {
        # Sample a candidatf
        candidate <- samples[i, ] + normal_shift[i, ]
        # calculate log target of candidate and store it in case it gets accepted
		r <- target(candidate) / target(samples[i, ])
        if (u[i] < r) theta <- candidate
		samples[i+1, ] <- theta
    }
    return(samples[-(1:burnin), ])
}

# Prior distribution for our parameter
prior <- function(theta) {
	beta <- theta["beta"]
	gamma <- theta["gamma"]
	dnorm(beta, mean=1) * dnorm(gamma, mean=1)
}

samples <- MCMCmetrop(c(beta=0, gamma=0), prior, 1000, 100)
layout(c(1,2,3), heights=c(3,1,1))
plot(samples[,1], samples[,2], xlab="beta", ylab="gamma")
plot(samples[,1], type="l", xlab="iteration", ylab="beta")
plot(samples[,2], type="l", xlab="iteration", ylab="gamma")
