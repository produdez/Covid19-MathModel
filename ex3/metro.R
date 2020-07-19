library(MASS)

rw.metro <- function(start, func, nsamples, pcov.scale=1, ...) {
	opt <- optim(start, func, control=list(fnscale=-1), hessian=TRUE, ...)
	if (opt$convergence != 0) stop(opt$message)
	V <- -pcov.scale^2 * solve(opt$hessian)
	# Dimension of state space
	dim <- length(start)
	# Initial state for the Markov chain
	theta <- opt$par
	ptheta <- func(theta, ...)
	# Generate matrix containing the samples.
	samples <- matrix(nrow=nsamples, ncol=dim)
	colnames(samples) <- names(start)
	# Generate uniform random numbers in advance, to save computation.
	log.u <- log(runif(nsamples))
	# Proposal is a multivariate normal distribution. Generate samples and
	# later on use linearity property of Gaussian distribution
	normal.shift <- mvrnorm(n=nsamples, mu=rep(0, dim), Sigma=V)
	for (i in seq_len(nsamples)) {
		# Sample a candidate
		candidate <- theta + normal.shift[i, ]
		# Calculate func of candidate and store it in case it gets accepted
		log.r <- func(candidate, ...) - ptheta
		if (log.u[i] < log.r) {
			theta <- candidate
			ptheta <- func(theta, ...)
		}
		samples[i, ] <- theta
	}
	return(samples)
}

samples <- rw.metro(
	start = c(beta=0, gamma=0),
	func = function(theta) sum(dnorm(theta, sd=10, log=TRUE)),
	nsamples = 10000,
	pcov.scale = 1
)
