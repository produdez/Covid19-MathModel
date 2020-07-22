# Initial state for the Markov chain
theta <- c(beta=0, gamma=0)
func <- function(theta) prod(dnorm(theta, sd=100))
niter <- 10000
burnin <- 10
sd <- 2.4 / sqrt(2) * 100

ptheta <- func(theta)
# Generate matrix containing the samples. Initialize first sample with the starting value
sample <- matrix(nrow=niter, ncol=2, dimnames=list(NULL, names(start)))
# Generate uniform random numbers in advance, to save computation.
u <- runif(niter)
# Proposal is a multivariate normal distribution. Generate samples and
# later on use linearity property of Gaussian distribution
normal.shift <- matrix(rnorm(niter*2,sd=sd), nrow=niter, ncol=2)
for (i in seq_len(niter)) {
	# Sample a candidate
	candidate <- theta + normal.shift[i, ]
	# Calculate func of candidate and store it in case it gets accepted
	r <- func(candidate) / ptheta
	if (u[i] < r) {
		theta <- candidate
		ptheta <- func(theta)
	}
	if (i > burnin) sample[i, ] <- theta
}

accept <- 1 - sum(apply(diff(samples) == 0, 1, FUN=all)) / (nrow(samples) - 1)
cat("Acception rate is", accept, "\n")
