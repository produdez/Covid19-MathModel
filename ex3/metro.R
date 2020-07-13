rk4 <- function(y0, times, func, parms) {
	h <- times[2] - times[1]
	tbl <- matrix(nrow=length(times), ncol=1+length(y0))
	colnames(tbl) <- c("time", names(y0))
	tbl[, 1] <- times
	tbl[1, -1] <- y0
	for (i in 1:(nrow(tbl)-1)) {
		ti <- tbl[i, 1]
		yi <- tbl[i, -1]
		k1 <- func(ti, yi, parms)
		k2 <- func(ti + h/2, yi + h/2 * k1, parms)
		k3 <- func(ti + h/2, yi + h/2 * k2, parms)
		k4 <- func(ti + h, yi + h * k3, parms)
		tbl[i+1, -1] <- yi + h/6 * (k1 + 2 * k2 + 2 * k3 + k4)
	}
	return(tbl)
}

SIR.model <- function(t, y, parms) with(as.list(c(y, parms)), {
	dS.dt <- -beta / N * S * I
	dI.dt <- beta / N * I * S - gamma * I
	dR.dt <- gamma * I
	c(dS.dt, dI.dt, dR.dt)
})

MCMCmetrop <- function(init.state, log.target, nsamples, burnin) {
	niter <- nsamples + burnin - 1
	# Dimension of state space
	dim <- length(init.state)
	# Current state of the Markov chain
	theta <- init.state
	ptheta <- log.target(theta)
	# Generate matrix containing the samples. Initialize first sample with the starting value
    samples <- matrix(nrow=niter+1, ncol=dim, dimnames=list(NULL, names(init.state)))
	samples[1, ] <- init.state
    # Generate uniform random numbers in advance, to save computation.
    log.u <- log(runif(niter))
    # Proposal is a multivariate standard normal distribution. Generate samples and
    # later on use linearity property of Gaussian distribution
	normal.shift <- matrix(rnorm(niter*dim, sd=.1), nrow=niter, ncol=dim)
    for (i in 1:niter) {
        # Sample a candidate
        candidate <- samples[i, ] + normal.shift[i, ]
		# print(candidate)
        # Calculate log target of candidate and store it in case it gets accepted
		log.r <- log.target(candidate) - ptheta
        if (log.u[i] < log.r) {
			# if (any(candidate < 0)) stop(str(candidate), log.u[i], "\n", log.target(candidate), "\n", ptheta)
			theta <- candidate
			ptheta <- log.target(theta)
		}
		samples[i+1, ] <- theta
    }
    return(samples[-(1:burnin), ])
}

# Prior distribution for our parameter
log.prior <- function(theta) dlnorm(theta, sdlog=100, log=TRUE)

log.likelihood <- function(theta, data, pop) {
	# The parameters must be non-negative
	if (any(theta < 0)) return(-Inf)
	# Number of intervals in each day
	n <- 20

	result <- rk4(
		y0 = data[1, -1],
		times = seq(1, nrow(data), by=1/n),
		func = SIR.model,
		parms=c(theta, N=pop)
	)
	estimated <- result[seq(n+1, nrow(result), by=n), "I"]
	sum(dpois(data$I[-1], estimated, log=TRUE))
}

dataset <- "https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_%s_global.csv"
infected = read.csv(sprintf(dataset, "confirmed"), check.names=FALSE)
recovered = read.csv(sprintf(dataset,"recovered"), check.names=FALSE)
death = read.csv(sprintf(dataset,"deaths"), check.names=FALSE)
# We analyze data of Afghanistan
infected = unlist(infected[1, -(1:4)])
removed = unlist(recovered[1, -(1:4)] + death[1, -(1:4)])
first.case = which(infected > 0)[1]
# We only analyze first 40 days
infected = infected[first.case:(first.case+40)]
removed = removed[first.case:(first.case+40)]
observed = data.frame(
	date = strptime(names(infected), "%m/%d/%y"),
	S = pop - infected - removed,
	I = infected,
	R = removed,
	row.names = NULL
)

# Population of Afghanistan
pop <- 38928346
nsamples <- 10000
burnin <- 100
samples <- MCMCmetrop(
	# init.state = c(beta=2.491471, gamma=2.392113),
	init.state = c(beta=9.70820698, gamma=9.56283767),
	burnin = burnin,
	nsamples,
	log.target = function(theta) log.prior(theta) + log.likelihood(theta, observed, pop)
)
R0 <- sum(samples[, "beta"] / samples[, "gamma"]) / nsamples
print(R0)
