library(MASS)

rk4 <- function(y0, times, func, parms=NULL) {
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

rw.metro <- function(start, func, nsamples, pcov.scale=1, ...) {
	opt <- optim(start, func, control=list(fnscale=-1), hessian=TRUE, ...)
	if (opt$convergence != 0) stop(opt$message)
	V <- -pcov.scale^2 * solve(opt$hessian)
	print(opt)
	# Dimension of state space
	dim <- length(start)
	# Current state of the Markov chain
	theta <- opt$par
	ptheta <- func(theta, ...)
	# Generate matrix containing the samples. Initialize first sample with the starting value
	samples <- matrix(nrow=nsamples, ncol=dim)
	colnames(samples) <- names(start)
	# Generate uniform random numbers in advance, to save computation.
	log.u <- log(runif(nsamples))
	# Proposal is a multivariate standard normal distribution. Generate samples and
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

# Prior distribution for our parameter
log.prior <- function(theta) sum(dnorm(theta, sd=100, log=TRUE))

log.likelihood <- function(theta, data, pop) {
	n <- 10
	ninterval <- nrow(data) - 1
	change <- data[1:ninterval, c("S.I", "I.R")]
	estimated <- t(apply(
		data[1:ninterval, c("I", "R")],
		MARGIN = 1,
		FUN = function(x) {
			sol <- rk4(
				y = c(S = pop - sum(x), x),
				times = seq(0, 1, by=1/n),
				func = SIR.model,
				parms=c(exp(theta), N=pop)
			)
			flow <- sol[n+1, ] - sol[1, ]
			c(S.I = -flow[["S"]], I.R = flow[["R"]])
		}
	))
	sum(dpois(change, estimated, log=TRUE))
}

log.posterior <- function(theta, obs, population) {
	log.prior(theta) + log.likelihood(theta, obs, population)
}

dataset <- "https://github.com/CSSEGISandData/COVID-19/raw/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_%s_global.csv"
confirmed <- read.csv(sprintf(dataset, "confirmed"), check.names=FALSE)
recovered <- read.csv(sprintf(dataset, "recovered"), check.names=FALSE)
death <- read.csv(sprintf(dataset, "deaths"), check.names=FALSE)
# We analyze data of Afghanistan
removed <- unlist(recovered[1, -(1:4)] + death[1, -(1:4)])
infected <- unlist(confirmed[1, -(1:4)]) - removed
observed <- data.frame(
	date = as.Date(names(infected), "%m/%d/%y"),
	I = infected,
	R = removed,
	row.names = NULL
)
observed$I.R <- c(diff(observed$R), NA)
observed$S.I <- c(diff(observed$I), NA) + observed$I.R
interval1 <- seq(as.Date("2020-05-01"), by = "day", length.out = 7)
obs1 <- subset(observed, date %in% interval1)

pop <- 38928346
nsamples <- 100
samples <- rw.metro(
	start = log(c(beta=.1, gamma=.01)),
	func = log.posterior,
	nsamples,
	pcov.scale = 2,
	obs = data.matrix(obs1[, -1]),
	population = pop
)
plot(samples, type="l")
R0 <- sum(exp(samples[, "beta"] - samples[, "gamma"])) / nrow(samples)
cat("R0 is", R0, "\n")
rej.rate <- sum(apply(diff(samples) == 0, 1, FUN=all)) / (nrow(samples) - 1)
cat("Rejection rate is", rej.rate, "\n")
