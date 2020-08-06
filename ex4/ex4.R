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

CR.model <- function(t, y, parms) with(as.list(c(y, parms)), {
	dC.dt <- beta / N * (N - C) * (C - R)
	dR.dt <- gamma * (C - R)
	c(dC.dt, dR.dt)
})

rw.metro <- function(start, func, sample.size, pcov.scale=1, ...) {
	# Dimension of state space
	dim <- length(start)

	opt <- optim(start, func, control=list(fnscale=-1), hessian=TRUE, ...)
	if (opt$convergence != 0) stop(opt$message)
	prop.cov <- -2.4^2 / dim * solve(opt$hessian)
	# Current state of the Markov chain
	theta <- mvrnorm(n=1, mu=opt$par, Sigma=prop.cov)
	ptheta <- func(theta, ...)
	# Generate matrix containing the sample. Initialize first sample with the starting value
	sample <- matrix(nrow=sample.size, ncol=dim)
	colnames(sample) <- names(start)
	sample[1, ] <- theta
	# Generate uniform random numbers in advance, to save computation.
	log.u <- log(runif(sample.size - 1))
	# Proposal is a multivariate standard normal distribution. Generate sample and
	# later on use linearity property of Gaussian distribution
	normal.shift <- mvrnorm(n=sample.size-1, mu=rep(0, dim), Sigma=prop.cov)
	for (i in seq_len(sample.size - 1)) {
		# Sample a candidate
		candidate <- theta + normal.shift[i, ]
		# Calculate func of candidate and store it in case it gets accepted
		log.r <- func(candidate, ...) - ptheta
		if (log.u[i] < log.r) {
			theta <- candidate
			ptheta <- func(theta, ...)
		}
		sample[i+1, ] <- theta
	}
	return(sample)
}

# Prior distribution for our parameter
log.prior <- function(theta) sum(dnorm(theta, sd=100, log=TRUE))

log.likelihood <- function(theta, data, pop) {
	n <- 20
	ninterval <- nrow(data) - 1
	change <- data[1:ninterval, c("S.I", "I.R")]
	estimated <- t(apply(
		data[1:ninterval, c("C", "R")],
		MARGIN = 1,
		FUN = function(x) {
			sol <- rk4(
				y0 = x,
				times = seq(0, 1, by=1/n),
				func = CR.model,
				parms=c(exp(theta), N=pop)
			)
			sol[n+1, -1] - sol[1, -1]
		}
	))
	sum(dpois(change, estimated, log=TRUE))
}

log.posterior <- function(theta, data, pop) {
	log.prior(theta) + log.likelihood(theta, data, pop)
}

dataset <- paste(
	"https://github.com/CSSEGISandData/COVID-19",
	"raw/master/csse_covid_19_data/csse_covid_19_time_series",
	"time_series_covid19_%s_global.csv",
	sep = "/"
)
confirmed <- read.csv(sprintf(dataset, "confirmed"), check.names=FALSE)
recovered <- read.csv(sprintf(dataset, "recovered"), check.names=FALSE)
deaths <- read.csv(sprintf(dataset, "deaths"), check.names=FALSE)
report.date <- as.Date(names(confirmed[, -(1:4)]), "%m/%d/%y")

confirmed.VN <- unlist(subset(confirmed, `Country/Region` == "Vietnam", select=-(1:4)))
recovered.VN <- unlist(subset(recovered, `Country/Region` == "Vietnam", select=-(1:4)))
deaths.VN <- unlist(subset(deaths, `Country/Region` == "Vietnam", select=-(1:4)))
observed <- data.frame(
	date = report.date,
	C = confirmed.VN,
	R = recovered.VN + deaths.VN,
	row.names = NULL
)
observed$I.R <- c(diff(observed$R), NA)
observed$S.I <- c(diff(observed$C), NA)
# Population of Vietnam
pop <- 97338579

interval1 <- seq(as.Date("2020-03-25"), by = "day", length.out = 7)
obs1 <- subset(observed, date %in% interval1)

init.beta <- (obs1$C[2] - obs1$C[1]) * pop / (pop - obs1$C[1]) / (obs1$C[1] - obs1$R[1])
init.gamma <- (obs1$R[2] - obs1$R[1]) / (obs1$C[1] - obs1$R[1])
sample.size <- 50000
sample1 <- rw.metro(
	start = log(c(beta=init.beta, gamma=init.gamma)),
	func = log.posterior,
	sample.size = sample.size,
	data = data.matrix(obs1[, -1]),
	pop = pop
)

sample1.R0 <- exp(sample1[, "beta"] - sample1[, "gamma"])
R0 <- mean(sample1.R0)
cat("Bayes estimate of R0 is", R0, "\n")

blen <- floor(sqrt(sample.size))
nbatch <- sample.size %/% blen
batch.means <- sapply(
    1:nbatch,
    FUN = function(k) mean(sample1.R0[((k-1)*blen + 1):(k * blen)])
)
var.hat <- blen / nbatch * sum((batch.means - R0)^2)
err1 <- sqrt(var.hat / sample.size)
cat("Monte Carlo standard error of above approximation is", err1, "\n")

interval2 <- seq(as.Date("2020-04-08"), by = "day", length.out = 7)
obs2 <- subset(observed, date %in% interval2)

init.beta <- (obs2$C[2] - obs2$C[1]) * pop / (pop - obs2$C[1]) / (obs2$C[1] - obs2$R[1])
init.gamma <- (obs2$R[2] - obs2$R[1]) / (obs2$C[1] - obs2$R[1])
sample2 <- rw.metro(
	start = log(c(beta=init.beta, gamma=init.gamma)),
	func = log.posterior,
	sample.size = sample.size,
	data = data.matrix(obs2[, -1]),
	pop = pop
)

sample2.R0 <- exp(sample2[, "beta"] - sample2[, "gamma"])
R0 <- mean(sample2.R0)
cat("Bayes estimate of R0 is", R0, "\n")

batch.means <- sapply(
    1:nbatch,
    FUN = function(k) mean(sample2.R0[((k-1)*blen + 1):(k * blen)])
)
var.hat <- blen / nbatch * sum((batch.means - R0)^2)
err2 <- sqrt(var.hat / sample.size)
cat("Monte Carlo standard error of above approximation is", err2, "\n")
