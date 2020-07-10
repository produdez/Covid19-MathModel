# Arguments:
#   y0: the initial values for the ODE system
#   times: times at which explicit estimates for y are desired, the first value must be the initial time
#   func: a function that computes the values of the derivatives in the ODE system at time t
#   parms: vector or list of parameters used in `func`.
euler <- function(y0, times, func, parms) {
	h <- times[2] - times[1]
	tbl <- matrix(nrow=length(times), ncol=1+length(y0))
	colnames(tbl) <- c("time", names(y0))
	tbl[, 1] <- times
	tbl[1, -1] <- y0
	for (i in 1:(nrow(tbl)-1)) {
		ti <- tbl[i, 1]
		yi <- tbl[i, -1]
		tbl[i+1, -1] <- yi + func(ti, yi, parms)
	}
	return(tbl)
}

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

result <- rk4(
	y0 = c(S=99990, I=10, R=0),
	times = seq(0, 200, by=1),
	func = SIR.model,
	parms = list(beta=.5, gamma=1/3, N=100000)
)

with(as.list(data.frame(result)), {
	plot(S ~ time, col="red", type="l", ylim=c(0, 100000), xlab="Day", ylab="People")
	lines(I ~ time, col="green")
	lines(R ~ time, col="blue")
	legend(150, 100000, legend=c("Susceptible", "Infected", "Removed"),
	    col=c("red", "green", "blue"), lty=1)
})
