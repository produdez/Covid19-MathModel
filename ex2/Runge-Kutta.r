library(reshape2)
library(ggplot2)

dSdt <- function(t, S, I) {
  return(-beta * S * I / N)
}

dIdt <- function(t, S, I) {
  return(beta * S * I / N - gamma * I)
}

dRdt <- function(t, I) {
  return(gamma * I)
}

RK4_SIR <- function(t, beta, gamma, S0, I0, R0 = 0, dt = 1) {
  N <<- S0 + I0 + R0  # fixed population
  
  S <- c(S0, rep(0, t))
  I <- c(I0, rep(0, t))
  R <- c(R0, rep(0, t))
  for (i in 1:t) {
    
    S.k1 <- dSdt(i, S[i], I[i])
    I.k1 <- dIdt(i, S[i], I[i])
    R.k1 <- dRdt(i, I[i])
    
    S.k2 <- dSdt(i + dt / 2, S[i] + dt / 2 * S.k1, I[i] + dt / 2 * I.k1)
    I.k2 <- dIdt(i + dt / 2, S[i] + dt / 2 * S.k1, I[i] + dt / 2 * I.k1)
    R.k2 <- dRdt(i + dt / 2, I[i] + dt / 2 * I.k1)
    
    S.k3 <- dSdt(i + dt / 2, S[i] + dt / 2 * S.k2, I[i] + dt / 2 * I.k2)
    I.k3 <- dIdt(i + dt / 2, S[i] + dt / 2 * S.k2, I[i] + dt / 2 * I.k2)
    R.k3 <- dRdt(i + dt / 2, I[i] + dt / 2 * I.k2)
    
    S.k4 <- dSdt(i + dt, S[i] + dt * S.k3, I[i] + dt * I.k3)
    I.k4 <- dIdt(i + dt, S[i] + dt * S.k3, I[i] + dt * I.k3)
    R.k4 <- dRdt(i + dt, I[i] + dt * I.k3)
    
    S[i + 1] <- S[i] + dt / 6 * (S.k1 + 2 * S.k2 + 2 * S.k3 + S.k4)
    I[i + 1] <- I[i] + dt / 6 * (I.k1 + 2 * I.k2 + 2 * I.k3 + I.k4)
    R[i + 1] <- R[i] + dt / 6 * (R.k1 + 2 * R.k2 + 2 * R.k3 + R.k4)
  }
  R <- N - S - I
  return(data.frame(n = 0:t, S = S, I = I, R = R))
}


r <- RK4_SIR(200, 0.5, 1/3, 99990, 10)


r.plot <- melt(r, id = "n", measure = c("S", "I", "R"))

p <- ggplot(r.plot, aes(x = n, y = value, group = variable, color = variable))
p + geom_line()