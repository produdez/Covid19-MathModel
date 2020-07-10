S <- c()
I <- c()
R <- c()
S[1] <- 995
I[1] <- 5
R[1] <- 0
dt <- 0.01
t <- seq(0, 22, by = dt)
week <- seq(1, 22 / dt + 1, 1 / dt)
beta <- 0.001407
gamma <- 0.6

dSdt <- function(S, I){
  return (- beta * S * I)
}
dIdt <- function(S, I){
  return (- gamma * I + beta * S * I)
}
dRdt <- function(I){
  return (gamma * I)
}
#Modified Euler's
for (i in 1:(length(t) - 1)) {
  S.k1 <- dSdt(S[i], I[i])
  I.k1 <- dIdt(S[i], I[i])
  R.k1 <- dRdt(I[i])
  
  S.k2 <- dSdt(S[i] + dt * S.k1, I[i] + dt * I.k1)
  I.k2 <- dIdt(S[i] + dt * S.k1, I[i] + dt * I.k1)
  R.k2 <- dRdt(I[i] + dt * I.k1)
  
  S[i + 1] <- S[i] + dt / 2 * (S.k1 + S.k2)
  I[i + 1] <- I[i] + dt / 2 * (I.k1 + I.k2)
  R[i + 1] <- R[i] + dt / 2 * (R.k1 + R.k2)
}

cbind(S[week], I[week], R[week])
plot(S~t, col="red")
lines(S~t, col="red")
lines(I~t, col="green")
points(I~t, col="green")
lines(R~t, col="blue")
points(R~t, col="blue")
which.max(I) * dt