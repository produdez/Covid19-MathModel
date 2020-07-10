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
#Euler's
for (i in 1:(length(t) - 1)) {
  S[i+1] = S[i] + dSdt(S[i], I[i]) * dt
  I[i+1] = I[i] + dIdt(S[i], I[i]) * dt
  R[i+1] = R[i] + dRdt(I[i])       * dt
}

cbind(S[week], I[week], R[week])
plot(S~t, col="red")
lines(S~t, col="red")
lines(I~t, col="green")
points(I~t, col="green")
lines(R~t, col="blue")
points(R~t, col="blue")
which.max(I) * dt