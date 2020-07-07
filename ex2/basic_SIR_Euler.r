S <- c()
I <- c()
R <- c()

S[1] <- 800
I[1] <- 7
R[1] <- 0
delta_t <- 1

t <- seq(0, 8, by = delta_t)
beta <- 0.002
gamma <- 0.5

for (i in 1:(length(t)-1)) {
  S[i+1] = S[i] - (beta*I[i]*S[i]) * delta_t
  
  I[i+1] = I[i] + (beta*I[i]*S[i] - gamma*I[i]) * delta_t
  
  R[i+1] = R[i] + (gamma*I[i]) * delta_t
}


plot(S~t, col="red")
lines(S~t, col="red")
lines(I~t, col="green")
points(I~t, col="green")
lines(R~t, col="blue")
points(R~t, col="blue")

text(7, 800, paste0("S(t) = ", round(S[length(t)], 4)), col = "red")
text(7, 750, paste0("I(t) = ", round(I[length(t)], 4)), col = "green")
text(7, 700, paste0("R(t) = ", round(R[length(t)], 4)), col = "blue")

