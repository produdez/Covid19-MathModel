Euler_SIR <- function(S0, I0, R0, beta, gamma, t, delta_t) {
  S <- c(); S[1] <- S0
  I <- c(); I[1] <- I0
  R <- c(); R[1] <- R0

  times <- seq(0, t, by = delta_t)
  
  for (i in 1:(length(times)-1)) {
    S[i+1] = S[i] - (beta * I[i] * S[i]) * delta_t
    
    I[i+1] = I[i] + (beta * I[i] * S[i] - gamma * I[i]) * delta_t
    
    R[i+1] = R[i] + (gamma * I[i]) * delta_t
  }
  
  result <- c(round(S[length(times)], 4),
              round(I[length(times)], 4),
              round(R[length(times)], 4))
  
  plot(S~times, col="red")
  lines(S~times, col="red")
  lines(I~times, col="green")
  points(I~times, col="green")
  lines(R~times, col="blue")
  points(R~times, col="blue")
  return(result)
}



Euler_SIR(800, 7, 0, 0.002, 0.5, 8, 1)
