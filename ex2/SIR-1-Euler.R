Euler_SIR <- function(S0, I0, R0, beta, gamma, last, dt) {
  S <- c(); S[1] <- S0
  I <- c(); I[1] <- I0
  R <- c(); R[1] <- R0
  
  t <- seq(0, last, by = dt)
  week <- seq(1, last / dt + 1, 1 / dt)
  
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
  
  plot(S~t, col="red")
  lines(S~t, col="red")
  lines(I~t, col="green")
  points(I~t, col="green")
  lines(R~t, col="blue")
  points(R~t, col="blue")
  which.max(I) * dt
  
  return(cbind(S[week], I[week], R[week]))
}

Euler_SIR(995, 5, 0, 0.001407, 0.6, 22, 0.5)
