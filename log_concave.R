rm(list = ls())


##### Preliminary functions #####
library(nloptr)
library(logcondens)

Phi <- function(eta) {
  #eta <- abs(eta)
  
  # vectorized
  out <- rep(1, length(eta))
  idx <- (eta != 0)
  out[idx] <- (1/eta[idx])*log(1 + eta[idx])
  return(out)
}

integral1 <- function(delta) {
  # vectorized
  out <- rep(1/2, length(delta))
  idx <- (delta != 0)
  out[idx] <- ((1 + delta[idx])*log(1 + delta[idx]) - delta[idx])/delta[idx]^2
  return(out)
}

integral2 <- function(delta) {
  # vectorized
  out <- rep(1/3, length(delta))
  idx <- (delta != 0)
  out[idx] <- (2*delta[idx] + (1 + delta[idx])*(log(1 + delta[idx]))^2 - 2*(1 + delta[idx])*log(1 + delta[idx]))/delta[idx]^3
  return(out)  
}

Qch <- function(C, h, u_seq) {
  # computes q_0, q_1, ...., q_K
  # h = (h_0, h_1, ..., h_K), length = 1 + K
  
  delta <- h[-1]/h[-length(h)] - 1  # delta_1, ..., delta_K
  Du <- diff(u_seq)  # Du_1, ..., Du_K
  q <- C + c(0, cumsum( (Du/h[-length(h)])*Phi(delta)))
  return(q)  # q_0, q_1, ..., q_K
}

W2_loss <- function(C, h, y, u_seq) {
  # W2 loss function
  
  Du <- diff(u_seq)  # Du_1, ..., Du_K
  delta <- h[-1]/h[-length(h)] - 1  # delta_1, ..., delta_K
  q <- Qch(C, h, u_seq)  # q_0, ..., q_K
  q <- q[-length(q)]  # q_0, ..., q_{K-1}
  
  I1 <- (q - y)^2  
  I2 <- ((2*(q - y)*Du)/h[-length(h)])*integral1(delta)
  I3 <- ((Du)^2/h[-length(h)]^2)*integral2(delta)
  W2 <- sum((I1 + I2 + I3)*Du)
  return(W2)
}

g <- function(x) {
  # Constraint (g(x) <= 0)
  # Here we assume u_seq is uniform
  C <- x[1]  # constant
  h <- x[2:length(x)]  # derivative of T
  return(c(diff(h, differences = 2)))
}


logcondens_quantile <- function(fit) {
  # Extract knots and log-density values
  x   <- fit$x
  phi <- fit$phi
  
  n <- length(x)
  if (n < 2)
    stop("Need at least two support points.")
  
  # Slopes and intercepts for log-density on each interval
  a <- diff(phi) / diff(x)
  b <- phi[-n] - a * x[-n]
  
  # Segment masses (unnormalized)
  seg_mass <- numeric(n - 1)
  
  for (i in seq_len(n - 1)) {
    if (abs(a[i]) > 1e-12) {
      seg_mass[i] <- (exp(a[i] * x[i + 1] + b[i]) -
                        exp(a[i] * x[i] + b[i])) / a[i]
    } else {
      # nearly constant density
      seg_mass[i] <- exp(b[i]) * (x[i + 1] - x[i])
    }
  }
  
  # Cumulative unnormalized mass
  cum_mass <- c(0, cumsum(seg_mass))
  total_mass <- cum_mass[length(cum_mass)]
  
  # Return quantile function
  function(p) {
    if (any(p < 0 | p > 1))
      stop("p must be in [0, 1]")
    
    sapply(p, function(pp) {
      if (pp == 0) return(x[1])
      if (pp == 1) return(x[n])
      
      target <- pp * total_mass
      
      # Find segment j such that
      # cum_mass[j] <= target < cum_mass[j+1]
      j <- max(which(cum_mass <= target))
      if (j == length(cum_mass))
        j <- j - 1
      
      base_mass <- cum_mass[j]
      delta <- target - base_mass
      
      # Invert within segment j
      if (abs(a[j]) > 1e-12) {
        fx0 <- exp(a[j] * x[j] + b[j])
        x[j] + (1 / a[j]) *
          log(1 + a[j] * delta / fx0)
      } else {
        x[j] + delta / exp(b[j])
      }
    })
  }
}




##### Mixture of two points #####
## Case 1
n <- 50  # K = n  # can be larger if we also implement the gradient
u_seq <- seq(0, 1, by = 1/n)  # uniform partition based on data
m <- n*0.4
y <- c(rep(-1, n - m), rep(1, m))

f <- function(x) {
  # wrapper of the objective
  C <- x[1]
  h <- x[2:length(x)]
  return(W2_loss(C, h, y, u_seq))
}  

# Let initial guesss = uniform distribution over [min x_i, max]
C0 <- mean(y)
h0 <- (1/(0.5*(max(y) - min(y))))*rep(1, 1 + n)
x0 <- c(C0, h0)

#pdf(file = "logconcave1.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))
# compute Wasserstein projection
opts <- list("algorithm" = "NLOPT_LN_COBYLA",
             "maxeval" = 50000,
             "xtol_rel" = 1e-6,
             "eval_g_ineq" = 1e-10)
result <- nloptr(x0 = x0, eval_f = f, eval_g_ineq = g,
                 opt = opts, lb = c(-Inf, rep(.Machine$double.eps, length(u_seq))))
C_opt <- result$solution[1]
h_opt <- result$solution[-1]
q_seq <- Qch(C_opt, h_opt, u_seq)
n_points <- 10  # ad hoc parameter used for plotting the density
Du <- diff(u_seq)
K <- length(u_seq) - 1
my_col <- grey.colors(1, start = 0.4, end = 0.4)
my_col <- paste(my_col, "70", sep = "")

plot(NULL, NULL, xlim = c(-1.7, 1.7), ylim = c(0, 2.5),
     xlab = "x", ylab = "", main = "Density")
legend(x = "topright",
       legend = c("Wasserstein", "MLE"), lty = c(1, 2), 
       col = c("grey", "blue"), lwd = c(3, 2))
abline(h = 0)
# plot Wasserstein
for (i in 1:K) {
  q_im <- q_seq[i]
  h_im <- h_opt[i]
  delta_i <- h_opt[i+1]/h_opt[i] - 1
  Du_i <- Du[i]
  x_seq <- seq(q_seq[i], q_seq[i + 1], length.out = n_points)
  f_seq <- h_im*exp(((h_im*delta_i)/Du_i)*(x_seq - q_im))
  polygon(c(q_seq[i], x_seq, q_seq[i + 1]), c(0, f_seq, 0), col = my_col, border = NA)
}

# plot MLE
MLE <- logConDens(y, smoothed = FALSE)
x_seq <- seq(-1, 1, by = 0.01)
points(x_seq, evaluateLogConDens(x_seq, MLE)[, 3], type = "l",
       lty = 2, lwd = 2, col = "blue")
segments(-1, 0, -1, evaluateLogConDens(-1, MLE)[, 3], lty = 2, lwd = 2, col = "blue")
segments(1, 0, 1, evaluateLogConDens(1, MLE)[, 3], lty = 2, lwd = 2, col = "blue")
points(-1, 0, pch = 4)
points(1, 0, pch = 4)

## Case 2
m <- n*0.2
y <- c(rep(-1, n - m), rep(1, m))

f <- function(x) {
  # wrapper of the objective
  C <- x[1]
  h <- x[2:length(x)]
  return(W2_loss(C, h, y, u_seq))
}  

C0 <- mean(y)
h0 <- (1/(0.5*(max(y) - min(y))))*rep(1, 1 + n)
x0 <- c(C0, h0)

result <- nloptr(x0 = x0, eval_f = f, eval_g_ineq = g,
                 opt = opts, lb = c(-Inf, rep(.Machine$double.eps, length(u_seq))))
C_opt <- result$solution[1]
h_opt <- result$solution[-1]
q_seq <- Qch(C_opt, h_opt, u_seq)
n_points <- 10
Du <- diff(u_seq)
K <- length(u_seq) - 1

my_col <- grey.colors(1, start = 0.4, end = 0.4)
my_col <- paste(my_col, "70", sep = "")

plot(NULL, NULL, xlim = c(-1.7, 1.7), ylim = c(0, 2.5),
     xlab = "x", ylab = "", main = "Density")
legend(x = "topright",
       legend = c("Wasserstein", "MLE"), lty = c(1, 2), 
       col = c("grey", "blue"), lwd = c(3, 2))
abline(h = 0)
for (i in 1:K) {
  q_im <- q_seq[i]
  h_im <- h_opt[i]
  delta_i <- h_opt[i+1]/h_opt[i] - 1
  Du_i <- Du[i]
  x_seq <- seq(q_seq[i], q_seq[i + 1], length.out = n_points)
  f_seq <- h_im*exp(((h_im*delta_i)/Du_i)*(x_seq - q_im))
  polygon(c(q_seq[i], x_seq, q_seq[i + 1]), c(0, f_seq, 0), col = my_col, border = NA)
}

MLE <- logConDens(y, smoothed = FALSE)
x_seq <- seq(-1, 1, by = 0.01)
points(x_seq, evaluateLogConDens(x_seq, MLE)[, 3], type = "l",
       lty = 2, lwd = 2, col = "blue")
segments(-1, 0, -1, evaluateLogConDens(-1, MLE)[, 3], lty = 2, lwd = 2, col = "blue")
segments(1, 0, 1, evaluateLogConDens(1, MLE)[, 3], lty = 2, lwd = 2, col = "blue")

points(-1, 0, pch = 4)
points(1, 0, pch = 4)
#dev.off()





##### Misspecified case #####

library(mistr)  # compute quantile of mixture distribution
my_mixture <- mixdist(gammadist(shape = 5, rate = 1),
                      gammadist(shape = 20, rate = 1),
                      weights = c(1/4, 3/4))
n <- 50
u_seq <- seq(0, 1, by = 1/n)
set.seed(2028)
y <- r(my_mixture, n)  # samples from mixture distribution
y <- sort(y, decreasing = FALSE)

f <- function(x) {
  # wrapper of the objective
  C <- x[1]
  h <- x[2:length(x)]
  return(W2_loss(C, h, y, u_seq))
}  

# Let initial guesss = uniform distribution over [min x_i, max]
C0 <- mean(y)
h0 <- (1/(max(y) - min(y)))*rep(1, 1 + n)
x0 <- c(C0, h0)

opts <- list("algorithm" = "NLOPT_LN_COBYLA", "maxeval" = 50000,
             "xtol_rel" = 1e-10, "xtol_abs" = 1e-10)
result <- nloptr(x0 = x0, eval_f = f, eval_g_ineq = g,
                 opt = opts, lb = c(-Inf, rep(.Machine$double.eps, length(u_seq))))
C_opt <- result$solution[1]
h_opt <- result$solution[-1]
MLE <- logConDens(y, smoothed = FALSE)

#pdf(file = "logconcave2.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))
Q_opt <- Qch(C_opt, h_opt, u_seq)
plot(u_seq, Q_opt, xlab = "u", ylab = "", type = "l", ylim = c(0, 32),
     col = "grey", lwd = 3, main = "Quantile function")
Q0 <- quantile(y, u_seq, type = 1)
points(u_seq, Q0, type = "S")
points(c(u_seq[-length(u_seq)], 0.999, 0.9999), mistr::q(my_mixture, c(u_seq[-length(u_seq)], 0.999, 0.9999)), type = "l")
q_MLE <- logcondens_quantile(MLE)
points(u_seq, q_MLE(u_seq), type = "l", lty = 2, lwd = 2, col = "blue")

legend(x = "topleft",
       legend = c("data/true", "Wasserstein", "MLE"),
       lwd = c(1, 3, 2),
       lty = c(1, 1, 2), col = c("black", "grey", "blue"))


plot(NULL, NULL, xlim = c(-2, 34), ylim = c(0, 0.07),
     main = "Density", xlab = "x", ylab = "")
q_seq <- Qch(C_opt, h_opt, u_seq)
n_points <- 5
Du <- diff(u_seq)
K <- length(u_seq) - 1
for (i in 1:K) {
  q_im <- q_seq[i]
  h_im <- h_opt[i]
  delta_i <- h_opt[i+1]/h_opt[i] - 1
  Du_i <- Du[i]
  x_seq <- seq(q_seq[i], q_seq[i + 1], length.out = n_points)
  f_seq <- h_im*exp(((h_im*delta_i)/Du_i)*(x_seq - q_im))
  points(x_seq, f_seq, type = "l", col = "darkgrey", lwd = 3)
  if (i == 1) {
    segments(q_seq[i], 0, q_seq[i], f_seq[1], col = "darkgrey", lwd = 3)
  }
  if (i == K) {
    segments(q_seq[i + 1], 0, q_seq[i + 1], f_seq[n_points], col = "darkgrey", lwd = 3)
  }
}

s_seq <- seq(0, 37, by = 0.01)
points(s_seq, d(my_mixture, s_seq), type = "l")
x_seq <- seq(0, 37, by = 0.01)
points(x_seq, evaluateLogConDens(x_seq, MLE)[, 3], type = "l",
       lty = 2, lwd = 2, col = "blue")
segments(-1, 0, -1, evaluateLogConDens(-1, MLE)[, 3], lty = 2, lwd = 2, col = "blue")
segments(1, 0, 1, evaluateLogConDens(1, MLE)[, 3], lty = 2, lwd = 2, col = "blue")
abline(h = 0)
y_matrix <- cbind(y, rep(0, length(y)))
points(y_matrix[, 1], y_matrix[, 2], pch = 4)
legend(x = "topleft",
       legend = c("true", "Wasserstein", "MLE"),
       lwd = c(1, 3, 2),
       lty = c(1, 1, 2), col = c("black", "grey", "blue"))
#dev.off()
