rm(list = ls())



##### Preliminary functions #####
library(quadprog) # quadratic programming
library(fdrtool) # for Grenander


construct_D_matrix <- function(u_seq, y) {
  # D matrix in solve.QP()
  # u_seq = (u0 = 0, u1, ..., uK = 1)  # given partition
  # y = (y1, ..., yK)  # data points (ordered)
  
  K <- length(y)
  Du <- c(diff(u_seq), 0)  # length K + 1
  
  # construct D matrix
  D <- matrix(0, nrow = K, ncol = K) 
  for (i in 1:K) {
    D[i, i] <- 2*(Du[i] + Du[i + 1])/3
  }
  for (i in 1:(K-1)) {
    D[i, i + 1] <- Du[i + 1]/3
    D[i + 1, i] <- Du[i + 1]/3
  }
  return(D)
}

construct_A_matrix <- function(u_seq) {
  # A matrix (for constraint (A^T)z >= 0) in solve.QP()
  
  K <- length(u_seq) - 1

  # monotonicity constraint
  B1 <- matrix(0, nrow = K, ncol = K)
  B1[1, 1] <- 1   # z_1 >= 0 = z_0
  for (i in 2:K) {
    # z_i - z_{i-1} >= 0
    B1[i, i - 1] <- -1
    B1[i, i] <- 1  
  }
  
  # convexity constraint
  Du <- diff(u_seq)  # (Du)_i = u_i - u_{i-1}, length = K
  B2 <- matrix(0, nrow = K - 1, ncol = K)
  for (i in 1:(K-1)) {
    B2[i, i] <- -1 
    B2[i, i + 1] <- Du[i]/(Du[i] + Du[i + 1])
  }
  for (i in 2:(K-1)) {
    B2[i, i-1] <- Du[i + 1]/(Du[i] + Du[i + 1])
  }
  
  AT <- rbind(B1, B2)
  A  <- t(AT)
  return(A)
}

construct_d_vector <- function(u_seq, y) {
  # d vector in solve.QP()
  
  K <- length(y)
  Du <- c(diff(u_seq), 0)  # length K + 1
  y <- c(y, y[length(y)])
  d <- numeric(K)
  for (i in 1:K) {
    d[i] <- y[i]*Du[i] + y[i+1]*Du[i + 1]
  }
  return(d)
}

construct_c_const <- function(u_seq, y) {
  # constant term in W2
  Du <- diff(u_seq)
  c <- sum((y^2)*Du)
  return(c)
}

W2_monotone <- function(x, u_seq = NULL) {
  # Wasserstein projection (p = 2) in the monotone case

  x <- sort(x, decreasing = FALSE)  # x in ascending order
  
  # set up data for quadratic programming problem
  if (is.null(u_seq)) {
    # use partition with mesh size 1/n if u_seq is not given
    n <- length(x) 
    u_seq <- seq(0, 1, length.out = 1 + n)
  } 
  y <- as.numeric(quantile(x, u_seq[-1], type = 1))  # inverse of empirical cdf
  
  D <- construct_D_matrix(u_seq, y)
  A <- construct_A_matrix(u_seq)
  d <- construct_d_vector(u_seq, y)
  b <- rep(0, dim(A)[2])
  
  # QP and return results
  QP_result <- solve.QP(D, d, A, b)
  Q_opt <- c(0, QP_result$solution)  # quantile w.r.t. u_seq
  
  c <- construct_c_const(u_seq, y)
  W2 <- QP_result$value + c  # sqaured Wasserstein distance
  
  return(list(Q = Q_opt, u = u_seq, W2 = W2))
}

monotone_density <- function(estimate, draw = TRUE) {
  # Compute monotone (piecewise constant) density
  # given output
  
  Q <- estimate$Q
  u <- estimate$u
  
  DQ <- diff(Q)
  Du <- diff(u)
  f <- Du/DQ
  x <- Q[-1]

  # plot density
  if (draw) {
    K <- length(u) - 1
    plot(NULL, NULL, xlim = c(0, max(x)), ylim = c(0, max(f)),
         xlab = "x", ylab = "", main = "Density")
    
    rect(0, 0, x[1], f[1], col = rgb(0, 0, 0, 0.3), border = NA)
    for (i in 2:K) {
      rect(x[i - 1], 0, x[i], f[i], col = rgb(0, 0, 0, 0.3), border = NA)
    }
  }
    
  return(list(x = x, f = f))
}


Grenander <- function(x) {
  # density of Grenander estimator
  
  x <- sort(as.numeric(x))
  n <- length(x)
  
  # empirical CDF including origin anchor
  ux <- unique(x)
  counts <- tabulate(match(x, ux))
  
  # force 0 at origin
  if (min(x) > 0) {
    X <- c(0, ux)  
    Y <- c(0, cumsum(counts)/n)
  } else {
    X <- ux
    Y <- cumsum(counts)/n
    Y[1] <- 0
  }
    
  # compute least concave majorant
  gcm <- gcmlcm(X, Y, type = "lcm")

  left <- head(gcm$x, -1)
  right <- tail(gcm$x, -1)
  height <- diff(gcm$y) / diff(gcm$x)
  
  # density function
  d <- function(s) {
    s <- as.numeric(s)
    out <- numeric(length(s))
    
    for (i in seq_along(height)) {
      out[s > left[i] & s <= right[i]] <- height[i]
    }
    
    out[s <= 0] <- height[1]
    out[s > max(right)] <- 0
    
    out
  }
  
  cdf_left <- c(0, cumsum(height * (right - left)))
  
  # quantile function
  q <- function(p) {
    p <- as.numeric(p)
    out <- numeric(length(p))
    
    for (i in seq_along(p)) {
      if (p[i] <= 0) {
        out[i] <- left[1]
      } else if (p[i] >= 1) {
        out[i] <- right[length(right)]
      } else {
        # find which interval p falls into
        j <- max(which(cdf_left <= p[i]))
        # linear interpolation
        out[i] <- left[j] + (p[i] - cdf_left[j]) / height[j]
      }
    }
    out
  }
  
  return(list(d = d, q = q))
}
#####







##### Comparision between Wasserstein projection and Grenander #####
##### using mixture distributions between two points           #####

u_seq <- seq(0, 1, by = 0.001)
mix_seq <- 0:10  # cases for plotting quantile functions
mix_seq2 <- c(2, 10)  # ad hoc (cases for plotting densities)

col_seq0 <- rev(grey.colors(length(mix_seq), start = 0.05, end = 0.6))
col_seq <- paste(col_seq0, "30", sep = "")
col_seq3 <- paste(col_seq0, "50", sep = "")
col_seq2 <- hcl.colors(2*length(mix_seq), palette = "Blues 3")[1:length(mix_seq)]

#pdf(file = "monotone1.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))

plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1.5),
     xlab = "u", ylab = "", main = "Quantile function (Wasserstein)")
for (i in 1:length(mix_seq)) {
  k <- mix_seq[i]
  x <- c(rep(0.2, 10 - k), rep(1, k))
  my_result <- W2_monotone(x, u_seq)
  n <- length(x)
  #u_seq <- seq(0, 1, length.out = 1 + n)
  Q_opt <- my_result$Q
  
  lwd <- 1
  if (k %in% mix_seq2) {
    lwd <- 3
  }
  points(u_seq, Q_opt, type = "l", lwd = lwd, col = col_seq[i])
}


plot(NULL, NULL, xlim = c(0, 1.5), ylim = c(0, 5),
     xlab = "x", ylab = "", main = "Density")
for (i in 1:length(mix_seq2)) {
  k <- mix_seq2[i]  
  x <- c(rep(0.2, 10 - k), rep(1, k))  # data
  
  # compute and plot Wasserstein projection estimator
  my_result <- W2_monotone(x, u_seq)
  my_density <- monotone_density(my_result, draw = FALSE)
  f <- my_density$f
  x_seq <- my_density$x
  K <- length(my_result$u) - 1

  for (j in 2:K) {
    rect(x_seq[j - 1], 0, x_seq[j], f[j], col = col_seq[k], border = NA)
  }
  segments(0, 0, 0, f[1], col = col_seq3[k], lwd = 2)
  points(x_seq, f, type = "S", col = col_seq[k], lwd = 2)
  segments(x_seq[K], 0, x_seq[K], f[K], col = col_seq3[k], lwd = 2)
}
for (i in 1:length(mix_seq2)) {
  k <- mix_seq2[i]  
  x <- c(rep(0.2, 10 - k), rep(1, k))  # data
  # compute and plot Grenander estimator
  MLE <- Grenander(x)$d
  points(c(0, x), MLE(c(0, x)), type = "S",
         col = col_seq2[k + 1], lty = 2, lwd = 2)
  segments(max(x), MLE(max(x)), max(x), 0,
           col = col_seq2[k + 1], lty = 2, lwd = 2)
  segments(0, 0, 0, MLE(0), lty = 2, lwd = 2, col = col_seq2[k + 1])
}
abline(h = 0) 
points(c(0.2, 1), c(0, 0), pch = 4)
legend(x = "topright",
       legend = c("Wasserstein (0)", "Wasserstein (0.8)", "MLE (0)", "MLE (0.8)"),
       lty = c(1, 1, 2, 2), col = c(col_seq3[10], col_seq3[2], col_seq2[11], col_seq2[3]),
       lwd = c(2, 2, 2, 2))
#dev.off()



##### Misspecified case 2 #####
library(mistr)  # compute quantile of mixture distribution
u_seq <- seq(0, 1, by = 0.001)
my_mixture <- mixdist(gammadist(shape = 5, rate = 1.5),
                      gammadist(shape = 20, rate = 1.5),
                      weights = c(3/5, 2/5))

x <- q(my_mixture, p = seq(0.001, 0.999, length.out = 1000))

#pdf(file = "monotone2.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))
plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 23),
     xlab = "u", ylab = "", main = "Quantile function")

Q0 <- quantile(x, u_seq, type = 1)
points(u_seq, Q0, type = "l")
my_result <- W2_monotone(x, u_seq)
points(u_seq, my_result$Q, type = "l", col = "darkgrey", lwd = 3)

Q_G <- Grenander(x)$q
points(u_seq, Q_G(u_seq), type = "l", lty = 2, lwd = 2, col = "blue")
legend(x = "topleft",
       legend = c("data", "Wasserstein", "Grenander"),
       lwd = c(1, 3, 2),
       lty = c(1, 1, 2), col = c("black", "darkgrey", "blue"))
plot(NULL, NULL, xlim = c(0, 26), ylim = c(0, 0.175),
     xlab = "x", ylab = "", main = "Density")
legend(x = "topright",
       legend = c("true", "Wasserstein", "Grenander"),
       lwd = c(1, 3, 2),
       lty = c(1, 1, 2), col = c("black", "darkgrey", "blue"))

# plot Wasserstein projection density
my_result <- W2_monotone(x, u_seq)
my_density <- monotone_density(my_result, draw = FALSE)
f <- my_density$f
x_seq <- my_density$x
K <- length(my_result$u) - 1
rect(0, 0, x_seq[1], f[1], col = "grey", border = NA)
for (j in 2:K) {
  rect(x_seq[j - 1], 0, x_seq[j], f[j], col = "grey", border = NA)
}
abline(h = 0)

# plot Grenander density
MLE <- Grenander(x)$d
points(c(0, x), MLE(c(0, x)), type = "S",
       col = "blue", lty = 2, lwd = 2)
segments(max(x), MLE(max(x)), max(x), 0,
         col = "blue", lty = 2, lwd = 2)
segments(0, 0, 0, MLE(0), lty = 2, lwd = 2, col = "blue")
# true density
s_seq <- seq(0, 28, by = 0.01)
points(s_seq, d(my_mixture, s_seq), type = "l")
#dev.off()




