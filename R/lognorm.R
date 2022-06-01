set.seed(1)
N <- 1000
x <- rnorm(N, 0, 1)
y <- exp(x)
denx <- dnorm(x, 0, 1)
plot(x, denx)
deny <- dlnorm(y, 0, 1)
plot(y, deny)

# estimate f(x) using y

# 1. KDE
r <- ks::hpi(y) # 0.14
ffixed <- ks::kde(y, eval.points = y)$estimate
plot(y, ffixed, cex = .2)

# 2. VKDE with riem
r <- 1
fhat <- rep(0, N)
for(i in 1:N){
  fi <- rep(0, N)
  bindex <- which(abs(y - y[i]) / r <= 1)
  fi[bindex] <- 1 / r * (2 * pi) ^ (-1/2) * exp(- 1 / 2 / r * ((y[bindex] - y[i]) ^ 2)) / (pnorm(1) - pnorm(-1)) * y[bindex]
  fhat <- fhat + fi
}
fhat <- fhat / N
plot(y, fhat, cex = .2)

# Plot together
par(mfrow=c(1,1))
plot(x, fhat, main = paste("Estimated density of y", "r=", round(r,3)), cex = .2, col = "red", lty = 3, xlim = c(min(y), max(y)),
     ylim = c(0, max(denx, ffixed, fhat))
)
points(x, denx, lty = 1, cex = .2, col = 1)
# text(x = 0.5, y = 0.5, paste("f(theta) = ", round(dentheta[1], 3)), col = "red")
points(x, ffixed, col = 3, lty = 2, cex = .2)
legend(x = "topright",          # Position
       legend = c("True density", "Estimates with riemannian", "Estimates with kde"),  # Legend texts
       # lty = c(1, 2, 3),           # Line types
       col = c(1, 2, 3),           # Line colors
       lwd = 2)                 # Line width
