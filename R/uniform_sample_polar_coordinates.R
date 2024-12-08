# How to sample uniformly from 3-D manifolds (mapping function) and transform the density after change of coordinates
library(tidyverse)
library(scatterplot3d)
library(plotly)
source("R/sources/mldata.R")
## --------------------------------------------------------
# 1. Sample uniformly from the unit disk
## --------------------------------------------------------
# The unit disk is given by the polar coordinates in the set [0,1] x [0,2pi)
# How to sample the polar coordinates so that the resulting point distribution is uniform on the disk
# Reference: Section 3.3 of https://www.cs.cornell.edu/courses/cs6630/2015fa/notes/pdf-transform.pdf
# (r, phi) -> (x, y)
# x = r*cos(phi); y = r*sin(phi)
# r := sqrt(e1); phi = 2 * pi * e2; where e1, e2 ~ U(0,1)
set.seed(1234)
N <- 2000
e1 <- runif(N, 0, 1)
e2 <- runif(N, 0, 1)
r <- sqrt(e1)
phi <- 2 * pi * e2
# e1, e2 ~ U(0,1)
# density of r: p_r(r) = p_{e1}(e1) * (2*r) = 2*r
# density of phi: p_phi(phi) = p_{e2}(e2) * 1 / (2 * pi) = 1/(2 * pi)
# CDF of r: P_r(r) = \int_0^r  p_r(x) dx = r^2
# CDF of phi: P_phi(phi) = \int_0^phi (2*pi) dx = 2*pi*phi
# Density of (r, phi): p_A(r, phi) = p_r(r) * p_phi(phi) = 2 * r * 1 / (2 * pi) = r / pi
# Density of (x, y): p_B(x, y) = p_A(r, phi) / |J| = (r / pi) / r = 1 / pi
# Therefore, (x, y) is uniformly distributed with density 1 / pi.
x <- r * cos(phi)
y <- r * sin(phi)
plot(x, y, asp = 1)

# True density of (x,y) is known as constant, 1 / pi
# Test if (x, y) is uniformly distributed in the set [0, 1] x [0, 1] 
# ks.test(cbind(x,y), "punif", 0, 1)


## --------------------------------------------------------
# 2. Sample uniformly from the sphere S^2
## --------------------------------------------------------
# (theta, phi) \in (0, pi) x [0, 2pi)
# P_B(w(theta, phi)) = P_A(theta, phi) / sin(theta)
# phi := 2*pi*e2; cos(theta) := 1 - 2*e1
# x:=sqrt(1 - cos(theta)^2) * cos(phi); y:=sqrt(1 - cos(theta)^2) * sin(phi); z:=cos(theta)
z <- 1 - 2 * e1
x <- sqrt(1 - z ^ 2) * cos(phi)
y <- sqrt(1 - z ^ 2) * sin(phi)
scatterplot3d(x, y, z,
              color = "blue",
              asp = 1,
              pch = 20)


## --------------------------------------------------------
# 3. Sample from twin peaks
## --------------------------------------------------------
set.seed(1)
mapping <- c("Swiss Roll", "semi-sphere", "Twin Peak", "S Curve")[3] # only run [1] and [3]
sr <- mldata(N = N, meta = "gaussian", mapping = mapping)
# sr <- mldata(N = N, meta = "uniform", mapping = mapping)
swissroll <- sr$data %>% as.data.frame()
preswissroll <- sr$metadata %>% as_tibble() %>% mutate(den = sr$den, label = c(rep(1:4, each = N / 4)))
colnames(preswissroll) <- c("X1", "X2", "den", "label")
plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = sr$den, # colored with 2d density
        type = "scatter3d", mode = "markers", size = 1, text = paste("density:", preswissroll$den))
summary(sr$den) # [0, 1]

# Density of 3d twin peaks using change of coordinates
# P_B(w(u,v)) = P_A(u,v) / sqrt(EG-F^2), where
# EG-F^2 = 1 + (pi * cos(pi*u) * tan(3*v)) ^ 2 + 9 * sin(pi * u) ^ 2 / (cos(3 * v)) ^ 4
# P_A(u,v) = sr$den
u <- preswissroll$X1
v <- preswissroll$X2
den_twinpeaks <- sr$den / sqrt(1 + (pi * cos(pi * u) * tan(3 * v)) ^ 2 + 9 * sin(pi * u) ^ 2 / (cos(3 * v)) ^ 4)
summary(den_twinpeaks)
plot_ly(data = swissroll, x = ~ x, y = ~ y, z = ~ z, color = den_twinpeaks, # colored with 3d density
        type = "scatter3d", mode = "markers", size = 1, text = paste("density:", preswissroll$den))

