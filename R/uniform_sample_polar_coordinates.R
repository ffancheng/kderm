# sample uniformly from the unit disk
# The unit disk is given by the polar coordinates in the set [0,1] x [0,2pi)
# How to sample the polar coordinates so that the resulting point distribution is uniform on the disk
# Reference: Section 3.3 of https://www.cs.cornell.edu/courses/cs6630/2015fa/notes/pdf-transform.pdf
# (r, phi) -> (x, y)
# x = r*cos(phi); y = r*sin(phi)
# r := sqrt(e1); phi = 2 * pi * e2; where e1, e2 ~ U(0,1)
set.seed(1234)
N <- 1000
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
