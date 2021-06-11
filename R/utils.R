library(reticulate)
np = import("numpy")
# adj_matrix <- load("data/adjacency_matrix_1id336tow.npy")
adj_matrix <- feather::read_feather("data/adjacency_matrix_1id336tow.feather")
# adj_matrix %>% view()
class(adj_matrix)
colnames(adj_matrix) <- paste0("p", 1:336)
# rownames(adj_matrix) <- paste0("p", 1:336)
# adj_matrix <- as.data.frame(adj_matrix)

g <- igraph::graph_from_data_frame(adj_matrix, directed = TRUE)
plot(g)
is_connected(g) # FALSE?
a <- as_adjacency_matrix(g)
dim(a)
head(a)

graph.data.frame(adj_matrix, directed = FALSE)

get.adjacency(g)

# contour plots
filled.contour(x = 10*1:nrow(volcano),y = 10*1:ncol(volcano),
               z = volcano, color.palette = terrain.colors,
               plot.title = title(main = "The Topography of Maunga Whau",
                                  xlab = "Meters North",ylab = "Meters West"),
               plot.axes = {axis(1, seq(100, 800, by = 100))
                 axis(2, seq(100, 600, by = 100))},
               key.title = title(main="Height\n(meters)"),
               key.axes = axis(4, seq(90, 190, by = 10)))


embed_den <- as_tibble(cbind(x = x[,1], y = x[,2], z = f)) %>%
  drop_na() # %>% 
# arrange(x, y)

embed_den_list <- list(x = embed_den$x, y = embed_den$y, z = embed_den$z)

embed_den
hdrscatterplot(embed_den$x, embed_den$y)
hdrinfo <- hdr.2d(embed_den$x, embed_den$y, prob = c(50, 95, 99)) 
plot.hdr2d(hdrinfo)
# If density values are specified
hdrinfo1 <- hdr.2d(embed_den$x, embed_den$y, prob = c(50, 95, 99), den = embed_den_list)
plot.hdr2d(hdrinfo)
# filled.contour(x = embed_den$x, y = embed_den$y, z = embed_den$z)

# interpolation
library(akima)
?interp()
akima.spl <- interp(embed_den$x, embed_den$y, embed_den$z, nx=100, ny=100, linear=FALSE)
p.zmin <- min(embed_den$z,na.rm=TRUE)
p.zmax <- max(embed_den$z,na.rm=TRUE)
breaks <- pretty(c(p.zmin,p.zmax),10)
colors <- heat.colors(length(breaks)-1)
hdrscatterplot(embed_den$x, embed_den$y)
plot.hdr2d(hdrinfo)
contour(akima.spl, main = "smooth interp(*, linear = FALSE)", col = colors, levels=breaks, add=TRUE)
points(embed_den, pch = 20)

filled.contour(akima.spl, color.palette = viridis,
               plot.axes = { axis(1); axis(2);
                 title("smooth  interp(*, linear = FALSE)");
                 points(embed_den, pch = 3, col= hcl(c=20, l = 10))})



# approx and approxfun
x <- 1:10
y <- rnorm(10)
par(mfrow = c(2,1))
plot(x, y, main = "approx(.) and approxfun(.)")
points(approx(x, y), col = 2, pch = "*")
points(approx(x, y, method = "constant"), col = 4, pch = "*")

f <- approxfun(x, y)
curve(f(x), 0, 11, col = "green2")
points(x, y)
is.function(fc <- approxfun(x, y, method = "const")) # TRUE
curve(fc(x), 0, 10, col = "darkblue", add = TRUE)
## different extrapolation on left and right side :
plot(approxfun(x, y, rule = 2:1), 0, 11,
     col = "tomato", add = TRUE, lty = 3, lwd = 2)






# ellipse plot
# different packages to produce ellipse plot
library(psych)
data(sat.act)
ellipses(sat.act)  #shows the pairs.panels ellipses
minkowski(2,main="Minkowski circles")
minkowski(1,TRUE)
minkowski(4,TRUE)

### The logic to plot an ellipse based on data centroid and covariance matrix
ctr    <- c(0, 0)                               # data centroid -> colMeans(dataMatrix)
A      <- matrix(c(2.2, 0.4, 0.4, 2.8), nrow=2) # covariance matrix -> cov(dataMatrix)
RR     <- chol(A)                               # Cholesky decomposition
angles <- seq(0, 2*pi, length.out=200)          # angles for ellipse
ell    <- 1 * cbind(cos(angles), sin(angles)) %*% RR  # ellipse scaled with factor 1
ellCtr <- sweep(ell, 2, ctr, "+")               # center ellipse to the data centroid
plot(ellCtr, type="l", lwd=2, asp=1)            # plot ellipse
points(ctr[1], ctr[2], pch=4, lwd=2)            # plot data centroid

library(car)  # verify with car's ellipse() function
ellipse(c(0, 0), shape=A, radius=1, col="red", lty=2, draw=T, add=F, fill = T)

# embed_den
# x = metric_isomap$embedding
# plot(x[,1], x[,2])
# ellipse(c(x[1,1], x[1,2]), shape = riem_isomap[,,1], radius=5, col="red", lty=2, fill = FALSE, fill.alpha = 0.1)
# riemannian matrix not positive definite: it has negative eigenvalues and isn't a possible covariance matrix


library(ellipse)
plot(ellipse::ellipse(0), type = "l")

for(i in 1:N){
  mat <- riem_isomap[,,i]
  center <- c(x[i,1], x[i,2])
  add <- ifelse(i == 1, F, T)
  car::ellipse(center, mat, radius = .002, 
               col = blues9[5], asp = 1, pty = "s", lwd = 1, center.pch = 19, center.cex = 0,
               fill = T, fill.alpha = 0.2, add = add, grid = T,
               xlim = c(-.2, .25), ylim = c(-.2, .2))
}

