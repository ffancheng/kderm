# neighborhood graph
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

# contour plots
filled.contour(x = 10*1:nrow(volcano),y = 10*1:ncol(volcano),
               z = volcano, color.palette = terrain.colors,
               plot.title = title(main = "The Topography of Maunga Whau",
                                  xlab = "Meters North",ylab = "Meters West"),
               plot.axes = {axis(1, seq(100, 800, by = 100))
                 axis(2, seq(100, 600, by = 100))},
               key.title = title(main="Height\n(meters)"),
               key.axes = axis(4, seq(90, 190, by = 10)))



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




## Turn x and y into coordinates, z as cell values
# coordinates into grids

embed_den <- as_tibble(cbind(x = x[,1], y = x[,2], z = f)) %>%
  drop_na() # %>% 
data.df <- embed_den

ji <- function(xy, origin=c(0,0), cellsize=c(1,1)) {
  t(apply(xy, 1, function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
}
JI <- ji(cbind(data.df$x, data.df$y))
data.df$X <- JI[, 1]
data.df$Y <- JI[, 2]
data.df$Cell <- paste(data.df$X, data.df$Y)

counts <- by(data.df, data.df$Cell, function(d) c(d$X[1], d$Y[1], nrow(d)))
counts.m <- matrix(unlist(counts), nrow=3)
rownames(counts.m) <- c("X", "Y", "Count")
# write.csv(as.data.frame(t(counts.m)), "f:/temp/grid.csv")

count.max <- max(counts.m["Count",])
colors = sapply(counts.m["Count",], function(n) hsv(sqrt(n/count.max), .7, .7, .5))
plot(counts.m["X",] + 1/2, counts.m["Y",] + 1/2, cex=sqrt(counts.m["Count",]/100),
     pch = 19, col=colors,
     xlab="Longitude of cell center", ylab="Latitude of cell center",
     main="Event counts within one-degree grid cells")




library(sp)
library(mapview)
station <- data.frame(lat = c(41.997946, 41.960669, 41.960669, 41.960669,41.909269,41.931841,41.909269,41.910561,41.866129,41.866129), long = c(-87.654561, -87.747456, -87.67459, -87.646438,-87.747456,-87.67459,-87.67459,-87.619112,-87.747456,-87.691617),station = 1:10)
coordinates(station) = ~ long+lat
proj4string(station) <- CRS("+proj=longlat +datum=WGS84")
stp <- spTransform(station, CRSobj = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m"))

mapview(stp)

library(raster)
r <- raster(stp, res=250)
rr <- setExtent(r, round(extent(r)+10000,-3), keepres=TRUE)
plot(rr)

data.df <- embed_den
coordinates(data.df) = ~x+y
proj4string(data.df) <- CRS("+proj=longlat +datum=WGS84")
stp <- spTransform(data.df, CRSobj = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m"))
mapview(stp)


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

# Ellipse covariance matrix to semi-major and semi-minor axis length a, b
x <- c(1.798805,2.402390,2.000000,3.000000,1.000000)
y <- c(0.3130147,0.4739707,0.2000000,0.8000000,0.1000000)
d <- cbind( x, y )
library(cluster)
r <- ellipsoidhull(d)
plot( x, y, asp=1, xlim=c(0,4) )
lines( predict(r) )
e <- sqrt(eigen(r$cov)$values)
a <- sqrt(r$d2) * e[1]  # semi-major axis
b <- sqrt(r$d2) * e[2]  # semi-minor axis

theta <- seq(0, 2*pi, length=200)
lines( r$loc[1] + a * cos(theta), r$loc[2] + a * sin(theta) )
lines( r$loc[1] + b * cos(theta), r$loc[2] + b * sin(theta) )

# Works with for loop
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



# Convert covariance matrix to a, b, angle, radius
library(tibble)
library(MASS)
sigma=matrix(c(1.4,-0.5,-.5,0.4), nrow = 2, ncol = 2)
samples<-mvrnorm(n=100, c(1,2), sigma)
colnames(samples)<-c("X","Y")
samples_df = as_tibble(samples)

cov2elipse<-function(my_cov, alpha=0.05) {
  aa = my_cov[1,1]
  bb = my_cov[1,2]
  cc = my_cov[2,2]
  a = (aa+cc)/2 + sqrt(((aa-cc)/2)^2 + bb^2) 
  b = (aa+cc)/2 - sqrt(((aa-cc)/2)^2 + bb^2) 
  # angle = atan(bb/(a-aa))
  angle = atan((a-aa)/bb)
  r1=sqrt(qchisq(1 - alpha, 2))
  return(list(a=sqrt(a)*r1, b=sqrt(b)*r1, angle=angle))
}

my_cov<-cov(samples_df)
par_elipsy<-cov2elipse(my_cov)

ggplot(samples_df, aes(x=X, y=Y)) +
  geom_point() + 
  geom_ellipse(mapping = aes(x0=1, y0=2, a=par_elipsy$a, b=par_elipsy$b, angle=par_elipsy$angle, fill="red", alpha=0.01))  






# Rtsne
# https://github.com/jkrijthe/Rtsne/pull/39

set.seed(101)
library(Rtsne)
iris_matrix <- matrix(rnorm(150*4), nrow=150) # See note 1

library(FNN)
NN <- get.knn(iris_matrix, 90)
D <- as.matrix(dist(iris_matrix))
re.nn <- NN$nn.dist
for (i in seq_len(ncol(re.nn))) {
  re.nn[,i] <- D[cbind(seq_len(nrow(re.nn)), NN$nn.index[,i])] 
}
range(re.nn - NN$nn.dist) # check it gives 0 0; no numerical precision problems here.

Y_in <- matrix(runif(nrow(iris_matrix)*2), ncol=2) # See note 2

par(mfrow=c(1,2))
set.seed(42) 
out <- Rtsne:::Rtsne_nn_cpp(t(NN$nn.index - 1L), t(NN$nn.dist), 
                            # origD=ncol(iris_matrix), 
                            no_dims=2,
                            Y_in=Y_in, 
                            init=TRUE, perplexity = 30, theta = 0.5, max_iter = 1000,
                            verbose = TRUE,
                            stop_lying_iter = 250L, mom_switch_iter = 250L, 
                            momentum = 0.5, final_momentum = 0.8, eta = 200, exaggeration_factor = 12, num_threads = 1)
plot(t(out$Y))

set.seed(42) 
blah <- Rtsne(D, is_distance=TRUE, verbose=TRUE, max_iter=1000,
              stop_lying_iter=250L, mom_switch_iter=250L, Y_init=Y_in) # See note 3
plot(blah$Y)
