# library(RcppCNPy)
# fmat <- npyLoad("data/hmatrix.npy", dotranspose=FALSE)
# fmat
# imat <- npyLoad("imat.npy", "integer")
# imat
# The fault for this is primarily on the Rcpp data types that are unable to scale above $N$-D array greater than or equal to 4. However, there is no object export inplace for an Rcpp object with 3 dimensions. Hence, there is no $N$-D arrays greater than 2 that can be loaded into R or written to a .npy binary using RcppCNPy.

library(reticulate)
py_config()
# use_condaenv("manifold_env")
# conda_install("numpy")

np <- import("numpy", convert=FALSE)
# data reading
h_isomap <- np$load("data/hmatrix.npy")
h_isomap
class(h_isomap)

# they look quite different but the R and Python arrays are really the same
(y <- py_to_r(h_isomap))
r_to_py(y)
