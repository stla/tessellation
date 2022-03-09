
cpp <- '
arma::umat tcp(const arma::umat & A){
  arma::umat M(A*A.t());
  return M;
}
'
library(Rcpp)
cppFunction(
  cpp,
  depends = "RcppArmadillo"
)

A <- toeplitz(c(1L, 2L))
M <- tcp(A)
str(M)
# num [1:2, 1:2] 5 4 4 5


cpp <- '
arma::mat fitPlane(const arma::mat& points){
  int K(points.n_rows);
  arma::mat A(points.t());
  arma::vec3 p3 = (1.0/static_cast<double>(K))*arma::sum(A,1);
  A.each_col() -= p3; //[x1-p, x2-p, ..., xk-p]
  arma::mat33 M(A*A.t());
  arma::vec3 D;
  arma::mat33 V;
  if(arma::eig_sym(D,V,M)){
    // diagonalization succeeded
    arma::vec3 n3 = V.col(0); // in ascending order by default
    if(n3(2) < 0){
      n3 = -n3; // upward pointing
    }
    return arma::join_rows(n3, p3);
  }else{
    throw Rcpp::exception("Something bad occured.");
  }
}'

library(Rcpp)
cppFunction(
  cpp,
  depends = "RcppArmadillo"
)

set.seed(666)
points <- matrix(rgamma(30L, 10, 1), nrow = 10L, ncol = 3L)
fitPlane(points)

a <- 1
b <- 1
c <- 1
d <- 4 #* sqrt(c(crossprod(c(a, b, c))))
x <- rgamma(30L, 10, 1)
y <- rgamma(30L, 10, 1)
z <- d - (a*x + b*y) / c #+ rnorm(30L, 0, 0.1)
points <- cbind(x, y, z)
p <- fitPlane(points)
n <- p[, 1L]
M <- p[, 2L]
a <- n[1L]
b <- n[2L]
c <- n[3L]
d <- crossprod(n, M)
O <- c(0, 0, d/c)
OM <- (M-O) / sqrt(c(crossprod(M-O)))
OQ <- tessellation:::crossProduct(OM, n)

coordinates <- matrix(NA_real_, nrow = nrow(points), ncol = 2L)
for(i in 1L:nrow(points)){
  pt <- points[i, ]
  Opt <- pt - O
  coordinates[i, ] <- c(
    crossprod(OM, Opt),
    crossprod(OQ, Opt)
  )
}



##
set.seed(666)
points <- matrix(rgamma(30L, 10, 1), nrow = 10L, ncol = 3L)

f <- function(x, y){
  exp(-(x*x + y*y))
}
x <- seq(-6, 6, length.out = 100)
y <- seq(-6, 6, length.out = 100)
points <- as.matrix(transform( # data (x_i, y_i, z_i)
  expand.grid(x = x, y = y), z = f(x, y)
))

points0 <- points
points <- sweep(points, 2, colMeans(points))
p <- fitPlane(points)
n <- p[, 1L]
#M <- colMeans(points)
M <- c(0, 1, 0)
a <- n[1L]
b <- n[2L]
c <- n[3L]
d <- crossprod(n, M)
O <- c(0, 0, 0)
OM <- (M-O) / sqrt(c(crossprod(M-O)))
OQ <- tessellation:::crossProduct(n,OM)

coordinates <- matrix(NA_real_, nrow = nrow(points), ncol = 2L)
for(i in 1L:nrow(points)){
  pt <- points[i, ]
  Opt <- pt - O
  coordinates[i, ] <- c(
    crossprod(OM, Opt),
    crossprod(OQ, Opt)
  )
}

del <- delaunay(coordinates)
triangles <- do.call(rbind, lapply(del[["tiles"]], function(tile){
  indices <- tile[["vertices"]]
  if(tile[["orientation"]] == -1L){
    indices <- indices[c(1L, 3L, 2L)]
  }
  indices
}))
volumes <- apply(triangles, 1L, function(trgl){
  trgl <- points0[trgl, ]
  tessellation:::volume_under_triangle(trgl[, 1L], trgl[, 2L], trgl[, 3L])
})
sum(volumes)


#############
fit <- lm(points[, 3L] ~ points[, 1L] + points[, 2L])

x <- points[, 1L]
y <- points[, 2L]
z <- points[, 3L]
A <- rbind(
  c(crossprod(x), crossprod(x, y), sum(x)),
  c(crossprod(x,y), crossprod(y), sum(y)),
  c(sum(x), sum(y), nrow(points))
)
b <- c(crossprod(x,z), crossprod(y,z), sum(z))
