% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samHL.R
\name{samHL}
\alias{samHL}
\title{Training function of Sparse Additive Machine}
\usage{
samHL(
  X,
  y,
  p = 3,
  lambda = NULL,
  nlambda = NULL,
  lambda.min.ratio = 0.4,
  thol = 1e-05,
  mu = 0.05,
  max.ite = 1e+05,
  w = NULL
)
}
\arguments{
\item{X}{The \code{n} by \code{d} design matrix of the training set, where \code{n} is sample size and \code{d} is dimension.}

\item{y}{The \code{n}-dimensional label vector of the training set, where \code{n} is sample size. Labels must be coded in 1 and 0.}

\item{p}{The number of basis spline functions. The default value is 3.}

\item{lambda}{A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio. Supplying a value of lambda overrides this. WARNING: use with care. Do not supply a single value for lambda. Supply instead a decreasing sequence of lambda values. samHL relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.}

\item{nlambda}{The number of lambda values. The default value is 20.}

\item{lambda.min.ratio}{Smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default is 0.4.}

\item{thol}{Stopping precision. The default value is 1e-5.}

\item{mu}{Smoothing parameter used in approximate the Hinge Loss. The default value is 0.05.}

\item{max.ite}{The number of maximum iterations. The default value is 1e5.}

\item{w}{The \code{n}-dimensional positive vector. It is the weight of each entry in the weighted loss. The default value is 1 for all entries.}
}
\value{
\item{p}{
  The number of basis spline functions used in training.
}
\item{X.min}{
  A vector with each entry corresponding to the minimum of each input variable. (Used for rescaling in testing)
}
\item{X.ran}{
  A vector with each entry corresponding to the range of each input variable. (Used for rescaling in testing)
}
\item{lambda}{
  A sequence of regularization parameter used in training.
}
\item{w}{
  The solution path matrix (\code{d*p+1} by length of \code{lambda}) with each column corresponding to a regularization parameter. Since we use the basis expansion with the intercept, the length of each column is \code{d*p+1}.
}
\item{df}{
  The degree of freedom of the solution path (The number of non-zero component function)
}
\item{knots}{
  The \code{p-1} by \code{d} matrix. Each column contains the knots applied to the corresponding variable.
}
\item{Boundary.knots}{
  The \code{2} by \code{d} matrix. Each column contains the boundary points applied to the corresponding variable.
}
\item{func_norm}{
  The functional norm matrix (\code{d} by length of \code{lambda}) with each column corresponds to a regularization parameter. Since we have \code{d} input variables, the length of each column is \code{d}.
}
}
\description{
The classifier is learned using training data.
}
\details{
We adopt various computational algorithms including the block coordinate descent, fast iterative soft-thresholding algorithm, and newton method. The computation is further accelerated by "warm-start" and "active-set" tricks.
}
\examples{

## generating training data
n = 200
d = 100
X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)
y = sign(((X[,1]-0.5)^2 + (X[,2]-0.5)^2)-0.06)

## flipping about 5 percent of y
y = y*sign(runif(n)-0.05)

## Training
out.trn = samHL(X,y)
out.trn

## plotting solution path
plot(out.trn)

## generating testing data
nt = 1000
Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)

yt = sign(((Xt[,1]-0.5)^2 + (Xt[,2]-0.5)^2)-0.06)

## flipping about 5 percent of y
yt = yt*sign(runif(nt)-0.05)

## predicting response
out.tst = predict(out.trn,Xt)
}
\seealso{
\code{\link{SAM}},\code{\link{plot.samHL},\link{print.samHL},\link{predict.samHL}}
}
