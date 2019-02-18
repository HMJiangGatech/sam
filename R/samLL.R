#-----------------------------------------------------------------------#
# Package: SAM                                                          #
# Method: Sparse Additive Modelling using Logistic Loss                 #
#-----------------------------------------------------------------------#

#' Training function of Sparse Additive Logistic Regression
#'
#' The logistic model is learned using training data.
#'
#' We adopt various computational algorithms including the block coordinate descent, fast iterative soft-thresholding algorithm, and newton method. The computation is further accelerated by "warm-start" and "active-set" tricks.
#'
#' @param X The \code{n} by \code{d} design matrix of the training set, where \code{n} is sample size and \code{d} is dimension.
#' @param y The \code{n}-dimensional label vector of the training set, where \code{n} is sample size. Labels must be coded in 1 and 0.
#' @param p The number of basis spline functions. The default value is 3.
#' @param lambda A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio. Supplying a value of lambda overrides this. WARNING: use with care. Do not supply a single value for lambda. Supply instead a decreasing sequence of lambda values. samLL relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.
#' @param nlambda The number of lambda values. The default value is 20.
#' @param lambda.min.ratio Smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default is 0.1.
#' @param thol Stopping precision. The default value is 1e-5.
#' @param max.ite The number of maximum iterations. The default value is 1e5.
#' @param regfunc A string indicating the regularizer. The default value is "L1". You can also assign "MCP" or "SCAD" to it.
#' @return
#' \item{p}{
#'   The number of basis spline functions used in training.
#' }
#' \item{X.min}{
#'   A vector with each entry corresponding to the minimum of each input variable. (Used for rescaling in testing)
#' }
#' \item{X.ran}{
#'   A vector with each entry corresponding to the range of each input variable. (Used for rescaling in testing)
#' }
#' \item{lambda}{
#'   A sequence of regularization parameter used in training.
#' }
#' \item{w}{
#'   The solution path matrix (\code{d*p+1} by length of \code{lambda}) with each column corresponding to a regularization parameter. Since we use the basis expansion with the intercept, the length of each column is \code{d*p+1}.
#' }
#' \item{df}{
#'   The degree of freedom of the solution path (The number of non-zero component function)
#' }
#' \item{knots}{
#'   The \code{p-1} by \code{d} matrix. Each column contains the knots applied to the corresponding variable.
#' }
#' \item{Boundary.knots}{
#'   The \code{2} by \code{d} matrix. Each column contains the boundary points applied to the corresponding variable.
#' }
#' \item{func_norm}{
#'   The functional norm matrix (\code{d} by length of \code{lambda}) with each column corresponds to a regularization parameter. Since we have \code{d} input variables, the length of each column is \code{d}.
#' }
#' @seealso \code{\link{SAM}},\code{\link{plot.samLL},\link{print.samLL},\link{predict.samLL}}
#' @examples
#'
#' ## generating training data
#' n = 200
#' d = 100
#' X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)
#' y = sign(((X[,1]-0.5)^2 + (X[,2]-0.5)^2)-0.06)
#'
#' ## flipping about 5 percent of y
#' y = y*sign(runif(n)-0.05)
#' y = sign(y==1)
#'
#' ## Training
#' out.trn = samLL(X,y)
#' out.trn
#'
#' ## plotting solution path
#' plot(out.trn)
#'
#' ## generating testing data
#' nt = 1000
#' Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)
#'
#' yt = sign(((Xt[,1]-0.5)^2 + (Xt[,2]-0.5)^2)-0.06)
#'
#' ## flipping about 5 percent of y
#' yt = yt*sign(runif(nt)-0.05)
#' yt = sign(yt==1)
#'
#' ## predicting response
#' out.tst = predict(out.trn,Xt)
#' @useDynLib SAM grpLR
#' @export
samLL = function(X, y, p=3, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.1, thol=1e-5, max.ite = 1e5, regfunc="L1"){

  gcinfo(FALSE)
  fit = list()
  fit$p = p

  fit = list()

  X = as.matrix(X)
  y = as.vector(y)

  n = nrow(X)
  d = ncol(X)
  m = d*p

  n1 = sum(y==1)
  n0 = sum(y==0)

  if((n1+n0)!=n){
    cat("Please check the labels. (Must be coded in 1 and 0)")
    fit = "Please check the labels."
    return(fit)
  }

  fit$p = p

  X.min = apply(X,2,min)
  X.max = apply(X,2,max)
  X.ran = X.max - X.min
  X.min.rep = matrix(rep(X.min,n),nrow=n,byrow=T)
  X.ran.rep = matrix(rep(X.ran,n),nrow=n,byrow=T)
  X = (X-X.min.rep)/X.ran.rep

  fit$X.min = X.min
  fit$X.ran = X.ran

  Z = matrix(0,n,m)
  fit$nkots = matrix(0,p-1,d)
  fit$Boundary.knots = matrix(0,2,d)
  for(j in 1:d){
    tmp = (j-1)*p + c(1:p)
    tmp0 = ns(X[,j],df=p)
    Z[,tmp] = tmp0
    fit$nkots[,j] = attr(tmp0,'knots')
    fit$Boundary.knots[,j] = attr(tmp0,'Boundary.knots')
  }

  L0 = norm(Z,"f")^2

  #if(n>d) L0 = max(eigen(Z%*%t(Z), symmetric=T, only.values = T)$values) else L0 = max(eigen(t(Z)%*%Z, symmetric=T, only.values = T)$values)


  Z = cbind(Z,rep(1,n))

  a0 = -log(n/n1-1)

  z = colSums(matrix(rep(y,m+1),n,m+1)*Z)

  if(is.null(lambda)){
    g = -z + colSums(matrix(rep(n1/n,m+1),n,m+1)*Z)

    if(is.null(nlambda)) nlambda = 20;

    lambda_max=max(sqrt(colSums(matrix(g[1:(p*d)],p,d)^2)))
    lambda = exp(seq(log(1),log(lambda.min.ratio),length=nlambda))*lambda_max
  } else nlambda = length(lambda)

  out = .C("grpLR", A = as.double(Z), y = as.double(y), lambda = as.double(lambda), nlambda = as.integer(nlambda), LL0 = as.double(L0), nn = as.integer(n), dd = as.integer(d), pp = as.integer(p), xx = as.double(matrix(0,m+1,nlambda)), aa0 = as.double(a0), mmax_ite = as.integer(max.ite), tthol = as.double(thol), regfunc = as.character(regfunc), aalpha = as.double(0.5), z = as.double(z),df = as.integer(rep(0,nlambda)),func_norm = as.double(matrix(0,d,nlambda)), package="SAM")

  fit$lambda = out$lambda
  fit$w = matrix(out$xx,ncol=nlambda)
  fit$df = out$df
  fit$func_norm = matrix(out$func_norm,ncol=nlambda)

  rm(out,X,y,Z,X.min.rep,X.ran.rep)

  class(fit) = "samLL"
  return(fit)
}

#' Printing function for S3 class \code{"samLL"}
#'
#' Summarize the information of the object with S3 class \code{samLL}.
#'
#' The output includes length and d.f. of the regularization path.
#'
#' @param x An object with S3 class \code{"samLL"}
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{samLL}}
#' @export
print.samLL = function(x,...){
  cat("Path length:",length(x$df),"\n")
  cat("d.f.:",x$df[1],"--->",x$df[length(x$df)],"\n")
}


#' Plot function for S3 class \code{"samLL"}
#'
#' This function plots the regularization path (regularization parameter versus functional norm)
#'
#' The horizontal axis is for the regularization parameters in log scale. The vertical axis is for the functional norm of each component.
#'
#' @param x An object with S3 class \code{"samLL"}
#' @param \dots System reserved (No specific usage)
#' @seealso \code{\link{samLL}}
#' @export
plot.samLL = function(x,...){
  par = par(omi = c(0.0, 0.0, 0, 0), mai = c(1, 1, 0.1, 0.1))
  matplot(x$lambda[length(x$lambda):1],t(x$func_norm),type="l",xlab="Regularization Parameters",ylab = "Functional Norms",cex.lab=2,log="x",lwd=2)
}

#' Prediction function for S3 class \code{"samLL"}
#'
#' Predict the labels for testing data.
#'
#' The testing dataset is rescale to the samLLe range, and expanded by the samLLe spline basis functions as the training data.
#'
#' @param object An object with S3 class \code{"samLL"}.
#' @param newdata The testing dataset represented in a \code{n} by \code{d} matrix, where \code{n} is testing sample size and \code{d} is dimension.
#' @param thol The decision value threshold for prediction. The default value is 0.5
#' @param \dots System reserved (No specific usage)
#' @return
#' \item{probs}{
#'   Estimated Posterior Probability for Prediction also represented in a \code{n} by the length of \code{lambda} matrix, where \code{n} is testing sample size.
#' }
#' \item{labels}{
#'   Predicted labels also represented in a \code{n} by the length of \code{lambda} matrix, where \code{n} is testing sample size. }
#' @seealso \code{\link{samLL}}
#' @export
predict.samLL = function(object, newdata, thol = 0.5 ,...){
  gcinfo(FALSE)
  out = list()
  nt = nrow(newdata)
  d = ncol(newdata)
  X.min.rep = matrix(rep(object$X.min,nt),nrow=nt,byrow=T)
  X.ran.rep = matrix(rep(object$X.ran,nt),nrow=nt,byrow=T)

  newdata = (newdata-X.min.rep)/X.ran.rep
  newdata = pmax(newdata,0)
  newdata = pmin(newdata,1)

  m = object$p*d

  Zt = matrix(0,nt,m)

  for(j in 1:d){
    tmp = (j-1)*object$p + c(1:object$p)
    Zt[,tmp] = ns(newdata[,j],df=object$p,knots=object$knots[,j],Boundary.knots=object$Boundary.knots[,j])
  }

  out$probs = exp(cbind(Zt,rep(1,nt))%*%object$w)
  out$probs = out$probs/(1+out$probs)
  out$labels = sign(out$probs>thol)

  rm(Zt,newdata)

  return(out)
}
