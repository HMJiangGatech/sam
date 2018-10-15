set.seed(19970224)
nt = 10000
n = 10000
d = 2000
p = 3
m = d * p
t = 3

library(splines)

X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)
y = sign(((X[,1]-0.5)^2 + (X[,2]-0.5)^2)-0.06)
## flipping about 5 percent of y
y = y*sign(runif(n)-0.05) 
y = sign(y==1)

genZ_t = proc.time()
Z = matrix(0,n,m)
for(j in 1:d){
  tmp = (j-1)*p + c(1:p)
  tmp0 = ns(X[,j],df=p)
  Z[,tmp] = tmp0
}

Z <-cbind(rep(1, n), Z)
colnames(Z) <- c("Intercept", paste("X", 1:m, sep = ""))

genZ_t = proc.time() - genZ_t

Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)
yt = sign(((Xt[,1]-0.5)^2 + (Xt[,2]-0.5)^2)-0.06)
## flipping about 5 percent of y
yt = yt*sign(runif(nt)-0.05) 
yt = sign(yt==1)

Zt = matrix(0,nt,m)
for(j in 1:d){
  tmp = (j-1)*p + c(1:p)
  tmp0 = ns(Xt[,j],df=p)
  Zt[,tmp] = tmp0
}

Zt <-cbind(rep(1, n), Zt)
colnames(Zt) <- c("Intercept", paste("X", 1:m, sep = ""))

index <- NA
for (i in 1:d) {
  index <- c(index, rep(i, p))
}

source("tests/testthat/logis_sam.R")
#source("tests/testthat/logis_gglasso.R")
source("tests/testthat/logis_grplasso.R")
source("tests/testthat/logis_grpreg.R")


