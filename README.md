# SAM
Sparse Additive Modelling

## Installation

SAM uses OpenMP to enables faster matrix multiplication. So, to use SAM, you must correctly enables OpenMP for the compiler.

For Windows and Linux users, newest version of GCC has fully support of OpenMP.

But for MAC OS users, things are a little tricky since the default llvm on MAC OS does not support OpenMP. But the solution is easy. You can simply install llvm with full OpenMP support and tells R to use this version of llvm.

```
brew install llvm
```

Then append the following lines into `~/.R/Makevars` to enable llvm with OpenMP support to be the compiler for R packages.

```
CC = /usr/local/bin/clang-omp
CXX = /usr/local/bin/clang-ompclang-ompclang-omp++
CXX98 = /usr/local/bin/clang-omp++
CXX11 = /usr/local/bin/clang-omp++
CXX14 = /usr/local/bin/clang-omp++
CXX17 = /usr/local/bin/clang-omp++
OBJC = /usr/local/bin/clang-omp
OBJCXX = /usr/local/bin/clang-omp++
```

### Installing manually

Download the repo to a folder, open an R terminal in the folder, and then install devtools by typing

```
install.packages(devtools)
```

Then load the devtools package with

```
library(devtools)
```

Then, to install the package, type

```
install()
```

Finally load SAM with

```
library(SAM)
```


### Installing by CRAN (not available yet)

Ideally you can just install and enable SAM using with the help of CRAN on an R console, but we haven't uploaded the package to CRAN yet, so this method is not available so far.

```
install.packages("SAM")
library(SAM)
```

## Usage

With SAM, you can run linear regression, logistic regression and poisson regression.

### Examples for linear regression

```R
## generating training data
n = 100
d = 500
X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)

## generating response
y = -2*sin(X[,1]) + X[,2]^2-1/3 + X[,3]-1/2 + exp(-X[,4])+exp(-1)-1

## Training
out.trn = samQL(X,y)
out.trn

## plotting solution path
plot(out.trn)

## generating testing data
nt = 1000
Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)

yt = -2*sin(Xt[,1]) + Xt[,2]^2-1/3 + Xt[,3]-1/2 + exp(-Xt[,4])+exp(-1)-1

## predicting response
out.tst = predict(out.trn,Xt)
```

### Examples for logistic regression

```R
## generating training data
n = 200
d = 100
X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)
y = sign(((X[,1]-0.5)^2 + (X[,2]-0.5)^2)-0.06)

## flipping about 5 percent of y
y = y*sign(runif(n)-0.05) 
y = sign(y==1)

## Training
out.trn = samLL(X,y)
out.trn

## plotting solution path
plot(out.trn)

## generating testing data
nt = 1000
Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)

yt = sign(((Xt[,1]-0.5)^2 + (Xt[,2]-0.5)^2)-0.06)

## flipping about 5 percent of y
yt = yt*sign(runif(nt)-0.05) 
yt = sign(yt==1)

## predicting response
out.tst = predict(out.trn,Xt)
```

### Examples for Poisson Regression

```R
## generating training data
n = 200
d = 100
X = 0.5*matrix(runif(n*d),n,d) + matrix(rep(0.5*runif(n),d),n,d)
u = exp(-2*sin(X[,1]) + X[,2]^2-1/3 + X[,3]-1/2 + exp(-X[,4])+exp(-1)-1+1)
y = rep(0,n)
for(i in 1:n) y[i] = rpois(1,u[i])

## Training
out.trn = samEL(X,y)
out.trn

## plotting solution path
plot(out.trn)

## generating testing data
nt = 1000
Xt = 0.5*matrix(runif(nt*d),nt,d) + matrix(rep(0.5*runif(nt),d),nt,d)
ut = exp(-2*sin(Xt[,1]) + Xt[,2]^2-1/3 + Xt[,3]-1/2 + exp(-Xt[,4])+exp(-1)-1+1)
yt = rep(0,nt)
for(i in 1:nt) yt[i] = rpois(1,ut[i])

## predicting response
out.tst = predict(out.trn,Xt)
```

To get the complete documentation of SAM, please type `?SAM` in an R terminal.

## Experiments

The experiment scripts are in the folder `tests/testthat/`, to run the experiments, you should open the R terminal in the root of the project folder, and type

```
source("tests/testthat/test_linear.R")
```

or

```
source("tests/testthat/test_logis.R")
```

The machine we ran experiments on is a PC with

```
RAM: 31.3GB
Processor: Intel Core I7-6700T @ 2.80GHZ x8
Operating System: Ubuntu 16.04 LTS
R version: 3.2.3
GCC version: 5.4.0
```

We compared our results on linear regression and logistic regression with other packages, namely, grplasso, grpreg and gglasso. Also, for SAM and grpreg which support MCP regularization function, we run tests on them both with MCP regularizer and with L1 regularizer.

### Linear Regression


|      | SAM with L1 | SAM with MCP | grplasso | grpreg with L1 | grpreg with MCP | gglasso|
| ---- | ----------- | ------------ | -------- | -------------- | --------------- | ------ |
|time/s| 42.591      | 43.787       | 53.680   | 44.329         | 43.719          | 46.777 |
|loss  | 2.08e-5     | 1.66e-5      | 2.26e-4  | 5.18e-4        | 2.12e-5         | 2.86e-3|


### Logistic Regression

|         | SAM with L1 | SAM with MCP | grplasso | grpreg with L1 | grpreg with MCP|
| ------- | ----------- | ------------ | -------- | -------------- | -------------- |
|time/s   | 18.700      | 43.787       | 22.229   | 121.701        | 39.9247        |
|accuracy | 0.9054      | 0.9226       | 0.9102   | 0.9052         | 0.8843         |



## Reference
[1] [Zhao, Tuo, and Han Liu. Sparse Additive Machine. AISTATS. 2012.](http://proceedings.mlr.press/v22/zhao12/zhao12.pdf)  
[2] [Ravikumar, Pradeep, et al. Sparse additive models. 2009](https://rss.onlinelibrary.wiley.com/doi/epdf/10.1111/j.1467-9868.2009.00718.x)   
[3] [Picasso: A Sparse Learning Library for High Dimensional Data Analysis in R and Python](https://cran.r-project.org/web/packages/picasso/vignettes/vignette.pdf)



## GSOC TODO List

### Code Refactoring and Benchmarking
- [x] Use `roxygen2` to manage Documents.
- [x] Benchmark the code.
- [x] Use `Eigen` to accelerate the code.
- [x] Use `OpenMP` to enable multi-processing.

### Solver Upgrade
- [x] Implement Active Set Newton solver.

### Functional Penalties
- [x] Add MCP
- [x] Add SCAD

### Polish and submit
- [x] Polish `readme.md`
- [x] Polish R Documents
- [ ] Polish `vignette`
- [ ] Submit to CRAN
