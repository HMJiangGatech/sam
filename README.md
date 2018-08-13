# SAM
Sparse Additive Modelling

## Installation

SAM uses OpenMP to enables faster matrix multiplication. So, to use SAM, you must correctly enables OpenMP for the compiler.

For Windows and Linux users, newest version of GCC has fully support of OpenMP.

But for MAC OS users, things a little tricky since the default llvm on MAC OS does not support OpenMP. But the solution is easy. You can simply installs llvm with full OpenMP support and tells R to use this version of llvm.

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

Then you can just install and enable SAM using with the help of CRAN on an R console.

```
install.packages("SAM")
library(SAM)
```

## Experiments

The machine we run experiments on is a PC with

```
RAM: 31.3GB
Processor: Intel Core I7-6700T @ 2.80GHZ x8
Operating System: Ubuntu 16.04 LTS
R version: 3.2.3
GCC version: 5.4.0
```

We compared our results on linear regression and logistic regression with other packages, namely, grplasso, grpreg and gglasso. Also, for SAM and grpreg which support MCP regularization function, we run tests on them both with MCP regularizer and with L1 regularizer.

### Linear Regression

|       SAM with L1 | SAM with MCP | grplasso | grpreg with L1 | grpreg with MCP | gglasso
----- |------------ | -------------| -------- | -------------- | --------------- | -------
time/s| 42.591      | 43.787       | 53.680   | 44.329         | 43.719          | 46.777
loss  | 2.08e-5     | 1.66e-5      | 2.26e-4  | 5.18e-4        | 2.12e-5         | 2.86e-3


### Logistic Regression

|          SAM with L1 | SAM with MCP | grplasso | grpreg with L1 | grpreg with MCP
-------- |------------ | -------------| -------- | -------------- | ---------------
time/s   | 18.700      | 43.787       | 22.229   | 121.701        | 39.9247
accuracy | 0.9054      | 0.9226       | 0.9102   | 0.9052         | 0.8843



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
