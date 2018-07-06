# SAM
Sparse Additive Modelling

## Installation

## Reference
[1] [Zhao, Tuo, and Han Liu. Sparse Additive Machine. AISTATS. 2012.](http://proceedings.mlr.press/v22/zhao12/zhao12.pdf)  
[2] [Ravikumar, Pradeep, et al. Sparse additive models. 2009](https://rss.onlinelibrary.wiley.com/doi/epdf/10.1111/j.1467-9868.2009.00718.x)   
[3] [Picasso: A Sparse Learning Library for High Dimensional Data Analysis in R and Python](https://cran.r-project.org/web/packages/picasso/vignettes/vignette.pdf)



## GSOC TODO List

### Code Refactoring and Benchmarking
- [x] Use `roxygen2` to manage Documents.
- [ ] Benchmark the code.
- [ ] Use `RcppEigen` to accelerate the code.
- [ ] Use `OpenMP` to enable multi-processing.

### Solver Upgrade
- [ ] Implement Active Set Newton solver.

### Functional Penalties
- [ ] Add MCP
- [ ] Add SCAD

### Polish and submit
- [ ] Polish `readme.md`
- [ ] Polish R Documents
- [ ] Polish `vignette`
- [ ] Submit to CRAN
