## Test environments
* Mac OS X, R 4.1.0
* ubuntu 18.04, R 4.1.0
* windows 10, R 4.1.0

## News
Compared to 1.1.2, we add weighted loss for SVM algorithm.

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE on Mac OS:
* checking installed package size ... NOTE
  installed size is 21.0Mb
  sub-directories of 1Mb or more:
    libs  20.9Mb
    
  I think this is due to the Rcpp/RcppEigen dependency.

There was 1 NOTE on windows:
* Note: information on .o files for i386 is not available
It seems not a big problem https://stackoverflow.com/questions/64402688/information-on-o-files-for-x64-is-not-available-note-on-r-package-checks-using .

## revdepcheck results

We checked 3 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

