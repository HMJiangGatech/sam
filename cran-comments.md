## Test environments
* Mac OS X, R 3.5.2
* ubuntu 16.04, R 3.5.2
* windows 10, R 3.5.2

## R CMD check results
0 errors | 0 warnings | 0 notes are found under windows 10 and Mac OS X

One NOTE under ubuntu environment:
```

* checking installed package size ... NOTE
    installed size is 22.5Mb
    sub-directories of 1Mb or more:
      libs  22.4Mb
```
It seems that on LINUX architectures, the CHECK returns one NOTE because the libs subdirectory is then above the 1MB threshold. However, it seems that this NOTE only appears under LINUX, but not under Windows or OSX. My understanding is that this inflation of the libs subdirectory is due to the use of Rcpp.


## Downstream dependencies
All dependecies have been check with 0 errors | 0 warnings | 0 notes