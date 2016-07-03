# scphaser
R package for haplotype phasing using single-cell RNA-seq data.

## Installation

Installation can be done via the devtools package:
```R
library('devtools')
devtools::install_github('edsgard/scphaser')
```

Upon publication scphaser will also be submitted to CRAN, the R
package repository, where installation can be done by:
```R
install.packages('scphaser')
```

## Tutorial and Reference Manual
Once you've installed scphaser you will be able to follow the vignette
which will walk you through a tutorial of how to phase variants using
scphaser. You can open it by:
```R
vignette('scphaser')
```

To get help for specific functions you can use ?fcn, for example:
```R
library('scphaser')
?phase
```

The complete function reference manual for all functions can be found
at "doc/refman.pdf" within the installed library directory. You can
also view the latest version by:
```R
browseURL('https://github.com/edsgard/scphaser/tree/master/inst/doc/refman.pdf')
```
