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

The function reference manual can be found at "doc/scphaser.pdf"
You can also use ?fcn to get help for specific functions, for example:
```R
library('scphaser')
?phase
```
