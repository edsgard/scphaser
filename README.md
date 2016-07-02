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
which will walk you through a tutorial of how to phase variants using scphaser, see "inst/doc/scphaser.html".

The function reference manual can be found at "inst/doc/scphaser.pdf"
or you can use the ?fcn to get help, for example:
```R
?new_acset
```
