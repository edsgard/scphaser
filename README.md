## About
scphaser is an R package for haplotype phasing using single-cell RNA-seq data.

scphaser is available free to use under the <a href="./inst/doc/LICENSE">GNU GPL version 3
license</a>.

## Installation
scphaser is available at the R package repository CRAN and can be installed by:
```R
source("http://www.bioconductor.org/biocLite.R")
biocLite("BiocParallel")
install.packages('scphaser')
```

Installation of the development version can be done via the devtools package:
```R
library('devtools')
devtools::install_github('edsgard/scphaser')
```

## Tutorial
Once you've installed scphaser you will be able to follow the
vignette-tutorial. You can open it by:
```R
vignette('scphaser')
```

## Minimal example
```R
library('scphaser')
data(marinov)
acset = new_acset(featdata = marinov[['featdata']], refcount = marinov[['refcount']],
	altcount = marinov[['altcount']], phenodata = marinov[['phenodata']])
acset = call_gt(acset)
acset = filter_acset(acset)
acset = phase(acset)
head(acset[['phasedfeat']])
```

## Function reference manual
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

## Citation

If you use scphaser, please cite it as follows:

Edsgärd D. <em>et al.</em>, scphaser: Haplotype Inference Using
Single-Cell RNA-Seq Data, <em>Bioinformatics</em>, 2016<br>
<a href="http://dx.doi.org/10.1093/bioinformatics/btw484">DOI: 10.1093/bioinformatics/btw484</a>
