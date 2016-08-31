


main <- function(){

    ##*###
    ##Installations
    ##*###
    install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
    devtools::install_github("hadley/devtools")
    install.packages("rmarkdown")
    ##pandoc: http://pandoc.org/installing.html
    ##after basictex installation:
    ##sudo tlmgr update --self
    ##sudo tlmgr install collection-fontsrecommended
    ##sudo tlmgr install inconsolata
    
    ##check devtools installation (restart R if you just installed it)
    library('devtools')
    has_devel()

    
    ##*###
    ##Setup pkg
    ##*###
    library('devtools')
    devtools::create('/Volumes/Data/cloud/btsync/work/rspd/projects/scphaser/git/code/pkg/scphaser')
    devtools::use_testthat()
    devtools::use_vignette("scphaser")
    devtools::use_data_raw() ##creates data-raw dir containing code for parsing and storing data and adds that dir to .Rbuildignore
    
    ##Ignore dirs
    devtools::use_build_ignore("ignore")    

    
    ##*###
    ##Test
    ##*###
    devtools::load_all()
    devtools::test()
    devtools::check()


    ##*###
    ##Documentation and namespace
    ##*###
    devtools::document()

    
    ##*###
    ##Reference Manual
    ##*###
    ##R CMD Rd2pdf git
    ##cp -p git.pdf git/inst/doc/refman.pdf
    ##I cant use scphaser.pdf since then scphaser.pdf overrides scphaser.html in the generated doc/index.html upon installation by install_github, which means that browseVignettes links to the reference manual instead.

    
    ##*##############
    ##Vignettes
    ##*##############

    ##*##
    ##Test vignette
    ##*##
    ##spin
    knitr::spin('ignore/code/vignettes/test.R', knit = FALSE)    
    
    ##add meta header
    file.copy('ignore/code/vignettes/meta.Rmd', 'vignettes/meta.Rmd')
    file.append('vignettes/meta.Rmd', 'ignore/code/vignettes/test.Rmd')
    file.rename('vignettes/meta.Rmd', 'vignettes/test.Rmd')

    ##*##
    ##scphaser vignette
    ##*##
    ##spin
    knitr::spin('ignore/code/vignettes/scphaser.R', knit = FALSE, format = 'Rmd', precious = TRUE)    
    
    ##add meta header
    file.copy('ignore/code/vignettes/meta.Rmd', 'vignettes/meta.Rmd')
    file.append('vignettes/meta.Rmd', 'ignore/code/vignettes/scphaser.Rmd')
    file.rename('vignettes/meta.Rmd', 'vignettes/scphaser.Rmd')

    ##knit
    devtools::build_vignettes() ##copied to inst/doc

    
    ##*###
    ##Pre-flight Readme.md -> html
    ##*###
    ##grip git/README.md --export README.html


    ##*###
    ##CRAN release
    ##*###
    ##cran-comments.md
    devtools::build()
    ##R CMD check --as-cran scphaser_1.0.0.tar.gz
    devtools::build_win(version = 'R-release')
    devtools::build_win(version = 'R-devel')
    devtools::revdep_check()

    
    ##*###
    ##Check citation
    ##*###
    ##doesnt work: install.packages(pkg = 'file:///mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/scphaser_1.0.0.tar.gz', repos = NULL, type = 'source')
    ##R CMD INSTALL scphaser_1.0.0.ubuntu.tar.gz
    citation('scphaser')
    vignette('scphaser')
    remove.packages('scphaser')

    
    ##*##
    ##Before submission
    ##*###
    ##rerun devtools::document (man/*.Rd newer than R/data.R)
    ##rerun devtools::build_vignettes() above
    ##rerun R CMD Rd2pdf above
    ##check if citation has been updated
    ##check vignette
    
    
}
