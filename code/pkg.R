

main <- function(){

    ##Installations
    install.packages(c("devtools", "roxygen2", "testthat", "knitr"))
    devtools::install_github("hadley/devtools")
    install.packages("rmarkdown")
    ##pandoc: http://pandoc.org/installing.html

    ##check devtools installation (restart R if you just installed it)
    library(devtools)
    has_devel()
    
    ##Setup pkg
    library('devtools')
    devtools::create('/Volumes/Data/cloud/btsync/work/rspd/projects/scphaser/git/code/pkg/scphaser')
    devtools::use_testthat()
    devtools::use_vignette("scphaser")
    devtools::use_data_raw() ##creates data-raw dir containing code for parsing and storing data and adds that dir to .Rbuildignore
    
    ##Ignore dir
    devtools::use_build_ignore("doc")    
    devtools::use_build_ignore("code")
    devtools::use_build_ignore("res")
    
    ##Workflow
    devtools::load_all()
    devtools::test()
    devtools::check()
    
}
