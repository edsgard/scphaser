


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

    ##Reference Manual
    ##R CMD Rd2pdf git
    ##cp -p git.pdf git/inst/doc/scphaser.pdf
    
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
    devtools::build_vignettes()
    devtools::build() ##this will also build vignettes

    
}
