

##Setup pkg
library('devtools')
devtools::create('/Volumes/Data/cloud/btsync/work/rspd/projects/scphaser/git/code/pkg/scphaser')
devtools::use_testthat()

##Workflow
devtools::load_all()
devtools::test()

