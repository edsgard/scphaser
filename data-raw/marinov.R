


data_dir = '../nogit/data/marinov'
phased_rds = file.path(data_dir, 'phased.acset.filt_vcfphased_homoase.rds')

main <- function(){

    phased_res = readRDS(phased_rds)
    phased_list = phased_res[['phased_list']]
    params = phased_res[['params']]

    acset = phased_list[[1]]

    
    ##*###
    ##Dump
    ##*###
    ##store in list
    marinov = acset[c('featdata', 'refcount', 'altcount', 'phenodata')]

    #dump
    devtools::use_data(marinov, overwrite = T)    
    
}
