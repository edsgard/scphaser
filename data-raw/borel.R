


data_dir = '../nogit/data/borel'
phased_rds = file.path(data_dir, 'phased.acset.rds') ##see code/analysis/borel/phase.R

main <- function(){

    phased_res = readRDS(phased_rds)
    phased_list = phased_res[['phased_list']]
    params = phased_res[['params']]
    head(params)
    
    acset = phased_list[[1]]

    
    ##*###
    ##Dump
    ##*###
    ##store in list
    borel = acset[c('featdata', 'refcount', 'altcount', 'phenodata')]

    #dump
    devtools::use_data(borel, overwrite = T)    
    ##1.7M (too big to be part of package, would need to create a separate dataset package in that case)    
}
