

##Params
sys = 'dna'
sys = 'rs13'

##*###
##Dirs
##*###

if(sys == 'dna'){
    cloud.dir = '/mnt/kauffman/edsgard/cloud/btsync/work/rspd'
}
if(sys == 'rs13'){
    cloud.dir = '/Volumes/Data/cloud/btsync/work/rspd'
}

##OUT
rds_dir = '../nogit/data/borel'

##*###
##Files
##*###
phased_rds = file.path(rds_dir, 'phased.acset.rds')

main <- function(){

    ##import data
    phase_res = readRDS(phased_rds)
    phased_list = phase_res[['phased_list']]
    phase_params = phase_res[['params']]

    ##get phasing for one parameter-setting
    print(phase_params)
    acset = phased_list[[1]]

    ##add phased sequence to featdata df
    ##NB: Between genes we don't know the phase!
    featdata = acset[['featdata']]
    varflip = acset[['varflip']]
    vars = rownames(featdata)
    hapA = featdata[, 'ref']
    hapB = featdata[, 'alt']
    names(hapA) = vars
    names(hapB) = vars
    hapA[varflip] = featdata[varflip, 'alt']
    hapB[varflip] = featdata[varflip, 'ref']
    featdata = cbind(featdata, hapA, hapB)

    ##dump
    readr::write_tsv(featdata, file.path(rds_dir, 'featdata_phased.tsv'))
}

