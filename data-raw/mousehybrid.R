
library('devtools')

##*###
##Read data
##*###
data_dir = '../nogit/data/mousehybrid'
files = list.files(data_dir, pattern = '.*\\.rds', full.names = TRUE)
data_list = lapply(files, readRDS)
names(data_list) = c('c57', 'cast', 'phenodata', 'featdata')
lapply(data_list, dim)

main <- function(data_list){

    ##*###
    ##featdata
    ##*###

    ##set row- and col-names
    featdata = data_list[['featdata']]
    var = featdata[['chrom:pos']]
    length(unique(var))
    rownames(featdata) = var
    colnames(featdata) = c('feat', 'chr', 'pos', 'c57', 'cast', 'var')

    ##factor -> chars
    factor.cols = unlist(lapply(featdata, is.factor))
    featdata[, factor.cols] = lapply(featdata[, factor.cols], as.character)


    ##*###
    ##Count data
    ##*###
    refcount = data_list[['c57']]
    altcount = data_list[['cast']]
    rownames(refcount) = featdata[, 'var']
    rownames(altcount) = featdata[, 'var']

    ##filter feats on number of vars
    ##get n variants per feature
    feat2nvar = table(featdata[, 'feat'])

    ##get passed feats
    nmin_var = 2
    pass_feat = names(feat2nvar)[which(feat2nvar >= nmin_var)]
    length(pass_feat)
    length(unique(featdata[, 'feat']))

    acset = list(featdata = featdata, refcount = refcount, altcount = altcount)
    acset = filter_feat(acset, pass_feat)
        
    ##filter feats on sum of counts across all vars and cells (minimum example: 2 vars, 2 cells)
    mincount = 3
    acset = filter_featcount(acset, mincount)
    
    
    ##*###
    ##phenodata
    ##*###

    pheno = data_list[['phenodata']]

    ##factor -> chars
    factor.cols = unlist(lapply(pheno, is.factor))
    pheno[, factor.cols] = lapply(pheno[, factor.cols], as.character)

    ##set row- and col-names
    colnames(pheno)[1] = 'sample'
    rownames(pheno) = pheno[, 'sample']

    ##ensure order is the same in pheno as in count matrixes
    nrow(pheno) == ncol(refcount)
    nrow(pheno) == length(which(pheno[, 'sample'] == colnames(refcount)))
    nrow(pheno) == length(which(colnames(altcount) == colnames(refcount)))
    pheno = pheno[colnames(refcount), ]


    ##*###
    ##Dump
    ##*###
    ##store in list
    mousehybrid = c(acset, list(phenodata = pheno))

    ##dump
    save(mousehybrid, file = '../nogit/data/mousehybrid/mousehybrid.rda')
    
    ##TODO: dump only one chr in pkg (since maximum 1M)
    devtools::use_data(mousehybrid_chr1)
    
}

filter_featcount <- function(acset, mincount){

    featdata = acset[['featdata']]
    refcount = acset[['refcount']]
    altcount = acset[['altcount']]
    
    feat2var = tapply(featdata[, 'var'], featdata[, 'feat'], unique)

    ##set na to 0
    refcount[is.na(refcount)] = 0
    altcount[is.na(altcount)] = 0
    totcount = refcount + altcount

    ##counts per feat
    var2count = apply(totcount, 1, sum)
    length(which(var2count == 0))
    feat2count = unlist(lapply(feat2var, function(vars, counts){sum(counts[vars, ])}, counts = totcount))
    ##TBD: slow, consider using dplyr
    pass_feat = names(feat2count)[which(feat2count > mincount)]

    ##filter on feats
    acset = filter_feat(acset, pass_feat)

    return(acset)
}

filter_feat <- function(acset, pass_feat){

    featdata = acset[['featdata']]
    refcount = acset[['refcount']]
    
    ##subset featdata on passed feats to get passed variants
    featdata = merge(featdata, as.matrix(pass_feat), by.x = 'feat', by.y = 1, stringsAsFactors = FALSE)
    rownames(featdata) = featdata[, 'var']

    ##subset on passed vars
    pass_var = featdata[, 'var']
    refcount = refcount[pass_var, ]
    altcount = altcount[pass_var, ]

    acset[['featdata']] = featdata
    acset[['refcount']] = refcount
    acset[['altcount']] = altcount

    return(acset)
}
