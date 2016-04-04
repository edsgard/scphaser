
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/humantcell/params2sens.R")

##Dirs
out_pdf_dir = './ignore/res/humantcell/pdf'

##Libs
source('./ignore/code/analysis/performance.R')
library('reshape2')

main <- function(){


    ##*###
    ##Read and prep data
    ##*###
    load('data/humantcell.rda')

    ##create acset
    acset = new_acset(humantcell$featdata, humantcell$refcount, humantcell$altcount, humantcell$phenodata)
    lapply(acset, dim) #872, 505

    ##Call gt
    min_acount = 3
    fc = 3 #75/25
    acset = call_gt(acset, min_acount, fc)

    n.cells = c(100, 200, 300, 400, 500, 505)
    nmincells = 3
    nminvar = 2    
    acset.list = list()
    for(j.cells in n.cells){
        print(j.cells)
        
        ##subsample n.cells    
        cells = acset[['phenodata']][, 'sample']
        pass_cells = sample(cells, j.cells)
        acset_filt = subset_cols(acset, pass_cells)

        ##Filter variants on n.cells monoallelic and feats with < 2 j.vars
        acset_filt = filter_acset(acset_filt, nmincells, nminvar)
        acset.list[[as.character(j.cells)]] = acset_filt
    }
    names(acset.list) = n.cells
    
    ##get number of genes post-filtering
    ncells2nfeats = lapply(acset.list, function(j.acset){n.genes = length(unique(j.acset[['featdata']][, 'feat'])); return(n.genes)})
    ncells2nfeats = melt(ncells2nfeats)
    colnames(ncells2nfeats) = c('n.genes', 'n.cells')
    
    ##plot
    cex.lab = 2
    cex.axis = 1.5
    j.pdf = file.path(out_pdf_dir, 'ncells2ngenes.pdf')
    pdf(j.pdf)
    plot(ncells2nfeats[, 'n.cells'], ncells2nfeats[, 'n.genes'], xlab = 'n.cells', ylab = 'n.genes', pch = 16, cex.axis = cex.axis, cex.lab = cex.lab, yaxt = 'n')
    axis(2, las =2, cex.axis = cex.axis, cex.lab = cex.lab) #at = ncells2nfeats[, 'n.genes']
    dev.off()
}
