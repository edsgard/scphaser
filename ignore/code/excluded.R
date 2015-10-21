
filter_nminmono <- function(acset, nmin_mono = 1){
###filter variants having less than nmin_mono mono-allelic calls of each allele across samples
    
    gt = acset[['gt']]
    pass_vars = which(rowSums(gt == 0, na.rm = TRUE) >= nmin_mono & rowSums(gt == 2, na.rm = TRUE) >= nmin_mono)
    acset = subset_rows(acset, pass_vars)

    ##store filter argument
    acset[['args']][['filter']]['nmin_mono'] = list(nmin_mono = nmin_mono)

    return(acset)
}

get_cv_ase <- function(ase){

    ##sum the cvs from every cell (assume cells are independent observations)
    cvs = apply(ase, 2, function(jase){sd(jase) / mean(jase)})
    cvtot = sum(cvs, na.rm = TRUE)

    return(cvtot)
}

get_wcv_ase <- function(ase, weights){

    ##sum the weighted cv2 from every cell (assume cells are independent observations)
    ncells = ncol(ase)
    cv2 = rep(NA, ncells)
    for(jcell in 1:ncells){
        jase = ase[, jcell]
        jw = weights[, jcell]
        cv2[jcell] = Hmisc::wtd.var(jase, jw) / (Hmisc::wtd.mean(jase, jw)^2)
    }
    cvtot = sum(cv2)

    return(cvtot)
}
