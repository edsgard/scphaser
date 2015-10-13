
get_initial_gt <- function(acset){
###We want to start with the matrix which induce an increase in the compactness score
###TBD: Don't think I'm actually doing something sensible here! Because, if the number of ref monoallelic calls are greater why can't we just as well start with that?!

    ##not finished...
    
    ##Calculate the fraction of alternative monoallelic calls, for each SNP
    n_samples = ncol(jfeat_gt)
    var2fracalt = rowSums(jfeat_gt == 2) / n_samples

    ##Calculate the fraction of SNPs with var2fracalt >0.5
    n_vars = nrow(jfeat_gt)
    frac_vars = length(which(var2fracalt > 0.5)) / n_vars
    
    ##If frac_vars > 0.5 start out with SNPs with var2fracalt < 0.5, otherwise the other set of SNPs.
    if(frac_vars > 0.5){
        jfeat_gt_phased = jfeat_gt
    }else{
        jfeat_gt_phased = jfeat_gt_comple
    }

    return(jfeat_gt_phased)
}
