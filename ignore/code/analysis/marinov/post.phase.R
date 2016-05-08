
###Assess complexity under varying parameters
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/marinov/post.phase.R")


##Params
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

##IN
rds_dir = './ignore/res/marinov/data/filt_vcfphased_randomase'

##OUT
pdf_dir = './ignore/res/marinov/pdf/filt_vcfphased_randomase'

##*###
##Files
##*###

##IN
phased_rds = file.path(rds_dir, 'phased.acset.rds')
vcf.phased.rsid = file.path(cloud.dir, 'data/external/Marinov_GenRes_2014', 'hg19.trio.het_phased.pass.rsid')

##Libs
source('./ignore/code/analysis/performance.R')
library('dplyr')
library('tidyr')

main <- function(){

    ##Read data
    phased_res = readRDS(phased_rds)
    phased_list = phased_res[['phased_list']]
    params = phased_res[['params']]
    
    
    ##*###
    ##Filter SNPs on phased VCF entry ("|")
    ##*###
    vcf.rsid = read.table(vcf.phased.rsid, stringsAsFactors = FALSE)[, 1]
    phased.filt = lapply(phased_list, function(j.acset, vcf.rs){rs = j.acset[['featdata']][, 'rsid']; pass.rs = intersect(rs, vcf.rs); j.acset = subset_rows(j.acset, pass.rs); return(j.acset);}, vcf.rsid)


    ##*###
    ##Filter SNPs on number of cells where ASE towards the same allele
    ##*###

    alpha = 0.1
    mono.ase = 0.2
    if(!mono.ase == 0){
        phased.filt2 = lapply(phased.filt, filter_homovars, alpha = alpha, mono.ase = mono.ase)
    }else{
        phased.filt2 = phased.filt
    }
        
    nvars.pre = unlist(lapply(phased_list, function(j.acset){nrow(j.acset[['altcount']])}))
    nvars.pre = unlist(lapply(phased.filt, function(j.acset){nrow(j.acset[['altcount']])}))
    nvars.post = unlist(lapply(phased.filt2, function(j.acset){nrow(j.acset[['altcount']])}))
    frac.left = nvars.post / nvars.pre
    summary(frac.left)
    ##mono.ase = 0.5: mean: 79, median: 83
    ##mono.ase = 0.2: 81, 83
    ##0.1 -> 82, 84
    ##0.05 -> 82, 84

    
    ##*###
    ##Get concordance
    ##*###
    
    perf.list = lapply(phased.filt2, check_phase)
    feat.perf = do.call(rbind, lapply(perf.list, '[[', 'feat'))
    var.perf = do.call(rbind, lapply(perf.list, '[[', 'var'))
    perf = list(feat = feat.perf, var = var.perf)
    param2perf = as.data.frame(perf)
    param2perf = cbind(params, perf)
    param2perf = tidyr::unite(param2perf, in.meth.w, input, method, weigh, remove = FALSE, sep = '.')
    weigh.num = as.numeric(param2perf[, 'weigh'])
    param2perf = cbind(param2perf, weigh.num);
    which(unlist(lapply(param2perf, is.factor)))
    ##TBD: check treatment of NAs in check_phase    

    ##Dump
    j.rds = file.path(rds_dir, paste('monoase_', mono.ase, sep = ''), 'param2perf.rds')
    dir.create(dirname(j.rds), recursive = TRUE)
    saveRDS(param2perf, file = j.rds)

    
    ##*###
    ##Plot perf
    ##*###
    library('ggplot2')

    j.pdf.dir = file.path(pdf_dir, paste('monoase_', mono.ase, sep = ''))
    dir.create(j.pdf.dir)
    feat.level = names(perf)
    y.resp = c('Number of', 'Fraction erronously phased')
    names(y.resp) = c('n', 'frac.err')
    
    for(j.y.resp.name in names(y.resp)){
        for(j.feat.level in feat.level){

            x.str = 'min_acount'
            group.col = 'in.meth.w'
            colour.col = 'method'
            lt.col = 'input'
            size.col = 'weigh.num'
            x.lab = 'Minimum allele count of at least one allele'
            y.lab = paste(y.resp[j.y.resp.name], j.feat.level, sep = ' ')
            y.str = paste(j.feat.level, j.y.resp.name, sep = '.')
            
            ##data and mappings
            gg = ggplot(param2perf, aes_string(x = x.str, y = y.str, colour = colour.col, group = group.col, linetype = lt.col, size = size.col))
            gg = gg + scale_size(range = c(0.5, 1), breaks = c(0, 1), labels = c('FALSE', 'TRUE'))

            ##glyph
            gg = gg + geom_line()

            ##facet
            gg = gg + facet_grid(fc ~ nmincells, scales = 'free_y')
            
            ##Background
            gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())

            ##ticks
            gg = gg + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))

            ##gg = gg + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
            
            ##axis lables
            gg = gg + xlab(x.lab) + ylab(y.lab)

            j.pdf = file.path(j.pdf.dir, paste(j.feat.level, j.y.resp.name, 'pdf', sep = '.'))
            pdf(file = j.pdf)
            print(gg)
            dev.off()
        }
    }    

    ##table
    j.rds = file.path(rds_dir, paste('monoase_', mono.ase, sep = ''), 'param2perf.rds')
    param2perf = readRDS(j.rds)
    param2perf.filt = dplyr::filter(param2perf, nmincells == 5, fc == 3, min_acount == 3)
    ##gt.exhaust.true:
    ##0.05 -> 0.045, n.err=24, n.feat=533, n.vars=1,367
    ##0.1 -> 0.043, n.err=23, n.feat=534, n.vars=1,363
    ##0.2 -> 0.045, n.err=24, n.feat=530, n.vars=1352
    ##ac.ex.false:
    ##0.05 -> 0.053, 28
    ##0.1 -> 0.0506, n.err=27
    ##0.2 -> 0.053, n.err=28

    ##*###
    ##mouse-hybrid data
    ##*###

    ##read data
    mouse.rds = './ignore/res/mousehybrid/perf/data/perf.ngenes_100.nparams_8.rds'
    m.param2perf = readRDS(mouse.rds)

    ##add weigh.num
    weigh.num = as.numeric(m.param2perf[, 'weigh'])
    m.param2perf = cbind(m.param2perf, weigh.num);
    which(unlist(lapply(m.param2perf, is.factor)))

    ##median
    which.quantile <- function (x, probs, na.rm = FALSE){
        if (! na.rm & any (is.na (x)))
            return (rep (NA_integer_, length (probs)))

        o <- order (x)
        n <- sum (! is.na (x))
        o <- o [seq_len (n)]

        nppm <- n * probs - 0.5
        j <- floor(nppm)
        h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
        j <- j + h

        j [j == 0] <- 1
        o[j]
    }
    
    a = split(m.param2perf, m.param2perf[, 'in.meth.w'])
    b = lapply(a, function(j.perf){median.ind = which.quantile(j.perf[, 'feat.frac.err'], probs = 0.5, na.rm = TRUE); return(j.perf[median.ind, ]);})
    b = lapply(a, function(j.perf){feat.frac.err = median(j.perf[, 'feat.frac.err']); var.frac.err = median(j.perf[, 'var.frac.err']); res = c(j.perf[1, 'in.meth.w'], feat.frac.err, var.frac.err);})
    m.summ = as.data.frame(do.call(rbind, b), stringsAsFactors = FALSE)
    colnames(m.summ) = c('in.meth.w', 'feat.frac.err', 'var.frac.err')
    num.cols = c('feat.frac.err', 'var.frac.err')
    m.summ[, num.cols] = lapply(m.summ[, num.cols], as.numeric)
    
    
    ##fix colnames to be the same
    m.cols = colnames(m.param2perf)
    m.cols[grep('nvars_postfilt', m.cols)] = 'var.n'
    m.cols[grep('nfeats_postfilt', m.cols)] = 'feat.n'    
    colnames(m.param2perf) = m.cols    
    
    ##bind human and mouse perf
    h.cols = colnames(param2perf)
    m.cols = colnames(m.summ)
    shared.cols = intersect(h.cols, m.cols)

    dataset = rep('mouse', nrow(m.summ))
    m.summ = cbind(m.summ, dataset, stringsAsFactors = FALSE)
    dataset = rep('human', nrow(m.summ))
    param2perf.filt = cbind(param2perf.filt, dataset, stringsAsFactors = FALSE)
    shared.cols = c(shared.cols, 'dataset')
    comb.perf = rbind(m.summ[, shared.cols], param2perf.filt[, shared.cols])

    
    ##*###
    ##Barplot
    ##*###
    gg.data = comb.perf
    var.frac.right = 1 - gg.data[, 'var.frac.err']
    gg.data = cbind(gg.data, var.frac.right)
    gg = ggplot(gg.data, aes_string(x = 'in.meth.w', y = 'var.frac.right'))
    gg = gg + geom_bar(stat = 'identity')
    gg = gg + facet_grid(~ dataset)
    gg = gg + xlab('') + ylab('Fraction correctly phased genes')
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))
    gg = gg + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    gg = gg + geom_hline(yintercept = 0.95, colour = 'gray', linetype = 'dashed')
    ##gg = gg + scale_y_continuous(breaks = c(seq(0, 1, 0.1), 0.95))
    gg = gg + scale_y_continuous(breaks = c(seq(0, 1, 0.25), 0.9, 0.95))
    j.pdf = file.path(pdf_dir, 'frac.vars.err.mouse_human.pdf')
    pdf.h = 3
    pdf.w = 5
    pdf(j.pdf, height = pdf.h, width = pdf.w)
    print(gg)
    dev.off()
    
    ##*###
    ##TODO: Add inferred haplotype column to acset (ADD TO R-PACKAGE)
    ##*###
    acset = phased_list[[1]]
    featdata = acset[['featdata']]
    varflip = acset[['varflip']]
    vars = rownames(featdata)
    hapA = featdata[, 'ref']
    hapB = featdata[, 'alt']
    names(hapA) = vars
    names(hapB) = vars
    hapA[varflip] = featdata[varflip, 'alt']
    hapB[varflip] = featdata[varflip, 'ref']
    ##NB: Between genes we don't know the phase!

}

check_phase <- function(acset){

    featdata = acset[['featdata']]
    varflip = acset[['varflip']]    

    ##binarize vars to be flipped
    var = featdata[, 'var']
    feat = featdata[, 'feat']
    obs = as.integer(var %in% varflip)
    exp = as.integer(featdata[, 'ref'] == featdata[, 'mat.allele'])
    exp_compl = exp
    exp_compl[exp == 0] = 1
    exp_compl[exp == 1] = 0
    varflip.df = as.data.frame(cbind(var, feat, obs, exp, exp_compl), stringsAsFactors = FALSE)
    int.cols = c('obs', 'exp', 'exp_compl')
    varflip.df[, int.cols] = lapply(varflip.df[, int.cols], as.integer)
    rownames(varflip.df) = varflip.df[, 'var']
                
    ##set the observed flip-vector to the complementary of obs, if obs is closer to exp_compl
    feat2vars = tapply(varflip.df[, 'var'], varflip.df[, 'feat'], unique)
    varflip.list = lapply(feat2vars, function(var, varflip.df){
        jflip = varflip.df[var, ]
        exp_error = length(which(jflip[, 'obs'] != jflip[, 'exp']))
        exp_compl_error = length(which(jflip[, 'obs'] != jflip[, 'exp_compl']))
        obs_adj = jflip[, 'obs']
        if(exp_compl_error < exp_error){
            obs = jflip[, 'obs']
            obs_adj[which(obs == 0)] = 1
            obs_adj[which(obs == 1)] = 0
        }
        jflip = cbind(jflip, obs_adj)

        return(jflip)
    }, varflip.df)
    varflip.df = do.call(rbind, varflip.list)
    rownames(varflip.df) = varflip.df[, 'var']
    varflip.df[, 'obs_adj'] = as.integer(varflip.df[, 'obs_adj'])
    
    ##var-level errors
    n.var.err = length(which(varflip.df[, 'exp'] != varflip.df[, 'obs_adj']))
    n.vars = nrow(varflip.df)
    frac.err = n.var.err / n.vars
    var.perf = c(n.vars, n.var.err, frac.err)
    names(var.perf) = c('n', 'n.err', 'frac.err')
    
    ##feat-level errors
    feat2vars = tapply(varflip.df[, 'var'], varflip.df[, 'feat'], unique)    
    feat2correct = unlist(lapply(feat2vars, function(vars, varflip.df){
        jflip = varflip.df[vars, ]
        identical(jflip[, 'exp'], jflip[, 'obs_adj'])
    }, varflip.df))
    n.feat.err = length(which(!feat2correct))
    n.feats = length(feat2vars)
    frac.err = n.feat.err / n.feats
    feat.perf = c(n.feats, n.feat.err, frac.err)
    names(feat.perf) = c('n', 'n.err', 'frac.err')

    ##store
    perf = list(varflip.df = varflip.df, var = var.perf, feat = feat.perf)

    return(perf)
}
