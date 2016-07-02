
###Get stats: Number of expr genes and seq-depth
###NB: The paths are relative to the root of the scphaser git repo
###
###SYNOPSIS
###library('devtools')
###devtools::load_all()
###source("ignore/code/analysis/borel/stats.R")


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

data.dir = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel')
mouse.data.dir = file.path(cloud.dir, 'projects/scphaser/nogit/data/mousehybrid')

##OUT
out_rds_dir = '../nogit/data/borel'
out_pdf_dir = './ignore/res/borel/pdf'
combo_pdf_dir = './ignore/res/combo/pdf'


##*###
##Files
##*###

##IN
mapstats.tab = file.path(data.dir, 'mapstats.tab')
ref.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', paste('ref', '.counts.rds', sep = ''))
alt.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', paste('alt', '.counts.rds', sep = ''))
snp.annot.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', 'snp.annot.rds')
rpkm.tab = file.path(data.dir, 'rpkm_refseq.txt')

mouse.rpkm.tab = file.path(mouse.data.dir, 'rpkm/rpkms_counts_rmnameoverlap.txt')
sample2id.tab = file.path(mouse.data.dir, 'old.new.cell_names.txt')

##OUT

##Libs
source('./ignore/code/analysis/performance.R')
library('BiocParallel')
library('dplyr')
library('tidyr')

main <- function(){

    ##*###
    ##Read and prep data
    ##*###
    altcount = readRDS(alt.counts.rds)
    refcount = readRDS(ref.counts.rds)
    
    ##totcount (dep on het.vars)
    tot = altcount + refcount
    cell2countsum = apply(tot, 2, sum)
    summary(cell2countsum) #mean: 505,300, median: 526,500, min: 40,870, max: 1,581,000


    
    ##*############################
    ##Sequencing depth
    ##*############################
    
    ##number of uniquely mapping reads (intergenic and genic)
    mapstats = read.table(mapstats.tab, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
    colnames(mapstats) = sub('\\.\\.\\.\\.', '', colnames(mapstats))

    n.uniq = mapstats[, 'Reads'] * mapstats[, 'Uniqely.mapping'] / 100
    mapstats = cbind(mapstats, n.uniq)

    ##Filter on the 163 cells
    rownames(mapstats) = mapstats[, 'Sample']
    mapstats = mapstats[colnames(altcount), ]

    summary(mapstats[, 'n.uniq'] / 1e6)
    ##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 0.9097  6.4010 12.1200 11.9900 15.8400 31.6100

    hist(mapstats[, 'n.uniq'])



    ##*########################
    ##BOREL
    ##*########################

    ##Read data
    rpkm.counts.list = read.rpkm(rpkm.tab)
    rpkm = rpkm.counts.list[['rpkm']]
    count = rpkm.counts.list[['count']]
    dim(rpkm)
    dim(count)

    ##filter on 163 single cells analyzed by scphaser
    pass.samples = colnames(altcount)
    rpkm = rpkm[, pass.samples]
    count = count[, pass.samples]

    
    ##*###
    ##Seqdepth at RefSeq genes
    ##*###

    ##reads aligned to RefSeq genes
    cell2nreads = apply(count, 2, sum)
    median(cell2nreads) #7,924,512
    
    hist(log10(cell2nreads))
    ##about 10M reads
    
    ##*###
    ##Number of expressed genes
    ##*###
        
    ##mean expr per gene
    gene2mean = apply(rpkm, 1, mean)
    hist(log10(gene2mean))

    ##filter genes on mean expr (quite strict filter)
    min.mean = 5
    pass.genes = names(gene2mean)[which(gene2mean >= min.mean)]
    length(pass.genes) #7,141

    
    ##*############################
    ##MOUSE
    ##*############################
    
    ##Read data
    mouse.rpkm.count.list = read.table(mouse.rpkm.tab, sep = '\t', stringsAsFactors = FALSE, header = FALSE)
    rpkm = rpkm.count.list[['rpkm']]
    count = rpkm.count.list[['count']]

    ##filter on 336 cells analyzed by scphaser
    sample2id = read.table(sample2id.tab, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
    pass.samples = sample2id[, 'old.name']
    setdiff(pass.samples, colnames(rpkm)) #BQx34_indG_EmbryoMEF_BxC
    colnames(rpkm)[grep('BQx34', colnames(rpkm))] #"BQx34_indG_EmbryoMEF_BxCM "
    colnames(rpkm) = sub('BQx34_indG_EmbryoMEF_BxCM ', 'BQx34_indG_EmbryoMEF_BxC', colnames(rpkm))
    colnames(count) = sub('BQx34_indG_EmbryoMEF_BxCM ', 'BQx34_indG_EmbryoMEF_BxC', colnames(count))
    rpkm = rpkm[, pass.samples]
    count = count[, pass.samples]

    ##reads aligned to RefSeq genes
    cell2nreads = apply(count, 2, sum)    
    hist(log10(cell2nreads))
    median(cell2nreads) #4,996,822
    
    ##mean expr per gene
    gene2mean = apply(rpkm, 1, mean)
    hist(log10(gene2mean))

    ##filter genes on mean expr (quite strict filter)
    min.mean = 5
    pass.genes = names(gene2mean)[which(gene2mean >= min.mean)]
    length(pass.genes) #8,545

    ##store
    mmu.cell2nreads = cell2nreads
    mmu.gene2mean = gene2mean


    ##*##############
    ##COMBO
    ##*##############

    ##*###
    ##cell2nreads
    ##*###

    ##vector -> df
    cell2nreads.mat = as.data.frame(cbind(cell2nreads, rep('Borel', length(cell2nreads))), stringsAsFactors = FALSE)
    colnames(cell2nreads.mat) = c('nreads', 'dataset')
    cell2nreads.mat[, 'nreads'] = as.integer(cell2nreads.mat[, 'nreads'])

    mmu.cell2nreads.mat = as.data.frame(cbind(mmu.cell2nreads, rep('Mouse-hybrid', length(mmu.cell2nreads))), stringsAsFactors = FALSE)
    colnames(mmu.cell2nreads.mat) = c('nreads', 'dataset')
    mmu.cell2nreads.mat[, 'nreads'] = as.integer(mmu.cell2nreads.mat[, 'nreads'])
    
    cell2nreads.df = rbind(cell2nreads.mat, mmu.cell2nreads.mat)
    dim(cell2nreads.df) #499

    ##construct gg
    gg = ggplot(cell2nreads.df, aes_string(x = 'nreads', color = 'dataset'))
    gg = gg + geom_density()
    gg = gg + scale_x_log10() #breaks = c(0, 10, 100, 1000, 10000, 100000), labels = c(0, 10, 100, 1000, 10000, 100000)

    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black", size = 10), axis.ticks = element_line(colour = 'black'))
    gg = gg + xlab('Number of reads aligned to RefSeq genes')
    gg = gg + ylab('Number of cells (density)')
    gg = gg + theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5))

    ##plot
    pdf.h = 4
    pdf.w = 5
    j.pdf = file.path(combo_pdf_dir, 'cell2nreads_refseq.dens.pdf')
    pdf(j.pdf, width = pdf.w, height = pdf.h)
    plot(gg)
    dev.off()


    ##*###
    ##gene2meanexpr
    ##*###

    ##vector -> df
    gene2mean.mat = as.data.frame(cbind(gene2mean, rep('Borel', length(gene2mean))), stringsAsFactors = FALSE)
    colnames(gene2mean.mat) = c('mean.expr', 'dataset')
    gene2mean.mat[, 'mean.expr'] = as.integer(gene2mean.mat[, 'mean.expr'])

    mmu.gene2mean.mat = as.data.frame(cbind(mmu.gene2mean, rep('Mouse-hybrid', length(mmu.gene2mean))), stringsAsFactors = FALSE)
    colnames(mmu.gene2mean.mat) = c('mean.expr', 'dataset')
    mmu.gene2mean.mat[, 'mean.expr'] = as.integer(mmu.gene2mean.mat[, 'mean.expr'])
    
    gene2mean.df = rbind(gene2mean.mat, mmu.gene2mean.mat)
    dim(gene2mean.df) #49137

    ##filter out genes with 0 mean (only want to compare the number of expressed genes)
    min.rpkm = 5
    gene2mean.df = gene2mean.df[which(gene2mean.df[, 'mean.expr'] > min.rpkm), ]
    
    ##add pseudocount
    ##pseudo.count = 0.01
    ##gene2mean.df[, 'mean.expr'] = gene2mean.df[, 'mean.expr'] + pseudo.count
    
    ##construct gg
    gg = ggplot(gene2mean.df, aes_string(x = 'mean.expr', color = 'dataset'))
    gg = gg + geom_freqpoly(bins = 17)
    gg = gg + scale_x_log10(breaks = c(0, 10, 100, 1000, 10000), labels = c(0, 10, 100, 1000, 10000))

    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black", size = 10), axis.ticks = element_line(colour = 'black'))
    gg = gg + xlab('Mean expression (RPKM)')
    gg = gg + ylab('Number of RefSeq genes')
    gg = gg + theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5))

    ##plot
    pdf.h = 4
    pdf.w = 6
    j.pdf = file.path(combo_pdf_dir, 'gene2mean_refseq.freqpoly.pdf')
    
    pdf(j.pdf, width = pdf.w, height = pdf.h)
    plot(gg)
    dev.off()

    
}

read.rpkm <- function(rpkm.tab, samplesfromcall = FALSE){

    ##read
    rpkm.count = read.table(rpkm.tab, sep = '\t', stringsAsFactors = FALSE, header = FALSE)
    header = read.table(rpkm.tab, sep = '\t', stringsAsFactors = FALSE, nrows = 4, comment.char = '', fill = TRUE)
    rownames(header) = sub('^#', '', header[, 1])
    header = header[, -1]
    samples = header['samples', ]

    if(samplesfromcall){
        arg.in = sub('.*-i', '', header['arguments', 1])
        arg.in = sub('-o.*', '', arg.in)
        arg.in = sub('^ ', '', arg.in)
        arg.in = sub(' $', '', arg.in)
        arg.in = strsplit(arg.in, ' ')[[1]]
        arg.in = basename(dirname(arg.in))
        samples = arg.in
    }
    
    ##get feats
    feat.annot = rpkm.count[, 1:2]
    colnames(feat.annot) = c('hugo', 'nm')

    ##Split rpkm.count -> rpkm, count
    rpkm.count = rpkm.count[, -(1:2)]
    n.cols = ncol(rpkm.count)
    rpkm = rpkm.count[, 1:(n.cols / 2)]
    count = rpkm.count[, (n.cols/2 + 1):n.cols]

    colnames(rpkm) = samples
    colnames(count) = samples
    
    ##filter ERCC
    ercc.ind = grep('^ERCC-', feat.annot[, 'nm'])
    length(ercc.ind) #92
    pass.ind = setdiff(1:nrow(feat.annot), ercc.ind)
    rpkm = rpkm[pass.ind, ]
    count = count[pass.ind, ]
    feat.annot = feat.annot[pass.ind, ]

    rpkm.count.list = list(feat.annot = feat.annot, rpkm = rpkm, count = count)
    return(rpkm.count.list)
}
