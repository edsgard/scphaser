
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

##OUT
out_rds_dir = '../nogit/data/borel'
out_pdf_dir = './ignore/res/borel/pdf'


##*###
##Files
##*###

##IN
mapstats.tab = file.path(data.dir, 'mapstats.tab')
ref.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', paste('ref', '.counts.rds', sep = ''))
alt.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', paste('alt', '.counts.rds', sep = ''))
snp.annot.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', 'snp.annot.rds')
rpkm.tab = file.path(data.dir, 'rpkm_refseq.txt')

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
    featdata = readRDS(snp.annot.rds)
    
    
    ##*#############################
    ##Distance between neighboring variants (within genes) at the DNA level
    ##*#############################

    ##create acset
    acset = new_acset(featdata, refcount, altcount)
    nrow(featdata) #231,075
    
    ##Filter vars with 0 counts
    acset = filter_zerorow(acset)
    lapply(acset, dim) #230,674    163

    featdata = acset[['featdata']]
    length(unique(acset[['featdata']][, 'feat'])) #13,688
    
    ##filter on feats with at least two vars
    nminvar = 2
    acset = filter_feat_nminvar(acset, nminvar)    
    lapply(acset, dim) #227,663 x 163
    length(unique(acset[['featdata']][, 'feat'])) #10,677

    featdata = acset[['featdata']]

    
    ##split by chr.gene
    chr.feat = apply(featdata[, c('chr', 'feat')], 1, paste, collapse = '.')
    featdata = cbind(featdata, chr.feat, stringsAsFactors = FALSE)
    
    ##for each gene get distance
    feat.list = split(featdata, featdata[, 'chr.feat'])
    feat2diff = unlist(lapply(feat.list, function(j.feat){diff(sort(j.feat[, 'pos']))}))

    ##exclude outliers (two)
    outliers = names(feat2diff)[which(feat2diff > 1e7)]
    feat2diff = feat2diff[setdiff(names(feat2diff), outliers)]

    ##exclude 0
    zero.len = names(feat2diff)[which(feat2diff == 0)]
    feat2diff = feat2diff[setdiff(names(feat2diff), zero.len)]
    
    ##hist
    hist(log10(feat2diff))

    ##ggplot hist
    library('ggplot2')
    feat2diff.mat = as.data.frame(as.matrix(feat2diff))
    colnames(feat2diff.mat) = 'neigh.dist'
    
    gg = ggplot(feat2diff.mat, aes_string(x = 'neigh.dist'))
    ##gg = gg + geom_histogram(bins = 30)    
    gg = gg + geom_density()
    gg = gg + scale_x_log10(breaks = c(0, 10, 100, 1000, 10000, 100000), labels = c(0, 10, 100, 1000, 10000, 100000))

    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))
    gg = gg + xlab('Distance between neighboring heterozygous variants within genes')
    gg = gg + ylab('Number of variants')
    gg = gg + theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5))

    j.pdf = file.path(out_pdf_dir, 'neighvar2dist.dna.dens.pdf')
    pdf(j.pdf)
    plot(gg)
    dev.off()


    ##*################
    ##Distances at CDS level
    ##*################
    
    ##http://finzi.psych.upenn.edu/library/GenomicFeatures/html/coordinate-mapping-methods.html
    ##source("https://bioconductor.org/biocLite.R")
    ##biocLite('GenomicFeatures')
    ##biocLite('TxDb.Hsapiens.UCSC.hg19.knownGene')    
    library("GenomicFeatures")

    ##gene ranges
    ##biocLite('TxDb.Hsapiens.UCSC.hg38.knownGene')
    ##library('TxDb.Hsapiens.UCSC.hg19.knownGene')
    library('TxDb.Hsapiens.UCSC.hg38.knownGene')
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
    hg38.genes.gr = genes(txdb)
    hg38.tx.gr = transcriptsBy(txdb, 'gene')
    hg38.cds.gr = cdsBy(txdb, 'gene')
    length(hg38.cds.gr) #19,741
    
    ##rsid -> genomic coordinate
    ##biocLite('SNPlocs.Hsapiens.dbSNP141.GRCh38')
    library('SNPlocs.Hsapiens.dbSNP141.GRCh38')
    rsids = c("rs4803846", "rs139322374", "rs7250736", "rs7250754", "rs9749185")
    rsids = featdata[, 'var']
    length(rsids) #227,663
    snps.gr = snpid2grange(SNPlocs.Hsapiens.dbSNP141.GRCh38, rsids)
    snps.gr = snpsById(SNPlocs.Hsapiens.dbSNP141.GRCh38, rsids, ifnotfound = "drop")
    length(snps.gr) #226,491

    ##fix annot styles
    seqlevelsStyle(snps.gr)
    seqlevelsStyle(hg38.genes.gr)
    genome(snps.gr)
    genome(hg38.genes.gr)
    
    seqlevelsStyle(snps.gr) = seqlevelsStyle(hg38.genes.gr)
    genome(snps.gr) = genome(hg38.genes.gr)
    
    ##genomic coord -> gene coord
    snp2cds.gr = mapToTranscripts(snps.gr, transcripts = unlist(hg38.cds.gr)) #extractor.fun = GenomicFeatures::genes
    snp2tx.gr = mapToTranscripts(snps.gr, transcripts = unlist(hg38.tx.gr)) #extractor.fun = GenomicFeatures::genes
    length(snp2cds.gr) #2,791
    
    ##several transcripts per gene. Best would be to use all exons and UTRs per gene
    ##a simpler but non-perfect solution to only choose one transcript per gene (the one with the most snp hits)
    snp2tx.df = as.data.frame(snp2tx.gr)
    snp2tx.df = as.data.frame(snp2cds.gr)
    snp2tx.list = split(snp2tx.df, snp2tx.df[, 'transcriptsHits'])
    nrow(snp2tx.df) #2,791
    length(unique(snp2tx.df[, 'xHits'])) #2,368
    length(snp2tx.list) #19,377; 2,304
    
    ##filter out tx with only one hit
    tx2nhits = unlist(lapply(snp2tx.list, nrow))
    hist(log10(tx2nhits))
    pass.tx = names(tx2nhits)[which(tx2nhits >= 2)]
    snp2tx.list = snp2tx.list[pass.tx] #18,386
    length(snp2tx.list) #305
    
    feat2diff = unlist(lapply(snp2tx.list, function(j.feat){diff(sort(j.feat[, 'start']))}))

    ##exclude outliers (two)
    max(feat2diff) #3,037
    max.dist = 1e6
    outliers = names(feat2diff)[which(feat2diff > max.dist)]
    outliers #0
    feat2diff = feat2diff[setdiff(names(feat2diff), outliers)]

    ##exclude 0
    zero.len = names(feat2diff)[which(feat2diff == 0)]
    length(zero.len) #0
    feat2diff = feat2diff[setdiff(names(feat2diff), zero.len)]
    
    ##hist
    summary(feat2diff)
    hist(log10(feat2diff))

    ##ggplot hist
    library('ggplot2')
    feat2diff.df = as.data.frame(as.matrix(feat2diff))
    colnames(feat2diff.df) = 'neigh.dist'
    neigh.dist.log10 = log10(feat2diff.df[, 'neigh.dist'])
    feat2diff.df = cbind(feat2diff.df, neigh.dist.log10)        
    
    gg = ggplot(feat2diff.df, aes_string(x = 'neigh.dist'))
    ##gg = gg + geom_histogram(bins = 12)
    gg = gg + scale_x_log10(breaks = c(0, 10, 100, 1000, 10000), labels = c(0, 10, 100, 1000, 10000))
    gg = gg + geom_density()
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))
    gg = gg + xlab('Distance between neighboring heterozygous variants within genes')
    gg = gg + ylab('Number of variants')
    gg = gg + theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5))

    j.pdf = file.path(out_pdf_dir, 'neighvar2dist.cds.dens.pdf')
    pdf(j.pdf)
    plot(gg)
    dev.off()

}
