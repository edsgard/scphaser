

sys = 'dna'


##*###
##Dirs
##*###

if(sys == 'dna'){
    cloud.dir = '/mnt/kauffman/edsgard/cloud/btsync/work/rspd'
}
if(sys == 'rs13'){
    cloud.dir = '/Volumes/Data/cloud/btsync/work/rspd'
}

snpdir = '/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009/dnaseq/vcf'


##*###
##Files
##*###

##IN
snp.meta.rds = file.path(snpdir, 'UCF1014.snv.het.dp.gene.dbsnp.rds') ##see ac2mat.R

main <- function(){

    snps = readRDS(snp.meta.rds)

    snp.filt = snps
    rownames(snp.filt) = snp.filt[, 'rsid']
    var = snp.filt[, 'rsid']
    feat = snp.filt[, 'gene']
    feat = sub('\\(.*', '', feat)
    snp.filt = cbind(snp.filt, var, feat, stringsAsFactors = FALSE)

    nrow(snp.filt) #787,883
    length(unique(feat)) ##18,150
    
    feat2nvars = table(feat)
    min.vars = 2
    pass.feat = names(feat2nvars)[which(feat2nvars >= min.vars)]
    n.feat = length(pass.feat)
    n.feat #15,556

    ##filter snps on being in genes with at least two vars
    snp.filt = merge(snp.filt, as.matrix(pass.feat), by.x = 'feat', by.y = 1)
    nrow(snp.filt)
    ##785,289

    785289 / 7833787 #10.0%
    
    ###
    ##Exclude intronic
    ###
    keep.annot = setdiff(unique(snp.filt[, 'gene.annot']), c('intronic', 'ncRNA_intronic'))
    snps.exonic = merge(snp.filt, as.matrix(keep.annot), by.x = 'gene.annot', by.y = 1)
    nrow(snps.exonic) #28,342

    feat = snps.exonic[, 'feat']
    length(unique(feat)) ##10,124
    
    feat2nvars = table(feat)
    min.vars = 2
    pass.feat = names(feat2nvars)[which(feat2nvars >= min.vars)]
    n.feat = length(pass.feat)
    n.feat #5,706

    ##filter snps on being in genes with at least two vars
    snps.exonic = merge(snps.exonic, as.matrix(pass.feat), by.x = 'feat', by.y = 1)
    nrow(snps.exonic)
    ##23,924

    23924 / 306574 #7.8%
    
}
