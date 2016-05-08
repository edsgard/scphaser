


##Params
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

acfiles.dir = file.path(cloud.dir, 'projects/scphaser/nogit/scripts/borel')
snpdir = '/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009/dnaseq/vcf'


##*###
##Files
##*###

##IN
snp.meta.rds = file.path(snpdir, 'UCF1014.snv.het.dp.gene.dbsnp.rds')
acfiles.txt = file.path(acfiles.dir, 'ac.files')

##OUT
ref.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', 'ref.counts.rds')
alt.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', 'alt.counts.rds')
snp.annot.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/borel', 'snp.annot.rds')

main <- function(){

    ##*###
    ##Read data
    ##*###
    
    ##snp info
    if(0){ #Once and for all
        snp.tab = file.path(snpdir, 'UCF1014.bcf.vcf.snv.het.dp.gene.dbsnp')
        snps = read.table(snp.tab, sep = '\t', stringsAsFactors = FALSE)
        colnames(snps) = c('chr.pos.ref.alt', 'chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 's1', 'rsid', 'gene', 'gene.annot')
        chr.pos = gsub(' ' , '', apply(snps[, c('chr', 'pos')], 1, paste, collapse = '.', sep = ''))
        snps = cbind(snps, chr.pos, stringsAsFactors = FALSE)
        saveRDS(snps, file = snp.meta.rds)
    }
    snps = readRDS(snp.meta.rds)
    
    ##allele counts
    ac.files = read.table(acfiles.txt)[, 1]
    ac.list = lapply(ac.files, function(j.file){
        ac.mat = try(read.table(j.file, sep = '\t', stringsAsFactors = FALSE), silent = TRUE)
        if(is.data.frame(ac.mat)){
            colnames(ac.mat) = c('chr', 'pos', 'A', 'C', 'G', 'T')
        }
        return(ac.mat)
    })
    names(ac.list) = ac.files
    mat2class = unlist(lapply(ac.list, class))
    empty.files = names(mat2class)[which(mat2class == 'try-error')]
    pass.files = setdiff(names(mat2class), empty.files)
    ac.list = ac.list[pass.files]
    
    
    ##*###
    ##create matrix
    ##*###
    ac.files = names(ac.list)
    base.files = basename(ac.files)
    sample2file = t(as.data.frame(strsplit(base.files, '\\.')))
    sample2file = cbind(sample2file, ac.files)
    sample2file.list = tapply(sample2file[, 'ac.files'], sample2file[, 1], unique)

    samples = names(sample2file.list)
    sample2mat.list = list()
    for(j.sample in samples){
        print(j.sample)
        sample.files = sample2file.list[[j.sample]]

        sample.mat = do.call(rbind, ac.list[sample.files])
        chr.pos = gsub(' ', '', apply(sample.mat[, c('chr', 'pos')], 1, paste, collapse = '.', sep = ''))        
        rownames(sample.mat) = chr.pos

        sample2mat.list[[j.sample]] = sample.mat
    }

    
    ##*###
    ##Get annot (ref, alt, rsid)
    ##*###
    sample2mat.annot.list = lapply(sample2mat.list, merge, y = snps, by.x = 'row.names', by.y = 'chr.pos')

    
    ##*###
    ##Get ref and alt counts
    ##*###
    samples = names(sample2mat.annot.list)
    n.samples = length(samples) #163
    rsids = unique(unlist(lapply(sample2mat.annot.list, '[[', 'rsid')))
    n.rsid = length(rsids) #231,075
    
    ref.mat = matrix(0, nrow = n.rsid, ncol = n.samples, dimnames = list(rsids, samples))
    alt.mat = ref.mat
    
    for(j.sample in samples){
        print(j.sample)
        j.mat = sample2mat.annot.list[[j.sample]]
        j.rsid = j.mat[, 'rsid']
        
        j.ref = as.integer(apply(j.mat, 1, function(j.var, allele.col){j.allele = j.var[allele.col]; j.var[j.allele];}, allele.col = 'ref'))        
        ref.mat[j.rsid, j.sample] = j.ref

        j.alt = as.integer(apply(j.mat, 1, function(j.var, allele.col){j.allele = j.var[allele.col]; j.var[j.allele];}, allele.col = 'alt'))
        alt.mat[j.rsid, j.sample] = j.alt
    }

    
    ##Dump
    saveRDS(ref.mat, file = ref.counts.rds)
    saveRDS(alt.mat, file = alt.counts.rds)


    ##*###
    ##snp info
    ##*###
    snp.filt = snps
    rownames(snp.filt) = snp.filt[, 'rsid']
    var = snp.filt[, 'rsid']
    feat = snp.filt[, 'gene']
    feat = sub('\\(.*', '', feat)
    snp.filt = cbind(snp.filt, var, feat, stringsAsFactors = FALSE)
    snp.filt = snp.filt[, c('feat', 'var', 'chr', 'pos', 'ref', 'alt', 'gene.annot')]

    snp.filt = snp.filt[rsids, ]
    
    saveRDS(snp.filt, file = snp.annot.rds)
    
}

test <- function(){
    ##simplify colnames
    if(0)
    a = sub('_.*', '', sub('_EGAR.*DEMO[0-9]_', '', samples))
    length(unique(a)) #93
    samples[grep('C77', samples)]
    ##"_EGAR00001248842_DEMO4_C77_130508_3" "_EGAR00001248843_DEMO4_C77_130508_6"
    ##probably technical replicate seqed twice
}
