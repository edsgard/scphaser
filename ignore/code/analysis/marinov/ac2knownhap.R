


##Params
sys = 'dna'
hap = 'pat'


##*###
##Dirs
##*###

if(sys == 'dna'){
    cloud.dir = '/mnt/kauffman/edsgard/cloud/btsync/work/rspd'
}
if(sys == 'rs13'){
    cloud.dir = '/Volumes/Data/cloud/btsync/work/rspd'
}

acfiles.dir = file.path(cloud.dir, 'projects/scphaser/nogit/scripts/marinov')


##*###
##Files
##*###

##IN
snp.meta.rds = file.path(cloud.dir, 'data/external/Marinov_GenRes_2014', 'hg19.ref2mat2pat.het.dbsnp.refseq.1based.rds')
acfiles.txt = file.path(acfiles.dir, paste(hap, 'ac.files', sep = '.'))

##OUT
hap.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov', paste(hap, '.counts.rds', sep = ''))


main <- function(){

    ##*###
    ##Read data
    ##*###
    
    ##snp info
    if(0){ #Once and for all
        snp.dir = file.path('/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/vars')
        snp.txt = file.path(snp.dir, 'ref2mat2pat.het.dbsnp.refseq.hap.1based') ##pileup is 1-based: http://samtools.sourceforge.net/pileup.shtml
        snps = read.table(snp.txt, sep = '\t', stringsAsFactors = FALSE)
        colnames(snps) = c('chr.pos', 'chr', 'ref.pos', 'ref', 'alt', 'rsid', 'gene', 'gene.annot', 'mat.pos', 'pat.pos', 'mat.allele', 'pat.allele')
        saveRDS(snps, file = snp.meta.rds)
    }
    snps = readRDS(snp.meta.rds)
    
    ##allele counts
    ac.files = read.table(acfiles.txt)[, 1]
    ac.list = lapply(ac.files, function(j.file){
        ac.mat = read.table(j.file, sep = '\t', stringsAsFactors = FALSE)
        colnames(ac.mat) = c('chr', 'pos', 'A', 'C', 'G', 'T')
        return(ac.mat)
    })
    names(ac.list) = ac.files

    
    ##*###
    ##get hap counts
    ##*###
    files = names(ac.list)
    hap.ac.list = list()
    j.file.it = 0
    for(j.file in files){
        j.file.it = j.file.it + 1
        print(j.file.it)
        hap.ac.list[[j.file]] = ntcount2hapcount(ac.list[[j.file]], snps = snps)
    }
    ##status: fin

    
    ##*###
    ##create matrix
    ##*###
    base.files = basename(ac.files)
    sample2file = t(as.data.frame(strsplit(base.files, '\\.')))
    sample2file = cbind(sample2file, ac.files)
    sample2file.list = tapply(sample2file[, 'ac.files'], sample2file[, 1], unique)

    samples = names(sample2file.list)
    sample2mat.list = list()
    for(j.sample in samples){
        sample.files = sample2file.list[[j.sample]]

        sample.mat = do.call(rbind, hap.ac.list[sample.files])
        rownames(sample.mat) = sample.mat[, 'rsid']

        sample2mat.list[[j.sample]] = sample.mat
    }

    ##Haplotype count matrix
    n.samples = length(samples)
    all.rsid = unique(unlist(lapply(sample2mat.list, '[[', 'rsid')))
    n.rsid = length(all.rsid) #29336
    all.genes = unique(unlist(lapply(sample2mat.list, '[[', 'gene')))
    n.genes = length(all.genes) #8027

    hap.counts.mat = matrix(nrow = n.rsid, ncol = n.samples, dimnames = list(all.rsid, samples))
    hapcount.col = paste(hap, 'count', sep = '.')
    for(j.sample in samples){
        j.mat = sample2mat.list[[j.sample]]
        j.rsid = j.mat[, 'rsid']
        hap.counts.mat[j.rsid, j.sample] = j.mat[, hapcount.col]
    }

    ##set NA to 0
    hap.counts.mat[is.na(hap.counts.mat)] = 0
    
    ##Dump
    saveRDS(hap.counts.mat, file = hap.counts.rds)
}

ntcount2hapcount <- function(ac.mat, snps){
    hap2ac.mat = merge(snps, ac.mat, by.x = c('chr', paste(hap, '.pos', sep = '')), by.y = c('chr', 'pos'))

    hap.col = paste(hap, '.allele', sep = '')
    hap.ac = as.integer(apply(hap2ac.mat, 1, function(j.var, hap.col){hap.allele = j.var[hap.col]; j.var[hap.allele];}, hap.col = hap.col))
    hap2ac.mat = cbind(hap2ac.mat, hap.ac)
    colnames(hap2ac.mat)[ncol(hap2ac.mat)] = paste(hap, 'count', sep = '.')

    return(hap2ac.mat)
}
