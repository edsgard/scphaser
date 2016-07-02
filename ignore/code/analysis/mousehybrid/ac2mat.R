


##Params
sys = 'dna'
n.downsamp.vars = 5000000

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
acfiles.dir = file.path(cloud.dir, 'projects/scphaser/nogit/scripts/mousehybrid')
snpdir = '/mnt/kauffman/edsgard/rsync/work/rspd/data/external/sanger_mouse_proj'

##OUT
ac.dir = file.path(cloud.dir, sprintf('projects/scphaser/nogit/data/mousehybrid/all_snps/downsampled_N%i', n.downsamp.vars))


##*###
##Files
##*###

##IN
snp.meta.rds = file.path(snpdir, 'mm9.cast.snps.refseq.rds')
acfiles.txt = file.path(acfiles.dir, sprintf('ac.N_%i.files', n.downsamp.vars))

##OUT
ref.counts.rds = file.path(ac.dir, 'c57.counts.rds')
alt.counts.rds = file.path(ac.dir, 'cast.counts.rds')
snp.annot.rds = file.path(ac.dir, 'snp.annot.rds')

main <- function(){

    ##*###
    ##Read data
    ##*###
    
    ##snp info
    if(0){ #Once and for all
        snp.tab = file.path(snpdir, 'mm9.snps.refseqgenes') ##see log.sh
        snps = read.table(snp.tab, sep = '\t', stringsAsFactors = FALSE)
        colnames(snps) = c('gene.annot', 'chr', 'pos', 'end', 'c57', 'cast', 'gene')
        chr.pos = gsub(' ' , '', apply(snps[, c('chr', 'pos')], 1, paste, collapse = '.', sep = ''))
        chr.pos.ref.alt = gsub(' ' , '', apply(snps[, c('chr', 'pos', 'c57', 'cast')], 1, paste, collapse = '.', sep = ''))
        snps = cbind(snps, chr.pos, chr.pos.ref.alt, stringsAsFactors = FALSE)
        rownames(snps) = chr.pos.ref.alt
        saveRDS(snps, file = snp.meta.rds)
    }
    snps = readRDS(snp.meta.rds)
    
    ##read allele counts
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
    names(ac.list) = basename(names(ac.list))
    
    ##*###
    ##add chr.pos
    ##*###
    ac.list = lapply(ac.list, function(sample.mat){chr.pos = gsub(' ', '', apply(sample.mat[, c('chr', 'pos')], 1, paste, collapse = '.', sep = '')); rownames(sample.mat) = chr.pos; return(sample.mat);})

    
    ##*###
    ##Get annot (ref, alt, rsid)
    ##*###
    sample2mat.annot.list = lapply(ac.list, merge, y = snps, by.x = 'row.names', by.y = 'chr.pos')
    sample2mat.annot.list = lapply(sample2mat.annot.list, function(j.mat){colnames(j.mat) = sub('Row.names', 'chr.pos', colnames(j.mat)); colnames(j.mat) = sub('chr.x', 'chr', colnames(j.mat)); colnames(j.mat) = sub('pos.x', 'pos', colnames(j.mat)); return(j.mat);})
    
    ##*###
    ##Get allele counts for the haplotype alleles
    ##*###
    samples = names(sample2mat.annot.list)
    n.samples = length(samples)
    n.samples #336
    
    snp.pos = unique(unlist(lapply(sample2mat.annot.list, '[[', 'chr.pos.ref.alt')))
    n.snps = length(snp.pos) 
    n.snps #5M: 1,990,209; 2.5M: 995,196; 1M: 398,399; 785k: 314,036
    ref.mat = matrix(0, nrow = n.snps, ncol = n.samples, dimnames = list(snp.pos, samples))
    alt.mat = ref.mat
    
    for(j.sample in samples){
        print(j.sample)
        j.mat = sample2mat.annot.list[[j.sample]]
        j.rsid = j.mat[, 'chr.pos.ref.alt']
        
        j.ref = as.integer(apply(j.mat, 1, function(j.var, allele.col){j.allele = j.var[allele.col]; j.var[j.allele];}, allele.col = 'c57'))        
        ref.mat[j.rsid, j.sample] = j.ref

        j.alt = as.integer(apply(j.mat, 1, function(j.var, allele.col){j.allele = j.var[allele.col]; j.var[j.allele];}, allele.col = 'cast'))
        alt.mat[j.rsid, j.sample] = j.alt
    }

    
    ##Dump
    dir.create(dirname(ref.counts.rds), recursive = TRUE)
    saveRDS(ref.mat, file = ref.counts.rds)
    saveRDS(alt.mat, file = alt.counts.rds)


    ##*###
    ##snp info
    ##*###
    snp.filt = snps    
    var = snp.filt[, 'chr.pos.ref.alt']
    feat = snp.filt[, 'gene']
    feat = sub('\\(.*', '', feat)
    snp.filt = cbind(snp.filt, var, feat, stringsAsFactors = FALSE)
    snp.filt = snp.filt[, c('feat', 'var', 'chr', 'pos', 'c57', 'cast', 'gene.annot')]
    colnames(snp.filt) = c('feat', 'var', 'chr', 'pos', 'ref', 'alt', 'gene.annot')

    ##filter
    snp.filt = snp.filt[snp.pos, ]
    dim(snp.filt) #314,036
    
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
