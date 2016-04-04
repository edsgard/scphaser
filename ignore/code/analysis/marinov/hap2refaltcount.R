


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


##*###
##Files
##*###

##IN
snp.meta.rds = file.path(cloud.dir, 'data/external/Marinov_GenRes_2014', 'hg19.ref2mat2pat.het.dbsnp.refseq.1based.rds')
mat.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov', paste('mat', '.counts.rds', sep = ''))
pat.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov', paste('pat', '.counts.rds', sep = ''))

##OUT
ref.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov', paste('ref', '.counts.rds', sep = ''))
alt.counts.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov', paste('alt', '.counts.rds', sep = ''))
snp.annot.rds = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov', 'snp.annot.rds')

main <- function(){

    ##*###
    ##Read data
    ##*###
    
    snp.annot = readRDS(snp.meta.rds)
    mat.counts = readRDS(mat.counts.rds)
    pat.counts = readRDS(pat.counts.rds)

    ##*###
    ##Create mat, pat matrices with the same vars
    ##*###
    mat.rsid = rownames(mat.counts)
    pat.rsid = rownames(pat.counts)
    all.rsid = unique(c(mat.rsid, pat.rsid))
    n.vars = length(all.rsid)
    samples = colnames(mat.counts)
    n.samples = length(samples)
    new.mat.counts = matrix(0, nrow = n.vars, ncol = n.samples, dimnames = list(all.rsid, samples))
    new.pat.counts = matrix(0, nrow = n.vars, ncol = n.samples, dimnames = list(all.rsid, samples))
    
    new.mat.counts[mat.rsid, ] = mat.counts
    new.pat.counts[pat.rsid, ] = pat.counts

    saveRDS(new.mat.counts, mat.counts.rds)
    saveRDS(new.pat.counts, pat.counts.rds)

    
    ##*###
    ##
    ##*###
    ref.counts = matrix(0, nrow = n.vars, ncol = n.samples, dimnames = list(all.rsid, samples))
    alt.counts = matrix(0, nrow = n.vars, ncol = n.samples, dimnames = list(all.rsid, samples))

    rownames(snp.annot) = snp.annot[, 'rsid']
    snp.filt = snp.annot[all.rsid, ]
    ref.ismat = snp.filt[which(snp.filt[, 'ref'] == snp.filt[, 'mat.allele']), 'rsid']
    ref.ispat = setdiff(all.rsid, ref.ismat)
    alt.ismat = ref.ispat
    alt.ispat = ref.ismat
    ref.counts[ref.ismat, ] = new.mat.counts[ref.ismat, ]
    ref.counts[ref.ispat, ] = new.pat.counts[ref.ispat, ]

    alt.counts[alt.ismat, ] = new.mat.counts[alt.ismat, ] 
    alt.counts[alt.ispat, ] = new.pat.counts[alt.ispat, ]
    
    saveRDS(ref.counts, file = ref.counts.rds)
    saveRDS(alt.counts, file = alt.counts.rds)
    saveRDS(snp.filt, file = snp.annot.rds)
}
