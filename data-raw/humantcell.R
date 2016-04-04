

library('devtools')
devtools::load_all()
source('data-raw/read_rpkmfiles.R')


##*###
##Read data
##*###
fem.Tcell <- read_rpkmfiles("../nogit/data/tcells/oct21/rpkmforgenes/allelichits/female_P1292_P1316_YFVNM_passedqc.fake_mouse_allelehits.txt")
fem.Tcell$c57.counts   <- fem.Tcell$counts[ , 1:(length(fem.Tcell$samples)/2)]
fem.Tcell$cast.counts  <- fem.Tcell$counts[ , (length(fem.Tcell$samples)/2+1):length(fem.Tcell$samples)]

##feature data
featdata = as.data.frame(cbind(fem.Tcell[['genes']], fem.Tcell[['tx']]), stringsAsFactors = FALSE)
colnames(featdata) = c('feat', 'var')
chr = sub('chr', '', sub('(^.*):.*', '\\1', featdata[, 'var']))
featdata = cbind(featdata, chr, stringsAsFactors = FALSE)
lapply(featdata, class)

##phenodata
samples = fem.Tcell[['samples']]
samples = sub('_c57only$', '', samples[grep('_c57only$', samples)])
phenodata = as.data.frame(samples, stringsAsFactors = FALSE)
colnames(phenodata) = 'sample'
rownames(phenodata) = phenodata[, 'sample']

##count matrixes
alt.counts = fem.Tcell[['c57.counts']]
rownames(alt.counts) = featdata[, 'var']
colnames(alt.counts) = samples
ref.counts = fem.Tcell[['cast.counts']]
rownames(ref.counts) = featdata[, 'var']
colnames(ref.counts) = samples

##fix duplicates of variant names
pass.ind = which(!duplicated(featdata[, 'var']))
alt.counts = alt.counts[pass.ind, ]
ref.counts = ref.counts[pass.ind, ]
featdata = featdata[pass.ind, ]
rownames(featdata) = featdata[, 'var']


##*###
##Filter
##*###
acset = list(featdata = featdata, refcount = ref.counts, altcount = alt.counts)

##filter vars on expression in at least n cells
ncells = 1
mincount = 3
acset = filter_var_mincount(acset, mincount, ncells)

##filter feats on number of vars
nmin_var = 2
acset = filter_feat_nminvar(acset, nmin_var)

##filter feats on counts. Require feat to have at least two vars that have expression in at least three cells, among cells that express at least two vars.
mincount = 3
ncells = 3
acset = filter_feat_counts(acset, mincount, ncells)

##filter feats on gt. Require feat to have at least two vars that have monoallelic calls in at least three cells, among cells that have at least two vars with monoallelic calls.
ncells = 3
acset = call_gt(acset)
acset = filter_feat_gt(acset, ncells)

##filter genes with variants on different chromosomes
featdata = acset[['featdata']]
feat2chr = tapply(featdata[, 'chr'], featdata[, 'feat'], table)
feat2nchr = unlist(lapply(feat2chr, length))
pass_feat = names(feat2nchr)[which(feat2nchr == 1)]
acset = subset_feat(acset, pass_feat)
lapply(acset, dim) #872, 505
length(unique(acset[['featdata']][, 'feat'])) #292


##ensure order is the same in pheno as in count matrixes
refcount = acset[['refcount']]
altcount = acset[['altcount']]
nrow(phenodata) == ncol(refcount)
nrow(phenodata) == length(which(phenodata[, 'sample'] == colnames(refcount)))
nrow(phenodata) == length(which(colnames(altcount) == colnames(refcount)))
phenodata = phenodata[colnames(refcount), , drop = FALSE]


##*###
##Dump
##*###

##store in list
humantcell = c(acset[c('featdata', 'refcount', 'altcount')], list(phenodata = phenodata))

##dump
devtools::use_data(humantcell, overwrite = TRUE)

##check compression algorithm
tools::checkRdaFiles('data')
