

library('devtools')
source('data-raw/read_rpkmfiles.R')


##*###
##Read data
##*###
fem.Tcell <- read_rpkmfiles("../../../../nogit/data/tcells/oct21/rpkmforgenes/allelichits/female_P1292_P1316_YFVNM_passedqc.fake_mouse_allelehits.txt")
fem.Tcell$c57.counts   <- fem.Tcell$counts[ , 1:(length(fem.Tcell$samples)/2)]
fem.Tcell$cast.counts  <- fem.Tcell$counts[ , (length(fem.Tcell$samples)/2+1):length(fem.Tcell$samples)]

##feature data
featdata = as.data.frame(cbind(fem.Tcell[['genes']], fem.Tcell[['tx']]), stringsAsFactors = FALSE)
colnames(featdata) = c('feat', 'var')

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

##store in list
humantcell = list(featdata = featdata, phenodata = phenodata, refcount = ref.counts, altcount = alt.counts)

##dump
devtools::use_data(humantcell)

##check compression algorithm
tools::checkRdaFiles('data')
