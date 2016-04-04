

##Dirs
cloud.dir = '/Volumes/Data/cloud/btsync/work/rspd'
datadir = file.path(cloud.dir, 'projects/scphaser/nogit/data/marinov')
db.dir = file.path(cloud.dir, 'data/external/sra')
scriptdir = file.path(cloud.dir, 'projects/scphaser/nogit/scripts/marinov')

##*###
##Files
##*###
##IN
meta.tab = file.path(datadir, 'GSE44618_series_matrix.stripped.tab')
sra.txt = file.path(scriptdir, 'sra.files')
libid.txt = file.path(datadir, 'singlecell.15.libid')

##OUT
pass.exp.txt = file.path(datadir, 'singlecell.sra.files')

##Libs
library(SRAdb)

main <- function(){

    ##read data
    meta.df = read.table(meta.tab, sep = '\t', stringsAsFactors = FALSE, header = FALSE)
    meta.df = t(meta.df)
    libids = read.table(libid.txt, stringsAsFactors = FALSE, header = FALSE)
    
    ##fix cols
    cols = meta.df[1, ]
    cols[11] = '!Sample_characteristics_ch1.2'
    cols[33] = '!Sample_relation.2'
    colnames(meta.df) = cols
    meta.df = meta.df[-1, ]
    cols = colnames(meta.df)
    cols = sub('^!Sample_', '', cols)
    cols = sub('title', 'id', cols)
    colnames(meta.df) = cols
    meta.df[, 'id']
    rownames(meta.df) = meta.df[, 'id']
    
    ##set labexpid to rownames
    labexpid = meta.df[, 'characteristics_ch1']
    labexpid = sub('labexpid: ', '', labexpid)
    colnames(meta.df)[which(colnames(meta.df) == 'characteristics_ch1')] = 'labexpid'
    meta.df[, 'labexpid'] = labexpid
    rownames(meta.df) = meta.df[, 'labexpid']

    ##single-cell labexpids of the 15 cells used in the paper
    libids = sub('-.*', '', libids[, 1])
    meta.single = meta.df[libids, ]
    
    ##get srx for what I believed was single-cell samples
    rownames(meta.df) = meta.df[, 'id']
    id = meta.df[, 'id']
    id.pass = setdiff(id, id[grep('pool_split', id)])
    id.pass = setdiff(id.pass, id.pass[grep('cells', id.pass)])
    id.pass = setdiff(id.pass, id.pass[grep('0ng_', id.pass)])
    id.pass = setdiff(id.pass, id.pass[grep('0pg_', id.pass)])
    id.pass = setdiff(id.pass, id.pass[grep('0[AB]', id.pass)])
    
    meta.filt = meta.df[id.pass, ]
    srx = basename(meta.filt[, 'supplementary_file_2'])

    ##check
    meta.diff = meta.filt[setdiff(id.pass, meta.single[, 'id']), ]
    meta.diff[, 'id']
    
    ##If using the subset of 15 cells among the 28 single cells
    ##srx = basename(meta.single[, 'supplementary_file_2'])

    
    ##*########################
    ##SRX to SRR mapping
    ##*#########################
    
    ##SRR ids are not in the meta-file
    grep('SRR', meta.df[1, ])
    
    ##Get DB
    db.file = file.path(db.dir, paste('2016-03-04', "SRAmetadb.sqlite", sep= '.'))
    if(0){ #download SRA meta-db
        db.file = file.path(db.dir, paste(Sys.Date(), "SRAmetadb.sqlite.gz", sep= '.'))
        sqlfile = getSRAdbFile(destfile = db.file)
    }
    sra_con = dbConnect(SQLite(), db.file)
    
    ##Get run accession for exp accession
    srx.str = paste("'", srx, "',", sep = '', collapse = '')
    srx.str = sub(',$', '', srx.str)
    query.str = sprintf("select run_accession, experiment_accession from sra where experiment_accession in (%s)", srx.str)
    run2exp = dbGetQuery(sra_con, query.str)
    srr.pass = run2exp[, 'run_accession']
    
    sra = read.table(sra.txt, stringsAsFactors = FALSE)
    srr = basename(dirname(sra[, 1]))
    rownames(sra) = srr
    sra.pass = sra[srr.pass, ]
    
    ##Dump run accessions
    write.table(sra.pass, quote = F, col.names = F, row.names = F, file = pass.exp.txt)
    
    
}
