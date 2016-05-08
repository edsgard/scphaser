

##allele counts: /mnt/kauffman/danielr/crick/Xandclones_BR/bulksamples_allelecounts/cast_c57/cellsums/Bf37_mRNA_Fibroblast_Total_fraction
##allele list
##SNP list: validated_cast_c57_snps.txt
##according to code in snp_stats2.py the columns should correspond to: c57, cast

sys = 'rs13'

##*###
##Dirs
##*###
if(sys == 'rs13'){
    cloud.dir = '/Volumes/Data/cloud/btsync/work/rspd'
}
data.dir = file.path(cloud.dir, 'projects/scphaser/nogit/data/mousehybrid/bulk')

##*###
##Files
##*###
val.snps.tab = file.path(data.dir, 'validated_cast_c57_snps.txt')
ac.tab = file.path(data.dir, 'Bf37_mRNA_Fibroblast_Total_fraction')

main <- function(){

    ##*###
    ##Read data
    ##*###
    val.snps = read.table(val.snps.tab, sep = '\t', stringsAsFactors = FALSE)
    colnames(val.snps) = c('chr', 'pos', 'c57', 'cast', 'd1', 'd2', 'd3', 'd4', 'd5', 'd6')
    ac = read.table(ac.tab, sep = '\t')
    colnames(ac) = c('chr', 'pos', 'A', 'C', 'G', 'T')
    nrow(ac) #507,338

    ##Filter ac counts on val.snps
    val.ac = merge(val.snps[, 1:4], ac, by = c('chr', 'pos')) 
    dim(val.ac) #203,449

    ##Get allele counts
    hap.ac = t(apply(val.ac, 1, function(j.row){c57.count = j.row[j.row['c57']]; cast.count = j.row[j.row['cast']]; new.row = c(j.row, c57.count, cast.count); return(new.row)}))
    colnames(hap.ac) = c(colnames(val.ac), c('c57.ac', 'cast.ac'))
    hap.ac = as.data.frame(hap.ac, stringsAsFactors = FALSE)
    num.cols = c('A', 'T', 'G', 'C', 'c57.ac', 'cast.ac')
    hap.ac[num.cols] = lapply(hap.ac[num.cols], as.integer)

    tot.ac = apply(hap.ac[, c('cast.ac', 'c57.ac')], 1, sum)
    min.count = 3
    hap.filt = hap.ac[which(tot.ac >= min.count), ]
    nrow(hap.filt) #123,318
    ase = hap.filt[, 'cast.ac'] / (hap.filt[, 'cast.ac'] + hap.filt[, 'c57.ac'])

    pdf(file = 'bulk.ase.pdf')
    plot(density(ase))
    dev.off()
}

