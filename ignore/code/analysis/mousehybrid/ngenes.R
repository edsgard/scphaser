

snp.dir = '/mnt/kauffman/edsgard/rsync/work/rspd/data/external/sanger_mouse_proj'
snp.file = file.path(snp.dir, '/split/mm9.snps.refseq.nointrons')

snp.file = file.path(snp.dir, '/split/mm9.snps.refseqgenes')

main <- function(){

    snps = read.table(snp.file, sep = '\t', stringsAsFactors = FALSE)
    colnames(snps) = c('gene.feat', 'chr', 'start', 'end', 'ref', 'alt', 'gene')
    nrow(snps) #7,834,497; intronic_no: 307,991
    
    gene2nvars = table(snps[, 'gene'])
    genes.pass = names(gene2nvars)[which(gene2nvars >= 2)]
    length(genes.pass) #22,181; intronic_no: 20,268

    snp.filt = merge(snps, as.matrix(genes.pass), by.x = 'gene', by.y = 1)

    nrow(snp.filt) #7,833,787; intronic_no: 306,574
    
    
}
