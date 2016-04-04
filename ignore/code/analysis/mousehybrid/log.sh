

#number of genes
wc -l /mnt/kauffman/edsgard/src/rnaseq/ase/danielr/annot/snps/mm9/cast/validated_mm9_refseq_snp2genes_nooverlap.txt
f='/mnt/kauffman/edsgard/src/rnaseq/ase/danielr/annot/snps/mm9/cast/validated_mm9_refseq_snp2genes_nooverlap.txt'
awk -F'\t' '$2 >1' $f | wc -l #14,088


##*###
##Annovar
##*###

f='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/data/mousehybrid/vars.twovargenes.anno'
scriptdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/scripts/mousehybrid'
snpdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/data/mousehybrid'

##*###
##Download annots
##*###
cd '/mnt/kauffman/edsgard/prg/annovar'

#download gene annot
annotate_variation.pl -downdb -buildver mm9 gene mousedb

#download genome fasta
annotate_variation.pl --buildver mm9 --downdb seq mousedb/mm9_seq

#genome fasta + gene annot -> gene fasta
retrieve_seq_from_fasta.pl mousedb/mm9_refGene.txt -seqdir mousedb/mm9_seq -format refGene -outfile mousedb/mm9_refGeneMrna.fa


#split input for parallel
cd $snpdir
mkdir split
cd split
split -n 70 $f mm9.snps.

cd $scriptdir
find ${snpdir}/split -name 'mm9.snps.*' >infiles.list

#Annovar
cat infiles.list | xargs -I% echo 'annotate_variation.pl -build mm9 --geneanno -dbtype refGene' % '/mnt/kauffman/edsgard/prg/annovar/mousedb/' >cmds.sh
screen -S snp.annovar
cat cmds.sh | parallel -n 1 -P 70
#status: fin

#Extract SNPs with gene annot
cd ${snpdir}/split
cat *.variant_function >mm9.snps.variant_function
wc -l mm9.snps.variant_function #174,772

awk -F'\t' '{print $1;}' mm9.snps.variant_function | sort -u
#downstream, exonic, exonic;splicing, intergenic, intronic, ncRNA_exonic, ncRNA_intronic, ncRNA_splicing, splicing, upstream, upstream;downstream, UTR3, UTR5, UTR5;UTR3
#keep: exonic, exonic;splicing, intronic, ncRNA_exonic, ncRNA_intronic, ncRNA_splicing, splicing, UTR3, UTR5, UTR5;UTR3
#rm: intergenic, downstream, upstream, upstream;downstream: As far as I remember 1kb up/down is included by annovar by default.
awk -F'\t' -v OFS='\t' '$1 != "intergenic" && $1 != "downstream" && $1 != "upstream" && $1 != "upstream;downstream" && $1 != "intronic" && $1 != "ncRNA_intronic" {print $0;}' mm9.snps.variant_function >mm9.snps.refseqgenes
wc -l mm9.snps.refseqgenes #145,287
cut -f2 mm9.snps.refseqgenes | sort -u | wc -l #91,262
cut -f2 mm9.snps.refseqgenes | sort -u | head #TODO: Need to strip parenthesis info from gene-name

awk -F'\t' -v OFS='\t' '{s=$2; sub(/\(.*/, "", s); print $0, s;}' mm9.snps.refseqgenes >mm9.snps.stripped.refseqgenes
cut -f8 mm9.snps.stripped.refseqgenes | sort -u | head
cut -f8 mm9.snps.stripped.refseqgenes | sort -u | wc -l #12,247
awk -F'\t' -v OFS='\t' '{print $1, $9, $3, $4, $5, $6, $7, $8;}' mm9.snps.stripped.refseqgenes | sort -u >mm9.snps.refseqgenes
wc -l mm9.snps.refseqgenes #
cut -f2 mm9.snps.refseqgenes | sort -u | wc -l #
