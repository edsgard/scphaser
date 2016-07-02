

#number of genes
wc -l /mnt/kauffman/edsgard/src/rnaseq/ase/danielr/annot/snps/mm9/cast/validated_mm9_refseq_snp2genes_nooverlap.txt
f='/mnt/kauffman/edsgard/src/rnaseq/ase/danielr/annot/snps/mm9/cast/validated_mm9_refseq_snp2genes_nooverlap.txt'
awk -F'\t' '$2 >1' $f | wc -l #14,088


##*###
##Annovar
##*###

f='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/data/mousehybrid/vars.twovargenes.anno' ##see data-raw/mousehybrid.R
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
cut -f8 mm9.snps.stripped.refseqgenes | sort -u | wc -l ##12,948 -> 12,247
awk -F'\t' -v OFS='\t' '{print $1, $9, $3, $4, $5, $6, $7, $8;}' mm9.snps.stripped.refseqgenes | sort -u >mm9.snps.refseqgenes
wc -l mm9.snps.refseqgenes #
cut -f2 mm9.snps.refseqgenes | sort -u | wc -l #




###
#Download CAST SNPs from SANGER
###
cd '/mnt/kauffman/edsgard/rsync/work/rspd/data/external/sanger_mouse_proj'
wget ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/CAST_EiJ.mgp.v5.snps.dbSNP142.vcf.gz
gunzip *.gz
vcf='CAST_EiJ.mgp.v5.snps.dbSNP142.vcf'
wc -l $vcf #22,613,883
#ref=ftp://ftp-mouse.sanger.ac.uk/ref/GRCm38_68.fa


###
#Filter
###
#homo (field: GT)
grep -v '^#' $vcf | awk -F'\t' '{split($10, a, ":"); if(a[1] == "1/1"){print $0};}' >homo.vcf
wc -l homo.vcf
#21,330,784

#pass
awk -F'\t' '$7=="PASS"' homo.vcf >homo.pass.vcf
wc -l homo.pass.vcf
#20,668,274

#gene (CSQ field: Ensembl 78)
awk -F'\t' '{split($8, a, ";"); split(a[3], b, "|"); print b[5];}' homo.pass.vcf | sort -u >csq.vals
#intergenic_variant
#upstream_gene_variant
awk -F'\t' '{split($8, a, ";"); split(a[3], b, "|"); if(b[5]!="intergenic_variant" && b[5]!="upstream_gene_variant"){print $0};}' homo.pass.vcf > homo.pass.gene.vcf
wc -l homo.pass.gene.vcf
#. Not sure I should rely on this annotation. Let's wait with this filter.

#possibly: dbSNP


###
#Liftover mm10 -> mm9
###
#m10 -> mm9
#liftOver oldFile map.chain newFile unMapped
wget 'http://hgdownload.cse.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm9.over.chain.gz'
#VCF is 1-based so remove 1 (CHAIN files are 0-based according to their def: http://genome.ucsc.edu/goldenPath/help/chain.html)
cat homo.pass.vcf | awk -F'\t' -v OFS='\t' '{print "chr"$1, $2-1, $2, NR;}' >mm10.homo.pass.0based.bed
liftOver mm10.homo.pass.0based.bed mm10ToMm9.over.chain mm9.homo.pass.0based.bed mm10_to_mm9.unmapped
wc -l mm9.homo.pass.0based.bed
#20,659,048

##add back annotation
cat homo.pass.vcf | awk -F'\t' -v OFS='\t' '{print $0, NR;}' | sort -k11,11 >file.tmp
mv file.tmp homo.pass.vcf

sort -k4,4 mm9.homo.pass.0based.bed >file.tmp
mv file.tmp mm9.homo.pass.0based.bed

join -1 4 -2 11 mm9.homo.pass.0based.bed homo.pass.vcf >file.tmp
mv file.tmp mm9.homo.pass.0based.bed
wc -l mm9.homo.pass.0based.bed
#20,659,048


###
#Annovar to get RefSeq annot
###
scriptdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/scripts/mousehybrid'
snpdir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/sanger_mouse_proj'

#Convert vcf to annovar format
#chromosome, start position, end position, the reference nucleotides and the observed nucleotides
#convert2annovar.pl -format vcf4 example/ex2.vcf > ex2.avinput
awk -v OFS='\t' '{print $2, $4, $4, $8, $9;}' mm9.homo.pass.0based.bed >mm9.snps.avformat

#split input for parallel
mkdir split
cd split
split -n 80 ../mm9.snps.avformat mm9.snps.
cd $scriptdir
find ${snpdir}/split -name 'mm9.snps.*' >infiles.list

#Annovar
cat infiles.list | xargs -I% echo 'annotate_variation.pl -build mm9 --geneanno -dbtype refGene' % '/mnt/kauffman/edsgard/prg/annovar/mousedb/' >cmds.sh
screen -S snp.annovar
cat cmds.sh | parallel -n 1 -P 80
#status: fin

#Extract SNPs with gene annot
cd ${snpdir}/split
cat *.variant_function >mm9.snps.variant_function
wc -l mm9.snps.variant_function #20,658,984

awk -F'\t' '{print $1;}' mm9.snps.variant_function | sort -u
#downstream, exonic, exonic;splicing, intergenic, intronic, ncRNA_exonic, ncRNA_intronic, ncRNA_splicing, splicing, upstream, upstream;downstream, UTR3, UTR5, UTR5;UTR3
#keep: exonic, exonic;splicing, intronic, ncRNA_exonic, ncRNA_intronic, ncRNA_splicing, splicing, UTR3, UTR5, UTR5;UTR3
#
#rm: intergenic, downstream, upstream, upstream;downstream: As far as I remember 1kb up/down is included by annovar by default.
awk -F'\t' -v OFS='\t' '$1 != "intergenic" && $1 != "downstream" && $1 != "upstream" && $1 != "upstream;downstream" {print $0;}' mm9.snps.variant_function >mm9.snps.refseqgenes
wc -l mm9.snps.refseqgenes #7,834,508
cut -f2 mm9.snps.refseqgenes | sort -u | wc -l #175,818

awk -F'\t' -v OFS='\t' '{s=$2; sub(/\(.*/, "", s); print $0, s;}' mm9.snps.refseqgenes >mm9.snps.stripped.refseqgenes
cut -f8 mm9.snps.stripped.refseqgenes | sort -u | head
cut -f8 mm9.snps.stripped.refseqgenes | sort -u | wc -l #22,891
awk -F'\t' -v OFS='\t' '{print $1, $3, $4, $5, $6, $7, $8;}' mm9.snps.stripped.refseqgenes | sort -u >mm9.snps.refseqgenes
wc -l mm9.snps.refseqgenes #7,834,497

#filter on genes with at least two vars
cut -f7 mm9.snps.refseqgenes | sort | uniq -c >gene2nvars
wc -l gene2nvars #22,891
cat gene2nvars | awk '$1 >1 {print $0;}' >gene2nvars.minvars_2
wc -l gene2nvars.minvars_2 #22,181. Including introns

#
awk -F'\t' -v OFS='\t' '$1 !~ /intronic/' mm9.snps.refseqgenes >mm9.snps.refseq.nointrons
wc -l mm9.snps.refseq.nointrons #307,991
cut -f7 mm9.snps.refseq.nointrons | sort | uniq -c >gene2nvars
wc -l gene2nvars #21,685
cat gene2nvars | awk '$1 >1 {print $0;}' >gene2nvars.minvars_2
wc -l gene2nvars.minvars_2 #20,268. Excluding introns
cat mm9.snps.refseq.nointrons | cut -f1 | sort -u #OK

#get number of vars in genes with at least two vars
#see ngenes.R

#Split by chr
cut -f2,3 mm9.snps.refseqgenes | tr '\t' ' ' >mm9.refseq.pos
pos=${snpdir}/split/mm9.refseq.pos
cut -d' ' -f1 $pos | sort -u >real.chroms
cat real.chroms | grep -v 'random' >file.tmp
mv file.tmp real.chroms

screen -S snp2chrsplit
snpdir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/sanger_mouse_proj'
pos=${snpdir}/split/mm9.refseq.pos
cat real.chroms | while read chr
do
    echo $chr
    awk -v jchr=${chr} '$1 == jchr' $pos >${snpdir}/split/${chr}.pos
done


###
#Merged BAM files
###
#/Volumes/Data/cloud/btsync/work/rspd/projects/support/bjorn_reinius/ase.bl6cast.sh
#Only liver here: /mnt/kauffman/sandberglab/pipeline3.0/rnaseq/mmu/bjorn_reinius_X0-project/star_merged_cast1_mm9
#/mnt/crick/danielr/Xandclones_BR/bjorn_reinius_complete_set/star_cast1mm9
#bamdir='/mnt/crick/danielr/Xandclones_BR/bjorn_reinius_complete_set/star_mm9'

scriptdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/scripts/mousehybrid'
bamdir='/mnt/crick/danielr/Xandclones_BR/bjorn_reinius_complete_set/star_merged_cast1_mm9'
sample2id='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/data/mousehybrid/old.new.cell_names.txt'
rpkm='/mnt/crick/danielr/Xandclones_BR/bjorn_reinius_complete_set/rpkmforgenes/star_merged_cast1_mm9/refseq/rpkms_counts_rmnameoverlap.txt'

#copy rpkmforgenes file to localx
cd '/Volumes/Data/cloud/btsync/work/rspd/projects/scphaser/nogit/data/mousehybrid/rpkm'
scp -p 130.237.142.53:'/mnt/crick/danielr/Xandclones_BR/bjorn_reinius_complete_set/rpkmforgenes/star_merged_cast1_mm9/refseq/rpkms_counts_rmnameoverlap.txt' .

#Filter on 336 cells used in the scphaser paper
cd $scriptdir
\ls -1 $bamdir | xargs -I% echo % ${bamdir}/% | sort -u -k1,1 >bam.files
cat $sample2id | tr '\t' ' ' | sort -k1,1 >sample2id
join -1 1 -2 1 sample2id bam.files >sample2id2bamfile
wc -l sample2id2bamfile #336. OK

bamdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/data/mousehybrid/bam/star_merged'
mkdir -p $bamdir
awk -v d=${bamdir} '{print "ln -s", $3, d"/."}' sample2id2bamfile >ln.cmds.sh

###
#Mapstats
###
make_summary_starlog.pl ${bamdir} >mapstats.tab #no file: Log.final.out. Of course, since this is the merged bam files...
#instead, do line counting

find -L ${bamdir} -name '*.bam' >bam.files
cat bam.files | xargs -I% echo 'samtools view' % '| grep -v "ERCC-" | wc -l >'%'.nlines' >cmds.sh

