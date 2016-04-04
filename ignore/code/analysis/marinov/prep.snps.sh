

#Filters: snps (6.2M -> 5.4M), genes (5.4M -> 2.2M), dbSNP (2.2M -> 2.1M), het (2.1M -> 1.0M)

scriptdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/marinov_SRP018838/scripts'
snpdir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/vars'

#############
#Filter out indels
#############
vcf='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/vars/CEUTrio.HiSeq.WGS.b37.bestPractices.phased.hg19.vcf'
snps='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/vars/hg19.snps.vcf'
awk 'length($4) == 1 && length($5) == 1' $vcf >$snps
#6,240,465 -> 5,353,397


###############
#Filter on genes
##############
#Annotate with Annovar

annodir='/mnt/kauffman/edsgard/prg/annovar'
cd $annodir

#Download hg19 refseq
annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene humandb/

cd $snpdir

#Convert vcf to annovar format
#chromosome, start position, end position, the reference nucleotides and the observed nucleotides
#convert2annovar.pl -format vcf4 example/ex2.vcf > ex2.avinput
awk -F'\t' -v OFS='\t' '{print $1, $2, $2, $4, $5, $3;}' $snps >hg19.snps.avformat

#split input for parallel
mkdir split
cd split
split -n 70 ../hg19.snps.avformat hg19.snps.
cd $scriptdir
find ${snpdir}/split -name 'hg19.snps.*' >infiles.list

#Annovar
cat infiles.list | xargs -I% echo 'annotate_variation.pl -build hg19 --geneanno -dbtype refGene' % '/mnt/kauffman/edsgard/prg/annovar/humandb/' >cmds.sh
screen -S snp.annovar
cat cmds.sh | parallel -n 1 -P 70

#Extract SNPs with gene annot
cd ${snpdir}/split
cat *.variant_function >hg19.snps.variant_function
wc -l hg19.snps.variant_function #5,353,352

awk -F'\t' '{print $1;}' hg19.snps.variant_function | sort -u
#downstream, exonic, exonic;splicing, intergenic, intronic, ncRNA_exonic, ncRNA_intronic, ncRNA_splicing, splicing, upstream, upstream;downstream, UTR3, UTR5, UTR5;UTR3
#keep: exonic, exonic;splicing, intronic, ncRNA_exonic, ncRNA_intronic, ncRNA_splicing, splicing, UTR3, UTR5, UTR5;UTR3
#rm: intergenic, downstream, upstream, upstream;downstream: As far as I remember 1kb up/down is included by annovar by default.
awk -F'\t' -v OFS='\t' '$1 != "intergenic" && $1 != "downstream" && $1 != "upstream" && $1 != "upstream;downstream" {print $0;}' hg19.snps.variant_function >hg19.snps.refseqgenes
wc -l hg19.snps.refseqgenes #2,171,090
cut -f2 hg19.snps.refseqgenes | sort -u | wc -l #65,753
cut -f2 hg19.snps.refseqgenes | sort -u | head #TODO: Need to strip parenthesis info from gene-name

awk -F'\t' -v OFS='\t' '{s=$2; sub(/\(.*/, "", s); print $0, s;}' hg19.snps.refseqgenes >hg19.snps.stripped.refseqgenes
cut -f9 hg19.snps.stripped.refseqgenes | sort -u | head
cut -f9 hg19.snps.stripped.refseqgenes | sort -u | wc -l #23,603
awk -F'\t' -v OFS='\t' '{print $1, $9, $3, $4, $5, $6, $7, $8;}' hg19.snps.stripped.refseqgenes | sort -u >hg19.snps.refseqgenes
wc -l hg19.snps.refseqgenes #2,171,090
cut -f2 hg19.snps.refseqgenes | sort -u | wc -l #23,603


##############
#Filter on dbSNPs
##############
#Even if we remove a proportion of genes this might be a bit more "fair" value to report since in cases with no DNA-seq and variants derived from RNA-seq variant-calling one would most probably want to filter on dbSNPs.
awk '$8 != "."' hg19.snps.refseqgenes >hg19.dbsnp.refseq
wc -l hg19.dbsnp.refseq #2,071,501

#check against vcf to ensure that nothing was lost
awk -F'\t' -v OFS='\t' '$3 == "chr1"' hg19.dbsnp.refseq | sort -k4,4 - >chr1.hg19.dbsnp.refseq.possorted
awk -F'\t' -v OFS='\t' '$1 == "chr1" {print $1, $2;}' $vcf | sort -k2,2 >chr1.vcf.sorted
join -1 4 -2 2 chr1.hg19.dbsnp.refseq.possorted chr1.vcf.sorted >join.test
wc -l chr1.hg19.dbsnp.refseq.possorted #176,297
wc -l join.test #176,297. OK


###########
#Filter on heterozygous, since the CEUtrio input file covers all snps that are just different from the reference genome!!
##########
#vcf genotype field:
#GT : genotype, encoded as allele values separated by either of / or |. The allele values are 0 for the reference allele (what is in the REF field), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1, 1 | 0, or 1/2, etc. For haploid calls, e.g. on Y, male non- pseudoautosomal X, or mitochondrion, only one allele value should be given; a triploid call might look like 0/0/1. If a call cannot be made for a sample at a given locus, ‘.’ should be specified for each missing allele in the GT field (for example ‘./.’ for a diploid genotype and ‘.’ for haploid genotype). The meanings of the separators are as follows (see the PS field below for more details on incorporating phasing information into the genotypes): / : genotype unphased; | : genotype phased
#PL : the phred-scaled genotype likelihoods rounded to the closest integer (and otherwise defined precisely as the GL field)
#GL : genotype likelihoods comprised of comma separated floating point log10-scaled likelihoods for all possible genotypes given the set of alleles defined in the REF and ALT fields. In presence of the GT field the same ploidy is expected and the canonical order is used; without GT field, diploidy is assumed. If A is the allele in REF and B,C,... are the alleles as ordered in ALT, the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j. In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc. For example: GT:GL 0/1:-323.03,-99.29,-802.53 (Floats)
#GQ : conditional genotype quality, encoded as a phred quality −10log10 p(genotype call is wrong, conditioned on the site’s being variant) (Integer)
#FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
#FORMAT=<ID=TP,Number=1,Type=Integer,Description="Phred score of the genotype combination and phase given that the genotypes are correct">
#GT:AD:DP:GQ:PL:TP	0/1:96,58:156:95:95,0,451:94
#Filter:
#GT: 0/1 or 0|1
#AD: any >5
cd $snpdir
mindp=5
awk -v mdp=${mindp} -F'\t' -v OFS='\t' '{f=$10; split(f, a, ":"); split(a[2], ad, ","); if((a[1]=="0/1" || a[1]=="0|1" || a[1] == "1/0" || a[1] == "1|0") && ad[1] >= mdp && ad[2] >= mdp){print $0};}' $vcf >hg19.trio.hetfilt.vcf

wc -l $vcf #6,240,465
wc -l hg19.trio.hetfilt.vcf #2,932,458

#NOTA BENE!
#Requiring a depth on the ref allele, basically excludes het.variants where two alt alleles (two 0/1 records at the same position). Or, possibly such a variant is denoted 1/2. That depends on if they have split multi-allelic calls in the vcf file.
#Let's skip these special cases, hopefully not too many...


########
#Filter on PASS (but I think they used all vars in the vcf file to create the genomes and it also says "bestpractices" so I guess it's already been filtered?)
########
awk -F'\t' '$7 == "PASS" {print $0;}' /mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/vars/hg19.trio.hetfilt.vcf >${snpdir}/hg19.trio.het.pass.vcf
#2,327,432


###############
#Liftover REF SNP positions to paternal and maternal positions
###############
#liftOver oldFile map.chain newFile unMapped
#oldFile and newFile are in bed format by default
#The map.chain file has the old genome as the target and the new genome as the query
#query of maternal.chain: maternal (see header)
#
#CHAIN files should be 0-based, so, input positions should also be 0-based.
#http://genome.ucsc.edu/goldenPath/help/chain.html

hg19_bed='hg19.dbsnp.refseq.0based.bed'
diplosnpdir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/NA12878_diploid_genome_dec16_2012'
matchain=${diplosnpdir}/maternal.chain
patchain=${diplosnpdir}/paternal.chain
mat_bed='mat.dbsnp.refseq.bed'
pat_bed='pat.dbsnp.refseq.bed'

#annot -> bed
awk -F'\t' '{print $3, $4, $5 + 1;}' hg19.dbsnp.refseq | sed 's/^chr//' >hg19.dbsnp.refseq.1based.bed

#VCF is 1-based so remove 1 (CHAIN files are 0-based according to their def: http://genome.ucsc.edu/goldenPath/help/chain.html)
cat hg19.dbsnp.refseq.1based.bed | awk '{print $1, $2 - 1, $3 - 1;}' >$hg19_bed

#print cmds
echo liftOver $hg19_bed $matchain $mat_bed mat.hg19tomat.unmapped >mat.lift.sh
echo liftOver $hg19_bed $patchain $pat_bed pat.hg19tomat.unmapped >pat.lift.sh
cat mat.lift.sh pat.lift.sh >lift.sh

#run
screen -S lift
sh lift.sh
#status: fin

#change chromnames
cat mat.dbsnp.refseq.bed | sed 's/_maternal//' | awk '{print "chr"$1, $2;}' >mat.dbsnp.refseq.pos
cat pat.dbsnp.refseq.bed | sed 's/_paternal//' | awk '{print "chr"$1, $2;}' >pat.dbsnp.refseq.pos
cat $hg19_bed | awk '{print "chr"$1, $2;}' >hg19.dbsnp.refseq.pos

sort mat.hg19tomat.unmapped >m.unmapped
sort pat.hg19tomat.unmapped >p.unmapped
comm -3 m.unmapped p.unmapped
#diffs, so need to filter one at a time

cat m.unmapped | grep -v '^#' | awk '{print $1"."$2;}' >m.unmapped.pos
cat p.unmapped | grep -v '^#' | awk '{print $1"."$2;}' >p.unmapped.pos
cat hg19.dbsnp.refseq.0based.bed | awk '{print $1"."$2, $1, $2;}' >hg19.dbsnp.refseq.keyed.pos

#rm snps not mappable from ref-file
srcdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/git/ignore/code/analysis/marinov'
python3 ${srcdir}/setdiff.py m.unmapped.pos hg19.dbsnp.refseq.keyed.pos --list_keyfield 1 --file_keyfield 1 >ref.mpass.pos
python3 ${srcdir}/setdiff.py p.unmapped.pos hg19.dbsnp.refseq.keyed.pos --list_keyfield 1 --file_keyfield 1 >ref.ppass.pos
wc -l ref.mpass.pos #2071489
wc -l ref.ppass.pos #2071489

#paste ref.mpass.pos and mat.dbsnp.refseq.pos
paste ref.mpass.pos mat.dbsnp.refseq.pos | awk '{print $1, "chr"$2, $3, $5;}' >ref2mat.pos
paste ref.ppass.pos pat.dbsnp.refseq.pos | awk '{print $1, "chr"$2, $3, $5;}' >ref2pat.pos

#join 
cat ref2mat.pos | sort -k1,1 >file.tmp
mv file.tmp ref2mat.pos
cat ref2pat.pos | sort -k1,1 >file.tmp
mv file.tmp ref2pat.pos
join -1 1 -2 1 ref2mat.pos ref2pat.pos | awk -v OFS='\t' '{print "chr"$1, $2, $3, $4, $7;}' | sort -k1,1 >ref2mat2pat.pos
cat ref2mat2pat.pos | wc -l #2,071,501 -> 2,071,487. OK!

#filter empty dbsnp field (also, 1-based -> 0-based, and add key-field)
awk -F'\t' '$8 !="" {print $3"."$4-1, $3, $4-1, $6, $7, $8, $2, $1;}' hg19.dbsnp.refseq | sort -k1,1 >hg19.dbsnp.refseq.0based
wc -l hg19.dbsnp.refseq.0based #2,071,500

#Add gene and dbsnp annot
join -1 1 -2 1 hg19.dbsnp.refseq.0based ref2mat2pat.pos | cut -d' ' -f1-8,11,12 | sort -k1,1 >ref2mat2pat.dbsnp.refseq
wc -l ref2mat2pat.dbsnp.refseq #2,071,485


######
#Filter on het.pos
######

#join ref2mat2pat.dbsnp.refseq and hg19.trio.hetfilt.vcf. NOTA BENE: Also shifting vcf pos from 1-based to 0-based.
cat hg19.trio.het.pass.vcf | awk -F'\t' '{print $1"."$2 - 1, $1, $2 - 1;}' | sort -k1,1 >hg19.trio.het.pass.pos
join -1 1 -2 1 hg19.trio.het.pass.pos ref2mat2pat.dbsnp.refseq | cut -d' ' -f1,4-12 | tr ' ' '\t' >ref2mat2pat.het.dbsnp.refseq

#check
awk -F'\t' '{print NF;}' ref2mat2pat.het.dbsnp.refseq | sort -u #10. OK
wc -l ref2mat2pat.het.dbsnp.refseq #830,200


########
#Get haplotype for given positions
########
#IN: ref2mat2pat.het.dbsnp.refseq (positions), mat.fa, pat.fa
#OUT: positions (ref, mat, pat) and phased alleles (can also keep the ref, alt columns to allow double-checking that usually the same two alleles)
#also double-check against read counts for a few het.snps in the mat and pat ac.files
#
#NB: snp2hapalleles.py use 0-based coordinates in the posfile
screen -S snp2hap
srcdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/git/ignore/code/analysis/marinov'
snpdir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/vars'
posfile=${snpdir}/ref2mat2pat.het.dbsnp.refseq
g1='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/NA12878_diploid_genome_dec16_2012/maternal.fa'
g2='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/NA12878_diploid_genome_dec16_2012/paternal.fa'
python3 ${srcdir}/snp2hapalleles.py $posfile $g1 $g2 --chr_fieldpos 2 --g1_fieldpos 9 --g2_fieldpos 10 >${snpdir}/ref2mat2pat.het.dbsnp.refseq.hap
#status: fin
cat ref2mat2pat.het.dbsnp.refseq.hap | head -100 #OK. ref,alt and hap alleles agree:)

#NOW, convert from 0- to 1-based since samtools view is 1-based (http://www.htslib.org/doc/samtools-0.1.19.html). samtools mpileup is therefore also 1-based.
mv ref2mat2pat.het.dbsnp.refseq.hap ref2mat2pat.het.dbsnp.refseq.hap.0based
awk -F'\t' -v OFS='\t' '{print $2"."$3+1, $2, $3+1, $4, $5, $6, $7, $8, $9+1, $10+1, $11, $12;}' ref2mat2pat.het.dbsnp.refseq.hap.0based >ref2mat2pat.het.dbsnp.refseq.hap.1based



######
#Split by chr
######
cut -f2 ref2mat2pat.het.dbsnp.refseq.hap.1based | sort -u >real.chroms
cat real.chroms | while read chr
do
    echo $chr
    awk -v jchr=${chr} '$2 == jchr {print $2, $9}' ref2mat2pat.het.dbsnp.refseq.hap.1based >mat.${chr}.pos
    awk -v jchr=${chr} '$2 == jchr {print $2, $10}' ref2mat2pat.het.dbsnp.refseq.hap.1based >pat.${chr}.pos
done


########
#List of unphased het.snps
########
snpdir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/vars'
cd $snpdir
cat hg19.trio.het.pass.vcf | awk -F '\t' -v OFS='\t' '$10 ~ /\|/' >hg19.trio.het_phased.pass.vcf
wc -l hg19.trio.het_phased.pass.vcf #2,327,432 -> 1,899,089
awk '$3 != "." && $3 != "" {print $3;}' hg19.trio.het_phased.pass.vcf | sort >hg19.trio.het_phased.pass.rsid
wc -l *.rsid #1,850,225
cp -p hg19.trio.het_phased.pass.rsid /mnt/kauffman/edsgard/cloud/btsync/work/rspd/data/external/Marinov_GenRes_2014/.

#join with ref2mat2pat.het.dbsnp.refseq.hap.1based
sort -k6,6 ref2mat2pat.het.dbsnp.refseq.hap.1based >file.tmp
join -1 1 -2 6 hg19.trio.het_phased.pass.rsid file.tmp >ref2mat2pat.het.dbsnp.refseq.trio_phased.hap.1based
wc -l ref2mat2pat.het.dbsnp.refseq.trio_phased.hap.1based #830,200 -> 689,360


#######
#Filter on genes with at least two variants
########
cut -d' ' -f7 ref2mat2pat.het.dbsnp.refseq.trio_phased.hap.1based | sort | uniq -c >gene2nvars
wc -l gene2nvars #18,447
cat gene2nvars | awk '$1 >1 {print $0;}' >gene2nvars.minvars_2
wc -l gene2nvars.minvars_2 #15,597. Including introns


########
#Filter out intronic vars
########
#original: exonic, exonic;splicing, intronic, ncRNA_exonic, ncRNA_intronic, ncRNA_splicing, splicing, UTR3, UTR5, UTR5;UTR3
#rm: intronic, ncRNA_intronic

awk -v '$8 !~ /intronic/ {print $0;}' ref2mat2pat.het.dbsnp.refseq.trio_phased.hap.1based >ref2mat2pat.het.dbsnp.refseq.trio_phased.nointrons.hap.1based
wc -l ref2mat2pat.het.dbsnp.refseq.trio_phased.nointrons.hap.1based #689,360 -> 28,995
#double-check
grep -v 'intronic' ref2mat2pat.het.dbsnp.refseq.trio_phased.hap.1based | wc -l #28,995

cut -d' ' -f7 ref2mat2pat.het.dbsnp.refseq.trio_phased.nointrons.hap.1based | sort | uniq -c >gene2nvars
wc -l gene2nvars #10,661
cat gene2nvars | awk '$1 >1 {print $0;}' >gene2nvars.minvars_2
wc -l gene2nvars.minvars_2 #6,065. Excluding introns.


#ref2mat2pat.het.dbsnp.refseq.trio_phased.hap.1based


##########
#Tests
############

#snp2hapalleles.py minimal example
cat /mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/NA12878_diploid_genome_dec16_2012/maternal.fa | head -3 >file.tmp
cat file.tmp file.tmp >test.fa
cat $posfile | head -2 >pos.test
python3 ${srcdir}/snp2hapalleles.py pos.test g1.test.fa g2.test.fa --chr_fieldpos 1 --g1_fieldpos 8 --g2_fieldpos 9 | head
#OK, the scripts work. So, probably either the position mapping (chaining/liftover) is wrong or there is a 0- vs 1-based issue somehere earlier.
#Are the positions in the original vcf-file 1-based? Yes, should be!
#Is the chain file wrong (wrong direction?)
#Were the provided haplo-fasta files generated with the same chain files as those provided?
#Were the positions in the vcf file using the same hg19 ref-genome as was used to create the haplo-fasta files?
