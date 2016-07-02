

#Study accession: EGAS00001001009
#dataset accession: EGAD00001001083, EGAD00001001084

#Variant calling workflow: http://www.htslib.org/workflow

wget https://www.ebi.ac.uk/ega/about/your_EGA_account/download_streaming_client
java -jar /mnt/kauffman/edsgard/prg/egademoclient/EgaDemoClient.jar
login daniel.edsgard@ki.se
datasets

request dataset EGAD00001001084 abc rnaseq
download rnaseq
#status: fin

request dataset EGAD00001001083 def dnaseq
download dnaseq
#status: fin

screen -S ega.meta
requestDownloadmetadata EGAD00001001084
requestDownloadmetadata EGAD00001001083

#Decrypt from command line
cd '/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009/rawdata'
ls -1 *.fastq.gz.cip >encrypted.files
cat encrypted.files | xargs -I% echo 'java -jar /mnt/kauffman/edsgard/prg/egademoclient/EgaDemoClient.jar -p daniel.edsgard@ki.se password -dc' % '-dck abc' >cmds.sh
screen -S ega_rna
sh cmds.sh
#status: fin

cd '/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009/dnaseq'
ls -1 *.fastq.gz.cip >encrypted.files
cat encrypted.files | xargs -I% echo 'java -jar /mnt/kauffman/edsgard/prg/egademoclient/EgaDemoClient.jar -p daniel.edsgard@ki.se password -dc' % '-dck def' >cmds.sh
screen -S ega_dna
sh cmds.sh



###################
#Align RNA-seq data
###################

###
#Paired end mapping not working in pipeline. Make star calls without pipeline:
###
projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009'
scriptdir=${projdir}'/scripts'
stardir='/mnt/kauffman/sandberglab/star_index/hg19_ERCC_eGFP'
cd $scriptdir

#suffix: _1.fastq.gz
find ${projdir}/rawdata -name '*_1.fastq.gz' | sed 's/_1.fastq.gz//' >samples.prefix
cat samples.prefix | sed 's/rawdata/rnaseq\/star_hg19/' >out.dirs
paste samples.prefix out.dirs >star.in
p=10

cat out.dirs | xargs -I% mkdir -p %
cat star.in | awk -v p=$p -v r=$stardir '{print "STAR --runThreadN", p, "--genomeDir", r, "--readFilesIn", $1"_1.fastq.gz", $1"_2.fastq.gz", "--readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix", $2"/", ">"$2"/star.out 2>"$2"/star.err";}' >star.chunk1.cmds.sh

#suffix: R1_001.fastq.gz
find ${projdir}/rawdata -name '*_R1_001.fastq.gz' | sed 's/_R1_001.fastq.gz//' >samples.prefix
cat samples.prefix | sed 's/rawdata/rnaseq\/star_hg19/' >out.dirs
paste samples.prefix out.dirs >star.in

cat out.dirs | xargs -I% mkdir -p %
cat star.in | awk -v p=$p -v r=$stardir '{print "STAR --runThreadN", p, "--genomeDir", r, "--readFilesIn", $1"_R1_001.fastq.gz", $1"_R2_001.fastq.gz", "--readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix", $2"/", ">"$2"/star.out 2>"$2"/star.err";}' >star.chunk2.cmds.sh

cat star.chunk1.cmds.sh star.chunk2.cmds.sh >star.cmds.sh
wc -l star.cmds.sh #209

#check 
cat star.cmds.sh | cut -d' ' -f7 | xargs -I% ls -lrth % | awk '{print $5;}'
cat star.cmds.sh | cut -d' ' -f8 | xargs -I% ls -lrth % | awk '{print $5;}'
#OK

screen -S star
cat star.cmds.sh | parallel -n 1 -P 5
#status: fin

#check
cat sam.files | xargs -I% ls -lrth % | awk '{print $5;}' 


###
#filter uniquely mapped
###
bamdir=${projdir}/rnaseq/star_hg19
find ${bamdir} -name 'Aligned.out.sam' >sam.files
cat sam.files | xargs -I% echo awk -F"'\t'" "'\$1 ~ /^@/ || \$5 == 255 {print \$0;}'" % '| samtools view -b - | samtools sort -O bam -T '%' - >'%'.unique.bam' >filter.multi.cmds.sh

screen -S filter.multi
cat filter.multi.cmds.sh | parallel -n 1 -P 100
#status: fin

#check
find $bamdir -name '*.bam' | grep -v 'TN2A' | grep -v -i 'library' | grep -v -i 'pool' >bamfiles.list
cat bamfiles.list | wc -l #163. OK
cat bamfiles.list | xargs -I% ls -lrth % | awk '{print $5;}' | sort -u #OK


########
#Check meta for number of single cells per sample
##########
cd 'meta/EGAD00001001084/delimited_maps'
cat Run_Sample_meta_info.map | grep 'UCF1014' | wc -l #169
cat Run_Sample_meta_info.map | grep 'UCF1014' | grep 'CAPT' | wc -l #163
cat Run_Sample_meta_info.map | grep 'T2N' | wc -l #40

cd 'meta/EGAD00001001083'
cat Sample_File.map
cat Run_Sample_meta_info.map



##############
#Align DNA-seq data
##############

#create bwa index
cat /mnt/kauffman/sandberglab/genomes/chrom/hg19_norandom/*.fa >/mnt/kauffman/edsgard/nobackup/data/annot/genomes/hg19.concat.fa
mkdir -p /mnt/kauffman/edsgard/nobackup/data/annot/bwa/hg19
cd /mnt/kauffman/edsgard/nobackup/data/annot/bwa/hg19
bwa index /mnt/kauffman/edsgard/nobackup/data/annot/genomes/hg19.concat.fa
#status: fin

projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009'
scriptdir=${projdir}'/scripts'
bwaindex='/mnt/kauffman/edsgard/nobackup/data/annot/bwa/hg19/hg19.concat.fa'
cd $scriptdir

find ${projdir}/dnaseq -name '*UCF*1014*R1*' | sed 's/_R1_001.fastq.gz//' >samples.prefix
cat samples.prefix | sed 's/dnaseq/dnaseq\/bwa_hg19/' >out.dirs
paste samples.prefix out.dirs >bwa.in
cat out.dirs | xargs -I% mkdir -p %
p=40
cat bwa.in | awk -v p=${p} -v r=${bwaindex} '{print "bwa mem -t", p, r, $1"_R1_001.fastq.gz", $1"_R2_001.fastq.gz", "| samtools view -b - >", $2".bam"}' >bwa.cmds.sh

screen -S bwa
cat bwa.cmds.sh | parallel -n 1 -P 2
#status: fin


#Merge bamfiles
bamdir=${projdir}/dnaseq/bwa_hg19
find $bamdir -name '*UCF*' >bam.files
echo 'samtools merge' -b bam.files ${bamdir}/UCF1014.bam >merge.cmds.sh

screen -S merge
sh merge.cmds.sh
#status: fin

#Sort
bamdir=${projdir}/dnaseq/bwa_hg19
bam=${bamdir}/UCF1014.bam
echo 'samtools sort -O bam -T /tmp/ucf1014.bamsort' $bam '>'${bam}.sorted >bamsort.cmds.sh
screen -S bam.sort
sh bamsort.cmds.sh

#Index
screen -S index
samtools index /mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009/dnaseq/bwa_hg19/UCF1014.bam.sorted


###############################
#DNA variant calling
################################

#TBD: -q255 -d2*mean cov | vcfutils.pl varFilter
projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009'
scriptdir=${projdir}'/scripts'
bamdir=${projdir}/dnaseq/bwa_hg19
vcfdir=${projdir}/dnaseq/vcf
ref='/mnt/kauffman/edsgard/nobackup/data/annot/bwa/hg19/hg19.concat.fa'
bam=${bamdir}/UCF1014.bam.sorted
vcf=${vcfdir}/UCF1014.bcf

cd $scriptdir
echo 'samtools mpileup -Q20 -ugf' $ref $bam '| bcftools call -vmO b -o' $vcf >dna.pileup.cmds.sh

mkdir $vcfdir
screen -S dna.pileup
sh dna.pileup.cmds.sh
bcftools view $vcf | wc -l #4,637,369

#Filter out indels
bcftools view $vcf | awk 'length($4) == 1 && length($5) == 1' >${vcf}.vcf.snvs
wc -l ${vcf}.vcf.snvs #3,991,604

#Filter on heterozygous
awk -F'\t' '{f=$10; split(f, a, ":"); if((a[1]=="0/1" || a[1]=="0|1" || a[1] == "1/0" || a[1] == "1|0")){print $0};}' ${vcf}.vcf.snvs >${vcf}.vcf.snv.het
wc -l ${vcf}.vcf.snv.het #2,528,388

#Filter on DP
mindp=5
awk -v mindp=${mindp} -F'\t' '{f=$8; split(f, a, ";"); split(a[13], ad, "="); split(ad[2], dp, ","); if(((dp[1] + dp[2])>=mindp && (dp[3] + dp[4])>=mindp)){print $0};}' ${vcf}.vcf.snv.het >${vcf}.vcf.snv.het.dp
wc -l ${vcf}.vcf.snv.het.dp
#2,278,823

#Filter on GQ
#no GQ output

#Intersect with dbSNP snps
dbsnpdir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/dbsnp/20141222/annovar'
dbsnp=${dbsnpdir}/hg19.sao.vld.refgene.sorted.tab #see embryodevel/code/hsa/general/dbsnp.sh
wc -l $dbsnp #19,162,728
awk -F'\t' -v OFS='\t' '{print $1"."$2"."$3"."$4, $0;}' ${dbsnp} | tr '\t' '%' | sort -t'%' -k1,1 >${dbsnp}.chrposrefalt
awk -F'\t' -v OFS='\t' '{print $1"."$2"."$4"."$5, $0;}' ${vcf}.vcf.snv.het.dp | tr '\t' '%' | sort -t'%' -k1,1 >${vcf}.vcf.snv.het.dp.chrposrefalt
join -t'%' -1 1 -2 1 ${vcf}.vcf.snv.het.dp.chrposrefalt ${dbsnp}.chrposrefalt | tr '%' '\t' >${vcf}.vcf.snv.het.dp.gene.dbsnp
cut -f1-11,16-18 ${vcf}.vcf.snv.het.dp.gene.dbsnp >file.tmp
mv file.tmp ${vcf}.vcf.snv.het.dp.gene.dbsnp
wc -l ${vcf}.vcf.snv.het.dp.gene.dbsnp #787,883

#double-check that both 1-based. I know that VCF is 1-based, but the dbSNP file may be 0-based
#final input list to samtools mpileup should be 1-based
awk '{print $2"."$3-1;}' ${vcf}.vcf.snv.het.dp.chrpos | sort >test.vcf
join -1 1 -2 1 test.vcf ${pos}.chrpos >test.joined
wc -l test.joined #29,243

#Print only chr,pos
awk -F'\t' '{print $2, $3;}' ${vcf}.vcf.snv.het.dp.gene.dbsnp >${vcf}.vcf.snv.het.dp.gene.dbsnp.pos

#Split by chr
pos=${vcf}.vcf.snv.het.dp.gene.dbsnp.pos
cut -d' ' -f1 $pos | sort -u >real.chroms
cat real.chroms | while read chr
do
    echo $chr
    awk -v jchr=${chr} '$1 == jchr' $pos >${vcfdir}/${chr}.pos
done

#Get number of genes with at least two variants
awk -F'\t' -v OFS='\t' '{s=$2; sub(/\(.*/, "", s); print $0, s;}' hg19.snps.refseqgenes >hg19.snps.stripped.refseqgenes

cut -f13 ${vcf}.vcf.snv.het.dp.gene.dbsnp | sort | uniq -c >gene2nvars
wc -l gene2nvars #18,447
cat gene2nvars | awk '$1 >1 {print $0;}' >gene2nvars.minvars_2
wc -l gene2nvars.minvars_2 #15,597. Including introns


###########################
#RNA-seq allele-counting
###########################
#SEE snpase.sh


##########
#Mapstats
##########

projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009'
scriptdir=${projdir}'/scripts'
cd $scriptdir

stardir=${projdir}/rnaseq/star_hg19
resdir=${projdir}/rnaseq/meta

mkdir $resdir

cd $resdir
make_summary_starlog.pl ${stardir} >mapstats.tab
wc -l mapstats.tab #210
head mapstats.tab
##nreads are all reads stated in Log.final.out, so, including multi-mapping reads.

cd '/Volumes/Data/cloud/btsync/work/rspd/projects/scphaser/nogit/data/borel'
scp -p 130.237.142.53:/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009/rnaseq/meta/mapstats.tab .
##See stats.R for seq-depth summary of uniquely mapping reads.


#############
#Rpkmforgenes
#############

###
#rnaseq.py pipeline
###
exp='borel_EGAS00001001009'
timestamp='20160627'

execdir='/mnt/kauffman/sandberglab/pipeline3.0'
conf=${HOME}'/cloud/btsync/work/rspd/data/annot/rpkmforgenes/rnaseq.kauffman.pipeline3.0.20150228.conf'
logdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/'${exp}'/logs'
projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/'${exp}
scriptdir=${projdir}'/scripts'
p='80'
ref='hg19'

#write pipecmd
cd $scriptdir
echo cd $execdir >rpkmforgenes.sh
echo python3 rnaseq.py -R -p $p -a $ref -e $exp -c $conf '>'${logdir}/rnaseqpipe.${timestamp}.out '2>' ${logdir}/rnaseqpipe.${timestamp}.err >>rpkmforgenes.sh

#filenames not sample-names so write commands explicitly
find -L ${projdir}/star_hg19 -name '*.bam' >bam.files
cat bam.files | tr '\n' ' ' >bam.sp
echo 'python /mnt/kauffman/sandberglab/pipeline3.0/scripts/rpkmforgenes.py -p' $p '-ulen -readcount -fulltranscript -mRNAnorm -rmnameoverlap -bothendceil -a /mnt/kauffman/sandberglab/pipeline3.0/additionaldata/rpkmforgenes/hg19refGene_eGFP.txt -bamu -i TMP -o /mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009/rpkmforgenes_star_hg19/refseq/rpkm_refseq.txt 1>'${logdir}/rpkmforgenes.out 2>${logdir}/rpkmforgenes.err >rpkmforgenes.sh

#mkdir
mkdir /mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009/rpkmforgenes_star_hg19/refseq

#run
screen -S ${exp}.rpkmforgenes
sh rpkmforgenes.sh
#status: fin

#copy
cd '/Volumes/Data/cloud/btsync/work/rspd/projects/scphaser/nogit/data/borel'
scp -p 130.237.142.53:/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009/rpkmforgenes_star_hg19/refseq/rpkm_refseq.txt .
#aligned reads to genes and expressed gene stats: SEE stats.R
