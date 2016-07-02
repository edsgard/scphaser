


#Bcftools doc: http://samtools.github.io/bcftools/
#Samtools doc: http://www.htslib.org; http://www.htslib.org/doc/samtools-1.1.html

##DEPENDENCIES: prep.snps.sh, pileup2ac.py, ac2knownhap.R


###############
#Mpileup on dbSNPs
##############

exp='borel_EGAS00001001009'
species='hsa'
projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/'${species}'/'${exp}

##IN/OUT
scriptdir=${projdir}'/scripts'

##IN
bamdir=${projdir}/'rnaseq/star_hg19'
snpdir=${projdir}/'dnaseq/vcf'

##OUT
bamfilelist=${scriptdir}/'bamfiles.list'
vcfbasedir=${projdir}'/rnaseq/ase'


###
#List bamfiles
###

mkdir $scriptdir
cd $scriptdir

#bamfiles
find $bamdir -name '*.bam' | grep -v 'TN2A' | grep -v -i 'library' | grep -v -i 'pool' >$bamfilelist
cat bamfiles.list | wc -l #163. OK

#check
cat $bamfilelist | xargs -I% dirname % | xargs -I% basename % >bam.base.files
paste $bamfilelist bam.base.files >bam.full.base



######
#Split bam files on chrom (resolves the mem leakage when total size of bam-files are too big)
######

###
#Each chr
###
bamsubdir=${vcfbasedir}'/bam/chrsplit'
mkdir -p $bamsubdir
cat bam.full.base | awk -F'\t' -v outdir=${bamsubdir} '{print "bamtools split -in", $1, "-refPrefix \"\" -reference -stub "outdir"/"$2;}' >cmds.sh

screen -S bamsplit
cat cmds.sh | parallel -n 1 -P 80
#status: fin

#one dir of bamfiles per chr
cd $scriptdir
\ls -1 $bamsubdir | sed 's/.bam//' | sed 's/_EGAR.*\.//' | sort -u >chroms
cat chroms | xargs -I% mkdir ${bamsubdir}/%
cat chroms | xargs -I% echo mv ${bamsubdir}'/*'%'.bam' ${bamsubdir}/%/. >cmds.sh
sh cmds.sh

#check
lt ${bamsubdir}/chr1 | wc -l
lt ${bamsubdir}/chr2 | wc -l

cp -p chroms sel.chroms

#files with reads in all chroms
#number of files may differ, probably since some empty
#only take subset of files with reads in all chroms
cat sel.chroms | xargs -I% echo 'ls -1' ${bamsubdir}/% '>>chr.files' >cmds.sh
sh cmds.sh
cat chr.files | sed 's/.chr.*.bam//' | sort | uniq -c | awk '$1 == "23" {print $2;}' >pass.chr.bams
wc -l pass.chr.bams #28

#for each chrom create a file with bam files
cat sel.chroms | while read chr
do
cat pass.chr.bams | awk -v chr=$chr -v bdir=${bamsubdir} '{print bdir"/"chr"/"$1"."chr".bam"}' >${chr}.bam.files
done

#for each chrom create a file with bam files
cat sel.chroms | while read chr
do
\ls -1 ${bamsubdir}/${chr} | xargs -I% echo ${bamsubdir}/${chr}/% >${chr}.bam.files
done


#####
#Down-sample the number of variants
#####
N='250000'
N='500000'
N='750000'
projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/borel_EGAS00001001009'
scriptdir=${scriptdir}/'N'${N}
vcfdir=${projdir}/dnaseq/vcf
pos=${vcfdir}/UCF1014.bcf.vcf.snv.het.dp.gene.dbsnp.pos
downsampdir=${vcfdir}'_N'${N}
snpfile=${downsampdir}/UCF1014.bcf.vcf.snv.het.dp.gene.dbsnp.pos

mkdir $scriptdir
mkdir $downsampdir
cd $scriptdir
shuf -n $N ${pos} >${snpfile}


###
#Split snp-file by chr
###
cut -d' ' -f1 $snpfile | sort -u >real.chroms
cat real.chroms | while read chr
do
    echo $chr
    awk -v jchr=${chr} '$1 == jchr' $snpfile >${downsampdir}/${chr}.pos
done


#Symlink bam-file lists
ln -s ../chr*.bam.files .

######
#mpileup per sample, per chr
######

#min baseQual: 30. Overly strict?
#https://github.com/samtools/samtools/wiki/FAQ:
#Q. In the pileup output, the base quality is occasionally smaller than the read quality in BAM. Bug?
#A. By default, samtools computes Base Alignment Quality (BAQ) for each base. The quality column in pileup gives the minimum between BAQ and the original base quality, which may be lower than the base quality if the base is likely to be affected by indels nearby. It is a feature.
cat real.chroms | while read chr
do    
    cat ${chr}.bam.files | xargs -I% echo 'samtools mpileup -O -Q 30 -l' ${downsampdir}/${chr}.pos % '1>'%.downsample_${N}.pileup '2>'%.downsample_${N}.pileup.err >${chr}.pileup.cmds.sh
done
#    cat ${chr}.bam.files | xargs -I% echo 'samtools mpileup -Q20 -q255 -d1000 -t DP,DPR,INFO/DPR,SP,DV,DP4 -uvf /mnt/crick/sandberglab/genomes/hg19_norandom.fa -l' ${snpdir}/${hap}.${chr}.pos % '| vcfutils.pl varFilter | vcf-sort | bcftools view -O b 1>' %.bcf '2>'%.pileup.err >${chr}.pileup.cmds.sh
cat chr*.pileup.cmds.sh >pileup.downsample_${N}.cmds.sh 
wc -l pileup.downsample_${N}.cmds.sh #3912 (OK: 24 * 163)

screen -S pileup.N_${N}
N='250000'
cat pileup.downsample_${N}.cmds.sh | parallel -n 1 -P 50
#status: fin


###
#pileup2allelecounts
###
exp='borel_EGAS00001001009'
species='hsa'
projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/'${species}'/'${exp}
vcfbasedir=${projdir}'/rnaseq/ase'
bamsubdir=${vcfbasedir}'/bam/chrsplit'
srcdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/code'
scriptdir=${projdir}'/scripts'

N=750000

scriptdir=${scriptdir}/'N'${N}
namestr='*.downsample_'${N}'.pileup'

cd $scriptdir
find $bamsubdir -name $namestr >pileup.downsample_${N}.files
cat pileup.downsample_${N}.files | xargs -I% echo python3 ${srcdir}/pileup2ac.py % %.downsample_${N}.ac >pileup2ac.downsample_${N}.cmds.sh
echo $N
cat pileup2ac.downsample_${N}.cmds.sh | wc -l #3912
cat pileup2ac.downsample_${N}.cmds.sh | head

screen -S pileup2ac.downsample_${N}
N=750000
cat pileup2ac.downsample_${N}.cmds.sh | parallel -n 1 -P 30
#status: fin (250k, 500k, 750k)

N=500000
namestr='*.downsample_'${N}'.ac'
find $bamsubdir -name $namestr >ac.downsample_${N}.files
wc -l ac.downsample_${N}.files
head ac.downsample_${N}.files
cp -p ac.downsample_${N}.files /mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/scripts/borel/.


###
#ac2hapac
###
#get counts for haplotype-allele. Also, create a matrix for the specific haplotype (rows: snps, columns: samples).
#SEE: ac2mat.R
#status: 
