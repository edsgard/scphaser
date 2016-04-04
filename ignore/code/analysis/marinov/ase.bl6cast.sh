


#################################
#Keep read with best mapping to the two different genomes
#################################

#Params
exp='marinov_SRP018838'
species='hsa'

projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/'${species}'/'${exp}
logdir=${projdir}'/logs/ase'
scriptdir=${projdir}'/scripts'
srcdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/code'

ref1bamdir=${projdir}/'star_hg19_na12878/maternal'
ref2bamdir=${projdir}/'star_hg19_na12878/paternal'
outbamdir1=${projdir}/'star_hg19_na12878/besthit/maternal'
outbamdir2=${projdir}/'star_hg19_na12878/besthit/paternal'
g1name='maternal'
g2name='paternal'
p=28

#Print cmds
mkdir -p $scriptdir
mkdir -p $logdir
cd $scriptdir
echo python ${srcdir}/filter.bestmatch.twogenomes.py $ref1bamdir $ref2bamdir $outbamdir1 $outbamdir2 --g1name $g1name --g2name $g2name --procs $p '>'${logdir}/besthit.out '2>' ${logdir}/besthit.err >filter.nm.cmds.sh

#Execute
screen -S ase.mergebams
sh filter.nm.cmds.sh
#status: fin

#fix filename (script now fixed such that using .sam as output-file suffix)
find ${outbamdir1} -name '*.bam' | awk '{f=$1; sub(/.bam/, ".sam", f); print "mv", $1, f;}' >cmds.sh
sh cmds.sh
find ${outbamdir2} -name '*.bam' | awk '{f=$1; sub(/.bam/, ".sam", f); print "mv", $1, f;}' >cmds.sh
sh cmds.sh

#check output
find ${outbamdir1} -name '*.sam' >mat.sam.files
find ${outbamdir2} -name '*.sam' >pat.sam.files

#index genome
fastadir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/NA12878_diploid_genome_dec16_2012'
samtools faidx ${fastadir}/maternal.fa
samtools faidx ${fastadir}/paternal.fa

#sam2bam
mat_fai=${fastadir}/'maternal.fa.fai'
pat_fai=${fastadir}/'paternal.fa.fai'
cat mat.sam.files | xargs -I% echo 'samtools view -bt' $mat_fai % '>'%'.bam' | sed 's/.sam.bam$/.bam/' >mat.sam2bam.sh
cat pat.sam.files | xargs -I% echo 'samtools view -bt' $pat_fai % '>'%'.bam' | sed 's/.sam.bam$/.bam/' >pat.sam2bam.sh
cat mat.sam2bam.sh pat.sam2bam.sh >sam2bam.sh
screen -S sam2bam
cat sam2bam.sh | parallel -n 1 -P 56
#status: fin

#sort
find ${outbamdir1} -name '*.bam' >mat.bam.files
find ${outbamdir2} -name '*.bam' >pat.bam.files
cat mat.bam.files pat.bam.files >mat.pat.bam.files
cat mat.pat.bam.files | xargs -I% echo 'samtools sort -O bam -T' % % '>'%.sorted | sed 's/.bam.sorted$/.sorted.bam/' >sort.sh
screen -S sort
cat sort.sh | parallel -n 1 -P 56
#status: fin

#index
find ${outbamdir1} -name '*.sorted.bam' >mat.bam.files
find ${outbamdir2} -name '*.sorted.bam' >pat.bam.files
cat mat.bam.files pat.bam.files >mat.pat.bam.files
cat mat.pat.bam.files | xargs -I% echo 'samtools index' % >index.sh
screen -S index
cat index.sh | parallel -n 1 -P 56

#remove originals
cat mat.sam.files pat.sam.files >sam.files
cat sam.files | xargs -I% rm %

find ${outbamdir1} -name '*unique.bam' >mat.bam.files
find ${outbamdir2} -name '*unique.bam' >pat.bam.files
cat mat.bam.files pat.bam.files >mat.pat.bam.files
cat mat.pat.bam.files | xargs -I% rm %

#check
find ${outbamdir1} -name '*.bam'
lt /mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/marinov_SRP018838/star_hg19_na12878/besthit/paternal/SRR764808


###############
#TESTs
###############
#Test best-hit filtering
samtools view -h /mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/marinov_SRP018838/star_hg19_na12878/maternal/SRR764782/Aligned.out.unique.bam | head -1000000 | samtools view -b - >file1.bam
samtools view -h /mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/marinov_SRP018838/star_hg19_na12878/paternal/SRR764782/Aligned.out.unique.bam | head -1000000 | samtools view -b - >file2.bam

file1=${PWD}/file1.bam
file2=${PWD}/file2.bam
sample='SRR764782'
echo python ${srcdir}/filter.minimal.py $file1 $file2 $outbamdir1 $outbamdir2 $sample --g1name $g1name --g2name $g2name '>'besthit.out '2>'besthit.err >cmds.sh


