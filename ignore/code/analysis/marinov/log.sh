


#################################
#Download diploid genome
#################################
http://sv.gersteinlab.org/NA12878_diploid/
##Checking the README, Rozowsky helped creating this so most prob based on his allele-seq paper

cd '/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878'
wget 'http://sv.gersteinlab.org/NA12878_diploid/NA12878_diploid_genome_2012_dec16.zip'
unzip *
wget 'http://sv.gersteinlab.org/NA12878_diploid/CEUTrio.HiSeq.WGS.b37.bestPractices.phased.hg19.vcf.gz'
gunzip *.gz


##################################
#Download single-cell data
##################################
#GSE44618
#FPKMs
cd '/mnt/kauffman/edsgard/rsync/work/rspd/data/external/marinov'
mkdir fpkm
cd fpkm
wget 'http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE44618&format=file'
tar -xvf GSE44618_RAW.tar

#sra files (Fastq)
wget -r 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP018/SRP018838'
#status: fin (total time: 15h)

#meta
mkdir meta
cd meta
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44618/soft/GSE44618_family.soft.gz'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44618/miniml/GSE44618_family.xml.tgz'
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44618/matrix/GSE44618_series_matrix.txt.gz'


###################
#sra2sam
###################
scriptdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/scripts/marinov'
cd $scriptdir
datadir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/marinov/sra/SRP018838'
find $datadir -name '*.sra' >sra.files

#subset on single-cell files
#see samplefilter.R
cp -p '/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/data/marinov/singlecell.sra.files' .

#print cmds
cat singlecell.sra.files | xargs -I% echo 'sam-dump -r' % '| samtools view -b - >'%.bam >samdump.cmds.sh

#run
screen -S samdump
cat samdump.cmds.sh | parallel -n 1 -P 34
#status: fin

#mv bamfiles
find sra/SRP018838 -name '*.bam' | xargs -I% mv % bam/.

#Test
#sam-dump -r /mnt/kauffman/edsgard/rsync/work/rspd/data/external/marinov/ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP018/SRP018838/SRR1066622/SRR1066622.sra | head | samtools view -b - >test.bam


###############
#sra2fastq
###############
scriptdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/scripts/marinov'
cd $scriptdir
fastqdir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/marinov/fastq'
cat singlecell.sra.files | xargs -I% echo 'fastq-dump --split-3 --gzip -O' $fastqdir % >fastqdump.cmds.sh

#run
screen -S fastqdump
cat fastqdump.cmds.sh | parallel -n 1 -P 34
#status: fin

#Test
#fastq-dump --split-3 /mnt/kauffman/edsgard/rsync/work/rspd/data/external/marinov/ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP018/SRP018838/SRR1066622/SRR1066622.sra


################
#Build STAR index for each haplo-genome
################
scriptdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/scripts/marinov'
cd $scriptdir

#catenate chroms
fastadir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/NA12878_diploid_genome_dec16_2012'
ls -1 ${fastadir}/*_paternal.fa >paternal.fa.list
ls -1 ${fastadir}/*_maternal.fa >maternal.fa.list
cat paternal.fa.list | xargs -I % cat % >${fastadir}/paternal.fa
cat maternal.fa.list | xargs -I % cat % >${fastadir}/maternal.fa

#sub chrom fasta-headers
cat ${fastadir}/paternal.fa | sed 's/>/>chr/' | sed 's/_paternal//' >p.fa
mv p.fa ${fastadir}/paternal.fa
cat ${fastadir}/maternal.fa | sed 's/>/>chr/' | sed 's/_maternal//' >m.fa
mv m.fa ${fastadir}/maternal.fa

#build star index
#erccfasta='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/data/annot/ercc/ERCC92.fa'
#arabidopsis was used as spike-in, lets skip

screen -S genstarindex
fastadir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/NA12878_diploid_genome_dec16_2012'
stardir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/starindex/paternal'
fasta=${fastadir}/paternal.fa
STAR --runMode genomeGenerate --genomeDir $stardir --genomeFastaFiles $fasta --runThreadN 23 >${stardir}/star.index.out &

stardir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/starindex/maternal'
fasta=${fastadir}/maternal.fa
STAR --runMode genomeGenerate --genomeDir $stardir --genomeFastaFiles $fasta --runThreadN 23 >${stardir}/star.index.out &
#status: fin


###################
#Align
###################

###
#Paired end mapping not working in pipeline. Make star calls without pipeline:
###
projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/marinov_SRP018838'
stardir='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/gerstein_NA12878/starindex'
cd $scriptdir
find ${projdir}/rawdata -name '*_1.fastq.gz' | sed 's/_1.*//' >paired.srr
cat paired.srr | sed 's/rawdata/star_hg19_na12878\/maternal/' | xargs -I% dirname % >out.dirs
paste paired.srr out.dirs >star.in
p=1

ref=${stardir}/'maternal'
cat out.dirs | xargs -I% mkdir -p %
cat star.in | awk -v p=$p -v r=$ref '{print "STAR --runThreadN", p, "--genomeDir", r, "--readFilesIn", $1"_1.fastq.gz", $1"_2.fastq.gz", "--readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix", $2"/", ">"$2"/star.out 2>"$2"/star.err";}' >star.mat.paired.cmds.sh
screen -S star.maternal.paired
cat star.mat.paired.cmds.sh | parallel -n 1 -P 9
#status: fin

ref=${stardir}/'paternal'
cat out.dirs | sed 's/maternal/paternal/' >pat.out.dirs
cat star.in | sed 's/maternal/paternal/' >pat.star.in
cat pat.out.dirs | xargs -I% mkdir -p %
cat pat.star.in | awk -v p=$p -v r=$ref '{print "STAR --runThreadN", p, "--genomeDir", r, "--readFilesIn", $1"_1.fastq.gz", $1"_2.fastq.gz", "--readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix", $2"/", ">"$2"/star.out 2>"$2"/star.err";}' >star.pat.paired.cmds.sh
screen -S star.paternal.paired
cat star.pat.paired.cmds.sh | parallel -n 1 -P 9
#status: fin

###
#Single end files
###
find ${projdir}/rawdata -name '*.fastq.gz' | grep -v '_1.fastq.gz' | grep -v '_2.fastq.gz' >single.srr
cat single.srr | sed 's/rawdata/star_hg19_na12878\/maternal/' | xargs -I% dirname % >mat.single.out.dirs
cat single.srr | sed 's/rawdata/star_hg19_na12878\/paternal/' | xargs -I% dirname % >pat.single.out.dirs
p=23

ref=${stardir}/'maternal'
paste single.srr mat.single.out.dirs >mat.single.star.in
cat mat.single.out.dirs | xargs -I% mkdir -p %
cat mat.single.star.in | awk -v p=$p -v r=$ref '{print "STAR --runThreadN", p, "--genomeDir", r, "--readFilesIn", $1, "--readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix", $2"/", ">"$2"/star.out 2>"$2"/star.err";}' >star.mat.single.cmds.sh

ref=${stardir}/'paternal'
paste single.srr pat.single.out.dirs >pat.single.star.in
cat pat.single.out.dirs | xargs -I% mkdir -p %
cat pat.single.star.in | awk -v p=$p -v r=$ref '{print "STAR --runThreadN", p, "--genomeDir", r, "--readFilesIn", $1, "--readFilesCommand zcat --outSAMstrandField intronMotif --outFileNamePrefix", $2"/", ">"$2"/star.out 2>"$2"/star.err";}' >star.pat.single.cmds.sh

screen -S star.single
sh star.mat.single.cmds.sh
sh star.pat.single.cmds.sh
#status: fin

cat singlecell.sra.files | xargs -I% dirname % | xargs -I% basename % | sort >singlecell.srr
comm -2 -3 singlecell.srr singlecell.15.srr >singlecell.left.srr

###
#filter uniquely mapped
###
cat singlecell.left.srr | xargs -I% echo ${projdir}'/star_hg19_na12878/maternal/'%'/Aligned.out.sam' >mat.sam.files
cat singlecell.left.srr | xargs -I% echo ${projdir}'/star_hg19_na12878/paternal/'%'/Aligned.out.sam' >pat.sam.files
cat mat.sam.files pat.sam.files >sam.files
cat sam.files | xargs -I% echo awk -F"'\t'" "'\$1 ~ /^@/ || \$5 == 255 {print \$0;}'" % '| samtools view -b - | samtools sort -O bam -T '%' - >'%'.unique.bam' >filter.multi.cmds.sh
cat filter.multi.cmds.sh | sed 's/Aligned.out.sam.unique.bam/Aligned.out.unique.bam/' >file.tmp
mv file.tmp filter.multi.cmds.sh
screen -S filter.multi
cat filter.multi.cmds.sh | parallel -n 1 -P 26
#status: fin

#rm six non-singlecell-samples
\ls -1 ${projdir}/star_hg19_na12878/maternal | sort >all.srr
comm -3 all.srr singlecell.srr >nonsingle.srr
cat nonsingle.srr | xargs -I% echo ${projdir}/star_hg19_na12878/maternal/% >mat.nonsingle.dirs
cat nonsingle.srr | xargs -I% echo ${projdir}/star_hg19_na12878/paternal/% >pat.nonsingle.dirs
cat mat.nonsingle.dirs pat.nonsingle.dirs >nonsingle.dirs
cat nonsingle.dirs | xargs -I% rm -r %
find ${projdir}/star_hg19_na12878 -name '*.bam' | xargs -I% ls -lrth %
find ${projdir}/star_hg19_na12878 -name '*.bam' | xargs -I% ls -lrth % | wc -l #56. OK.

#paired-end single-cells
comm paired.only.srr singlecell.srr


######################
#Merge alignments (keeping read with least mismatches)
######################
#SEE ase.bl6cast.sh


################################################################
#Allelecounts
################################################################
#SEE snpase.sh


##############
#OBSOLETE
##############

###
#rnaseq.py pipeline
###
exp='marinov_SRP018838'
timestamp='20160304'

execdir='/mnt/kauffman/sandberglab/pipeline3.0'
conf=${HOME}'/cloud/btsync/work/rspd/data/annot/rpkmforgenes/rnaseq.kauffman.pipeline3.0.20150228.conf'
projdir='/mnt/kauffman/sandberglab/pipeline3.0/rnaseq/hsa/'${exp}
logdir=${projdir}'/logs'
scriptdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/scripts/marinov'
p='34'
ref='hg19'

#change dir
cd $scriptdir

#symlink marinov fastq files to projdir
mkdir -p ${projdir}/rawdata
mkdir $logdir
ls -1 $fastqdir | sed 's/_.*//' | sed 's/.fast.*//' >files.srr
ls -1 $fastqdir | awk -v p=${fastqdir} '{print p"/"$0;}' >fastq.files
paste files.srr fastq.files >fastq2srr.files

#mkdir
cat files.srr | sort | uniq | xargs -I% mkdir ${projdir}/rawdata/%

#symlink
cat fastq2srr.files | awk -F'\t' -v p=${projdir}/rawdata '{print "ln -s", $2, p"/"$1"/."}' >ln.cmds.sh
sh ln.cmds.sh

#write pipecmd
echo cd $execdir >rpkmforgenes.sh
echo python3 rnaseq.py -MR -p $p -a $ref -e $exp -c $conf '>'${logdir}/rnaseqpipe.${timestamp}.out '2>' ${logdir}/rnaseqpipe.${timestamp}.err >>rpkmforgenes.sh

#run
screen -S ${exp}.rnaseqpipe
sh rpkmforgenes.sh
#status: fin
