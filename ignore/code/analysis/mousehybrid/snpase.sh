

#Params
projdir='/mnt/kauffman/edsgard/nobackup/data/internal/mmu_sc_fibro/scphaser_set'
scriptdir='/mnt/kauffman/edsgard/cloud/btsync/work/rspd/projects/scphaser/nogit/scripts/mousehybrid'
snpfile='/mnt/kauffman/edsgard/rsync/work/rspd/data/external/sanger_mouse_proj/split/mm9.refseq.pos'
ref='refseq'

#symlink 336 bamfiles
cd $scriptdir
cat sample2id2bamfile | cut -d' ' -f3 | xargs -I% ln -s % ${projdir}/star_merged_cast1_mm9/.


######
#Down-sample variants (mm9.refseq.pos)
######

N='787883'
N='2500000'
N='5000000'
N='1000000'
shuf -n $N $snpfile > ${snpfile}.nrand_${N}


#################################
#Get ASE
#################################

###############
#Setup files and folders
###############

#Params
rundir=${projdir}/'ase'/${ref}/downsampvars_${N}
srcdir='/mnt/crick/edsgard/src/rnaseq/ase/danielr/code'
confdir='/mnt/crick/edsgard/src/rnaseq/ase/danielr/conf'

#Mkdirs
mkdir -p $rundir
cd $rundir

mkdir c57
mkdir c57_cast
mkdir cast
mkdir cast_c57

#Symlink python scripts (not from srcdir since that one doesnt include -Q 30 in the mpileup command)
ln -s ${projdir}/'ase'/${ref}/snp_stats2.py

#Symlink conf files
ln -s ${confdir}/samples.conf

#symlink snp-file
ln -s ${snpfile}.nrand_${N} castsnps.txt

#bamfiles per genotype (dummies)
touch cast_bamfiles.txt
touch c57_bamfiles.txt
touch c57_cast_bamfiles.txt


###############
#Print cmds
###############

###
#Params
###
#IN
bamref='star_merged_cast1_mm9'

#OUT
logdir=${projdir}'/log/ase/'${ref}/N_${N}
cmds=ase.${ref}.N_${N}.cmds.sh
bamlist='cast_c57_bamfiles.txt' #out from make_bam_list.py
##allelecounts=${rundir}/'allelecounts.tab' #out from snp_stat2.py
##inf=${rundir}/'gene_by_cells_cast_c57_reads.txt' #out from snp_stat2.py


###
#Print cmds
###
#Mkdirs
mkdir -p $logdir
cd $scriptdir

echo cd $rundir >${cmds}

#Make bam list (same files as before, so just symlinking)
ln -s ${projdir}'/ase/'${ref}/${bamlist} ${rundir}/.
##echo python ${srcdir}/make_bam_list.py ${projdir}/${bamref}/* '>'${bamlist} >>${cmds}

#mpileup and nucleotide-count (but having manually added -Q 30 to mpileup to be the same as for the human data)
echo python3 snp_stats2.py -MS -p 30 '1>'${logdir}/snp_stats.out '2>'${logdir}/snp_stats.err >>${cmds}
#excluded -F


###############
#Execute
###############
screen -S ase.${ref}.N_${N}
ref='refseq'
N=5000000
sh ase.${ref}.N_${N}.cmds.sh
#status: sub: N=1M, 2.5M, 5M (13.28)

#check logs
N=1000000
N=2500000
N=5000000
rundir=${projdir}/'ase'/${ref}/downsampvars_${N}
logdir=${projdir}/log/ase/${ref}
lt ${logdir}/N_${N}
cat ${logdir}/N_${N}/*.err

#check mpileup
lt $rundir/cast_c57/mpileups
lt $rundir/cast_c57/mpileups | wc -l

#check allelecounts
lt $rundir/cast_c57/cellsums
lt $rundir/cast_c57/cellsums | wc -l


#############
#nucleotide counts 2 allele counts
#############

#list of ac-files
N=1000000
N=2500000
N=5000000
rundir=${projdir}/'ase'/${ref}/downsampvars_${N}
find ${rundir}/cast_c57/cellsums -type f -name '*' >ac.N_${N}.files
head ac.*.files
wc -l ac.*.files

#SEE ac2mat.R


##################
#CHECKS
##################
emacs -nw ${projdir}/ase/refseq/snp_stats2.py
cat ${confdir}/samples.conf

annotdir='/mnt/crick/edsgard/src/rnaseq/ase/danielr/annot/snps/mm9/cast'

#validated_mm9_refseq_snp2genes.txt
head ${annotdir}/validated_mm9_refseq_snp2genes_nooverlap.txt
#tab-sep cols: gene, n.snps, snps
#probably as easy to use script I already have made for the human data as to create the snp2genes.txt formatted annotation file.
#SEE: marinov::ac2knownhap.R or borel::ac2mat.R

#validated_cast_c57_snps.txt
head ${annotdir}/validated_cast_c57_snps.txt
#tab-sep cols: chr, pos, c57, cast 

emacs -nw ${srcdir}/summarize_cell_allelereads.py
#validated_cast_c57_snps.txt:
#-only first four cols used
#-third col is c57 and fourth cast

