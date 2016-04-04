

#Study accession: EGAS00001001009
#dataset accession: EGAD00001001083, EGAD00001001084

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

