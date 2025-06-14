#entrare nel percorso dove sono presenti le cartelle
##assicurarsi di essere in /home/gitpod
##pwd /home/gitpod

cd datiesame/
ls -l

##pwd /home/gitpod/datiesame

mkdir -p analysis
cd analysis


mkdir reads
cd reads

tar -xvzf /home/gitpod/datiesame/data_rnaseq.tar.gz -C .

### check if you can execute salmon by typing "salmon" on the terminal
## sometimes it fails on RStudio terminal, on CodeSpaces while it works on GitPod
## export PATH=${PATH}:/usr/local/bin

## the index for the transcriptome is located in
## /home/gitpod/datiesame/datasets_reference_only/trascriptome/chr21_transcripts_index

## now we can quantify all samples, by running a loop with salmon and the following
##essere in reads
###cd ..
###cd analisi/reads/
###con salmon vado a quantificare, cioè stimare l’abbondanza dei trascritti

for sample in `ls *_1.fasta.gz`
do
index="/home/gitpod/datiesame/datasets_reference_only/trascriptome/chr21_transcripts_index"
name=${sample%_1.fasta.gz}
echo "quantifying $name"
salmon quant \
 -p 2 \
 -i $index \
 -l IU \
 -1 "${name}_1.fasta.gz" -2 "${name}_2.fasta.gz" \
 --validateMappings \
 -o "${name}.quant"
echo -e "$name done now\n"
done

### let's inspect a quantification file

cd sample_01.quant
head quant.sf

## more information on the format of the output
## https://salmon.readthedocs.io/en/latest/file_formats.html
