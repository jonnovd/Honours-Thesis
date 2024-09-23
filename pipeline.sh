#!/bin/bash

#PBS -N profiler
#PBS -l walltime=12:00:00
#PBS -l ncpus=16
#PBS -l mem=150GB
#PBS -e pipe.err
#PBS -o pipe.out
#PBS -q bix
#PBS -M 25092936@sun.ac.za 
#PBS -m ae

cd $PBS_O_WORKDIR

# if [[ -e 'module_versions.txt' ]]; then 
#     rm module_versions.txt 
# fi
# module avail app/fastp >> module_versions.txt 2>&1
# module avail app/bowtie >> module_versions.txt 2>&1
# module avail app/kraken2 >> module_versions.txt 2>&1
# module avail app/bracken >> module_versions.txt 2>&1
# module avail app/multiqc >> module_versions.txt 2>&1


KRAKEN2_DB='/new-home/databases/kraken.db/standard'
HOST_FASTA=''
MQC_DIR='.'
THREADS=16
#READS_DIR=~/thesis/xiao/reads
READS_DIR=~/thesis/pipeline-test/reads

for dir in "$READS_DIR"/*; do
    echo "$dir" >> samples.txt
done

# FASTQC/FASTP - Quality control
# q = Phred Quality threshold
# l = Minimum read length
# Y = Minimum % sequence complexity
# w = Threads
echo ">> QC with fastp"
module load app/fastp/0.23.4
if [[ ! -d fastp ]]; then
    mkdir fastp
    mkdir fastp-out
fi
for dir in "$READS_DIR"/*;do
    R1=$(basename "$(ls $dir/*R1.*)")
    R2=$(basename "$(ls $dir/*R2.*)")
    SAMPLE=$(awk -F'/' '{print $NF}' <<< $dir)
    if [[ ! -d "fastp/$SAMPLE" ]]; then
        mkdir "fastp/$SAMPLE"
    fi

    fastp -i "$dir/$R1" \
        -I "$dir/$R2" \
        -o "fastp/$SAMPLE/qc_$R1" \
        -O "fastp/$SAMPLE/qc_$R2" \
        -q 30 \
        -l 50 \
        -y -Y 50 \
        -w $THREADS
        # --detect_adapter_for_pe \ Check if this is necessary because it always seems to fail to detect
    
    mv fastp.json fastp-out/fastp-$SAMPLE.json
    mv fastp.html fastp-out/fastp-$SAMPLE.html

done
module unload app/fastp/0.23.4

QC_READS_DIR=~/thesis/pipeline-test/fastp

echo ">> Host removal with Bowtie2"
module load app/bowtie/2.5.4
if [[ ! -d 'GRCh38_noalt_as' ]]; then
    mkdir GRCh38_noalt_as
    bowtie2-build $HOST_FASTA GRCh38_noalt_as --threads $THREADS
    mv *.bt2 GRCh38_noalt_as
fi

if [[ ! -d 'bowtie2' ]]; then
    mkdir bowtie2
fi

for dir in "$QC_READS_DIR"/*;do
    R1=$(ls $dir/*R1.*)
    R2=$(ls $dir/*R2.*)
    SAMPLE=$(awk -F'/' '{print $NF}' <<< $dir)
    if [[ ! -d "bowtie2/$SAMPLE" ]]; then
        mkdir "bowtie2/$SAMPLE"
    fi

    bowtie2 -p 16 -x GRCh38_noalt_as/GRCh38_noalt_as \
        -1 $R1 \
        -2 $R2 \
        --very-sensitive-local \
        --un-conc-gz \
        "bowtie2/$SAMPLE/${SAMPLE}_host_removed" > "bowtie2/$SAMPLE/${SAMPLE}_mapped_and_unmapped.sam"
    
    mv bowtie2/$SAMPLE/${SAMPLE}_host_removed.1 bowtie2/$SAMPLE/${SAMPLE}_host_removed_R1.fastq.gz
    mv bowtie2/$SAMPLE/${SAMPLE}_host_removed.2 bowtie2/$SAMPLE/${SAMPLE}_host_removed_R2.fastq.gz
done
module unload app/bowtie2/2.5.3

if [[ ! -d 'kraken2' ]]; then
    mkdir kraken2
fi

echo ">> Taxonomic Profiling with Kraken2"
module load app/kraken2/2.1.3
for dir in bowtie2/*;do
    R1=$(ls $dir/*R1.*)
    R2=$(ls $dir/*R2.*)
    SAMPLE=$(awk -F'/' '{print $NF}' <<< $dir)
    if [[ ! -d "kraken2/$SAMPLE" ]]; then
        mkdir "kraken2/$SAMPLE"
    fi

    kraken2 --db $KRAKEN2_DB \
        --paired \
        --classified-out $SAMPLE/${SAMPLE}_cseqs#.fq \
        --unclassified-out $SAMPLE/${SAMPLE}_useqs#.fq $R1 $R2 \
        --report kraken2/$SAMPLE/${SAMPLE}_kreport.txt \
        --threads $THREADS
done
if [[ ! -d 'bracken' ]]; then
    mkdir bracken
fi

echo ">> Abundance Estimation with Bracken"
module load app/bracken/2.7
for dir in kraken2/*;do
    SAMPLE=$(awk -F'/' '{print $NF}' <<< $dir)
    if [[ ! -d "bracken/$SAMPLE" ]]; then
        mkdir "bracken/$SAMPLE"
    fi

    bracken -d $KRAKEN2_DB \
        -i ${dir}${SAMPLE}_kreport.txt \
        -o bracken/$SAMPLE/report.bracken \
        -r 150
done
module unload app/bracken/2.7

# echo ">> De novo Assembly with SPAdes"
# module load app/SPAdes
# # SPADES - Assembly of unclassified reads
# spades.py \
#     -1 \
#     -2 \
#     -o
# module unload app/SPAdes

echo ">> Report with Multiqc"
module load app/multiqc/1.21
multiqc $MQC_DIR