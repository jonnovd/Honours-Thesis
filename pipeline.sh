#!/bin/bash

#PBS -N bracken
#PBS -l walltime=6:00:00
#PBS -l ncpus=32
#PBS -l mem=150GB
#PBS -e pipe.err
#PBS -o pipe.out
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
THREADS=32
PATIENT=p096
BASE_DIR=~/thesis/preterm/$PATIENT
READS_DIR=$BASE_DIR/reads

for dir in "$READS_DIR"/*; do
    echo "$dir" >> samples.txt
done

# FASTQC/FASTP - Quality control
# q = Phred Quality threshold
# l = Minimum read length
# Y = Minimum % sequence complexity
# w = Threads
echo ">> QC WITH FASTP"
module load app/fastp/0.23.4
if [[ ! -d fastp ]]; then
    mkdir fastp
    mkdir fastp-out
fi
for dir in "$READS_DIR"/*;do
    R1=$(basename "$(ls $dir/*R1*)")
    R2=$(basename "$(ls $dir/*R2*)")
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
    
    mv fastp.json fastp-out/fastp-$PATIENT-$SAMPLE.json
    mv fastp.html fastp-out/fastp-$PATIENT-$SAMPLE.html

done
module unload app/fastp/0.23.4

QC_READS_DIR=$BASE_DIR/fastp

echo ">> HOST REMOVAL WITH BOWTIE2"
module load app/bowtie/2.5.4
if [[ ! -d '../GRCh38_noalt_as' ]]; then
    mkdir ../GRCh38_noalt_as
    bowtie2-build $HOST_FASTA ../GRCh38_noalt_as --threads $THREADS
    mv *.bt2 ../GRCh38_noalt_as
fi
if [[ ! -d 'bowtie2' ]]; then
    mkdir bowtie2
fi
for dir in "$QC_READS_DIR"/*;do
    R1=$(ls $dir/*R1*)
    R2=$(ls $dir/*R2*)
    SAMPLE=$(awk -F'/' '{print $NF}' <<< $dir)
    if [[ ! -d "bowtie2/$SAMPLE" ]]; then
        mkdir "bowtie2/$SAMPLE"
    fi

    bowtie2 -p 16 -x ../GRCh38_noalt_as/GRCh38_noalt_as \
        -1 $R1 \
        -2 $R2 \
        --very-sensitive-local \
        --un-conc-gz \
        "bowtie2/$SAMPLE/${SAMPLE}_host_removed" > "bowtie2/$SAMPLE/${SAMPLE}_mapped_and_unmapped.sam"
    
    mv bowtie2/$SAMPLE/${SAMPLE}_host_removed.1 bowtie2/$SAMPLE/${SAMPLE}_host_removed_R1.fastq.gz
    mv bowtie2/$SAMPLE/${SAMPLE}_host_removed.2 bowtie2/$SAMPLE/${SAMPLE}_host_removed_R2.fastq.gz
done
module unload app/bowtie2/2.5.3

echo ">> TAXONOMIC PROFILING WITH KRAKEN2"
if [[ ! -d 'kraken2' ]]; then
    mkdir kraken2
fi
module load app/kraken2/2.1.3
for dir in bowtie2/*;do
    R1=$(ls $dir/*R1*)
    R2=$(ls $dir/*R2*)
    SAMPLE=$(awk -F'/' '{print $NF}' <<< $dir)
    if [[ ! -d "kraken2/$SAMPLE" ]]; then
        mkdir "kraken2/$SAMPLE"
    fi

    kraken2 --db $KRAKEN2_DB \
        --paired \
        --classified-out kraken2/$SAMPLE/hr_${SAMPLE}_cseqs#.fq \
        --unclassified-out kraken2/$SAMPLE/hr_${SAMPLE}_useqs#.fq $R1 $R2 \
        --report kraken2/$SAMPLE/hr-$PATIENT-$SAMPLE-kreport.txt \
        --threads $THREADS
done

echo ">> ABUNDANCE RE-ESTIMATION WITH BRACKEN"
if [[ ! -d 'bracken' ]]; then
    mkdir bracken
fi
module load app/bracken/2.7
for dir in kraken2/*;do
    SAMPLE=$(awk -F'/' '{print $NF}' <<< $dir)
    if [[ ! -d "bracken/$SAMPLE" ]]; then
        mkdir "bracken/$SAMPLE"
    fi

    bracken -d $KRAKEN2_DB \
        -i kraken2/$SAMPLE/hr-$PATIENT-$SAMPLE-kreport.txt \
        -o bracken/$SAMPLE/hr-$PATIENT-$SAMPLE-report.bracken \
        -r 150
done
module unload app/bracken/2.7

echo ">> REPORT WITH MULTIQC"
module load app/multiqc/1.21
multiqc kraken2/. -o kraken2
multiqc fastp-out/. -o fastp-out

echo ">> ARG analysis with ABRicate"
module load python/2.7.11
module load app/megahit/1.0.4-beta-58b0995
module load app/abricate/1.0.0
module load app/NCBI/2.15.0+
for dir in bowtie2/*;do
    R1=$(ls $dir/*R1*)
    R2=$(ls $dir/*R2*)
    SAMPLE=$(awk -F'/' '{print $NF}' <<< $dir)
    if [[ ! -d "abricate/$SAMPLE" ]]; then
        mkdir "abricate/$SAMPLE"
    fi
    megahit -o abricate/$SAMPLE/mega-out -1 $R1 -2 $R2
    
    abricate abricate/$SAMPLE/mega-out/*.fa > results.tab
    abricate --summary results.tab > summary.tab
done