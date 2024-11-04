# Honours-Thesis

## The Effect of Broad-Spectrum Antibiotics on the Nasopharyngeal Microbiome in Premature Infants

### Author: Jonathan van Druten

### PIPELINE
- The pipeline for runnning the analysis is written in bash in `pipeline.sh`
- This script is written for use on High Performance Clusters, specifically the Stellenbosch HPC
- Usage: `qsub pipeline.sh`
- File structure and the location from where the script is run, are essential for the program to work. The script must be run in the same directory as the reads directory, containing a forward and reverse read for each sample of the patient
    - For example the file structure should look like this:
    - `pipeline.sh`
    - reads
        - t1
            - `R1.fastq.gz`
            - `R2.fastq.gz`
        - t2
            - `R1.fastq.gz`
            - `R2.fastq.gz`
        - etc.

### DATA VISUALISATION
`visualisation.R` was used for visualising the results of the pipeline and performing statistical analyses 

### SUPPLEMENTARY MATERIALS
The supplementary tables for figure 5 can be found in the Supplementary Figures folder

### FULL QUALITY FIGURES
The full quality figures from the paper can be found in the Figures folder

### METADATA
Scripts for handling metadata, and metadata itself can be found in the Metadata-scripts folder.

### ORIGINAL RUN INFORMATION
This can be found in the supplementary materials of the paper by Dhariwal et al (2024).


#### REFERENCES
1. Dhariwal, A., Rajar, P., Salvadori, G., Åmdal, H.A., Berild, D., Saugstad, O.D., Fugelseth, D., Greisen, G., Dahle, U., Haaland, K. & Petersen, F.C. 2024. Prolonged hospitalization signature and early antibiotic effects on the nasopharyngeal resistome in preterm infants. Nature Communications, 15(1):6024, doi: 10.1038/s41467-024-50433-7

2. Lu, J., Rincon, N., Wood, D.E., Breitwieser, F.P., Pockrandt, C., Langmead, B., Salzberg, S.L. & Steinegger, M. 2022. Metagenome analysis using the Kraken software suite. Nature Protocols, 17(12):2815-2839, doi: 10.1038/s41596-022-00738-y