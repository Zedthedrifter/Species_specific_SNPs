#!/bin/bash
#SBATCH --job-name="SNP analysis" 
#SBATCH --export=ALL
#SBATCH --partition=medium
#SBATCH --mem=32G 
#SBATCH --cpus-per-task=8


function setup_evn {

conda create -n $ENV -y
conda install -c bioconda samtools -y --name $ENV
conda install -c bioconda bcftools -y --name $ENV
conda install -c bioconda vcftools -y --name $ENV
conda install -c bioconda pandas -y --name $ENV #FOR THE EXTRACTION STAGE
}

#STEP1
#MAKE VCF FROM BAM FILES
function bam_to_vcf {

BAMLIST=$1 #THE DIRECTORIES OF THE SAMPLE BAM FILES TO BE PROCESSED: TRAINING SET AS PARENT OR HYBRIDS ETC. 
OUTDIR=$2
ref=$3

#Prep: Converting bam files to VCF files which can then be used as input for the NucBarcoder pipeline
#COPY THE REFERENCE FASTA
rm $INPUTS/REFERENCE.fa* #YOU CAN COMMENT THIS OUT AFTER THE FIRST TIME
cp $ref $INPUTS/REFERENCE.fa #YOU CAN COMMENT THIS OUT AFTER THE FIRST TIME

bcftools mpileup -Ob -f $INPUTS/REFERENCE.fa --bam-list $BAMLIST | \
    bcftools call -Ob -mv| \
    bcftools filter -Ob -s LowQual -e 'QUAL<30 || DP<100' > $OUTDIR/var.flt1.bgzf #filter in the next step. no need to add redundant step?? (it's nice to see PASS/Low Qual though)

#FILTERING
#--bcf can read bgzf input as well
vcftools --bcf $OUTDIR/var.flt1.bgzf \
         --out $OUTDIR/var.flt2 \
         --recode --recode-INFO-all \
         --minQ 30 \
         --max-missing 0.95 \
         --minDP 7  \
         --min-alleles 2 \
         --max-alleles 2 \
         --remove-indels \
         --hwe 0.05
#
#--maxDP 100
}

#STEP2
#YOU NEED TO MAKE A CSV FILE FOR SAMPLE NAMES (JUST COPY THE bamlist.txt FILE IN YOUR CURRENT DIRECTORY AND $RESULT1) AND SPECIES IN $RESULT1, NAME IT ID_to_scientific_name.csv
function extract_loci {

INDIR=$1
OUTDIR=$2
NAME=$3
high=$4
low=$5

#./calculate_snp.freq.py specifi_SNPs \
./species_specific_allele.py specifi_SNPs \
        -v $INDIR/var.flt2.recode.vcf \
        -n $NAME \
        -o $OUTDIR #output directory

}


#FOR A GIVEN SAMPLE, CHECK SPECIES SPECIFIC SNPS FROM THE TRAINING DATASET
function species_classifier {

VCF=$1
OUTDIR=$2
PARENT_SNP=$3
high=$4
low=$5

#./calculate_snp.freq.py find_species \
./species_specific_allele.py find_species \
        -v $VCF \
        -n $PARENT_SNP \
        -o $OUTDIR #output directory

}
#===================================================================

function main {

#parameters and paths to adjust##########################################################################################################################################

#UNIVERSAL VARIABLES
ENV=snps 
WORKDIR=$SCRATCH/SNPS/hybrids #the first project: WORKDIR=$SCRATCH/SNPS
#high=50 #LOWEST FREQUENCY OF THE 'PRESENT' SNP
#low=0 #HIGHEST FREQEUENCY OF THE 'ABSENT' SNP


#CONSTANTS PROVIDED BY USER
REF=/mnt/shared/projects/rbge/A_projects_Markus/Araucaria/Lib2_mydata_extracted/Araucaria_input-seq_with400Ns_beginend.fas #the fasta file for mapping
NAME=/mnt/shared/projects/rbge/A_projects_Markus/Araucaria/bamlist_102.csv #SAMPLE BAM & CORRESPONDING SPECIES
PARENTS='FULL PATH TO THE LIST OF BAMS AS TRAINING SET FOR SPECIES SPECIFIC SNPS'#TRAINING SET FOR SPECIES SPECIFIC SNPS
SAMPLES=/mnt/shared/scratch/pholling/SNPS/hybrids/results/01_raw_data/hybrids.txt #UNCLASSIFIED SAMPLES

#End of parameters and paths to adjust#####################################################################################################################################

#subdirectories
INPUTS=$WORKDIR/results/00_inputs #JUST A PLACE TO PUT COPIES OF REF FILE AND THE FAI INDEXING
RESULT1=$WORKDIR/results/01_raw_data #WHERE THE VCF FILES GO
RESULT2=$WORKDIR/results/02_output #THE SNP STATS AND SELECTED SNPS
RESULT3=$WORKDIR/results/03_sample_vcf
RESULT4=$WORKDIR/results/04_sample_SNPs

function setup_workdir {

mkdir $WORKDIR
mkdir $WORKDIR/results
mkdir $INPUTS
mkdir $RESULT1
mkdir $RESULT2
mkdir $RESULT3
mkdir $RESULT4
}
#RUN COMMANDS
#SET UP YOUR ENV AND WORK DIRECTORIES
#setup_evn #IF YOU HAVE AN ENV WITH ALL THE REQUIRED TOOLS INSTALLED (SEE THE FUNCTION), YOU CAN SKIP THIS AND ACTIVATE THE CORRESPONDING ENV
setup_workdir

#CONDA ACTIVATE snps!!!!!!
#conda activate snps BEFORE PROCEEDING!!! THAT'S WHY I COMMENTED OUT THE REST OF THE CODES! 

#STEP1
#bam_to_vcf $PARENTS $RESULT1 $REF

#STEP 2: THIS IS FAST AND YOU DON'T EVEN NEED TO SBATCH IT
extract_loci $RESULT1 $RESULT2 $NAME 

#STEP 3: MAKE VCF FILES OF THE QUERIES
#bam_to_vcf $SAMPLES $RESULT3 $REF

#STEP 4: IDENTIFY SPECIES BASED ON SNPS
#list of inputs: sample_vcf, output dir, species specific SNPs csv, $high, $low
#species_classifier $RESULT3/var.flt2.recode.vcf $RESULT4 $RESULT2/hq_specific_allele_freq.csv $high $low
}

main


