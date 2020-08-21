#!/bin/bash

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

# Created on Wed May 06 13:45 2020
# @author: Maria del Carmen Martos Contreras

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

# VARIABLES

SIM=#Put references DNA or RNA simulation path
REF=#Put reference fasta genomes path
TRIMMOMATIC= #Put your trimmomatic file version with the path
BWA= #Put your bwa execute file with it path
SAMTOOLS= #Put the samtools execute file path
OUT=${PWD}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

# SIMULATION

SIM_FILE=$1
COVERAGE=$2
READS=$3

mkdir ART_SIMULATION
mkdir ${OUT}/ART_SIMULATION/ART_FILES
art_illumina -ss MSv1 -na -i ${SIM_FILE} -p -l 150 -f ${COVERAGE} -m 152 -s 10 -o ${OUT}/ART_SIMULATION/ART_FILES/$(basename ${SIM_FILE%.*})

mkdir ISS_SIMULATION
mkdir ${OUT}/ISS_SIMULATION/SINGLE_FASTA
mkdir ${OUT}/ISS_SIMULATION/ISS_FILES
python3 multifasta_converter.py -N 150 -o ${PWD}/ISS_SIMULATION/SINGLE_FASTA/$(basename ${SIM_FILE}) -f ${SIM_FILE}
iss generate -p 12 -g ${OUT}/ISS_SIMULATION/SINGLE_FASTA/$(basename ${SIM_FILE}) -n ${READS} -m novaseq -o ${OUT}/ISS_SIMULATION/ISS_FILES/$(basename ${SIM_FILE%.*})

wait

############################################################################################################################################################
# CHECK QUALITIES
############################################################################################################################################################

# fastqc tool to check fastq reads qualities.
#mkdir ./LIBRARIES_FASTQC

#${FASTQC} 01_ISS_REF/*.fastq
#mv ./01_ISS_REF/*.fastqc* ./REF_FASTQC/

############################################################################################################################################################
# FILTER SAMPLE DATA
############################################################################################################################################################

# Trim at 200bp length and remove reads shorter than 140bp while keeping paired-end reads
# Tool trimmomatic

##Make ART_TRIM directory
mkdir ${OUT}/ART_SIMULATION/02_ART_TRIM

NAME=$(basename ${SIM_FILE%.*})
RUN1=${OUT}/ART_SIMULATION/ART_FILES/${NAME}1.fq
##Call trimmomatic
java -jar ${TRIMMOMATIC} SE -threads 12 -phred33 ${RUN1} ${OUT}/ART_SIMULATION/02_ART_TRIM/${NAME}_R1_filtered.fastq CROP:150 MINLEN:140

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

##Make ISS_TRIM directory
mkdir ${OUT}/ISS_SIMULATION/02_ISS_TRIM
RUN1=${OUT}/ART_SIMULATION/ART_FILES/${NAME}_R1.fq
##Call trimmomatic
java -jar ${TRIMMOMATIC} SE -threads 12 -phred33 ${RUN1} ${OUT}/ISS_SIMULATION/02_ISS_TRIM/${NAME}_R1_filtered.fastq CROP:150 MINLEN:140

wait

############################################################################################################################################################
# MAPPING SAMPLE TO EACH REFERENCE GENOME
############################################################################################################################################################

#Aligne each read to the genome with BWA

${BWA}/bwa index ${REF}/${FILE}

##make ART_SAM directory
mkdir ${OUT}/ART_SIMULATION/03_ART_SAM_${NAME}

##Aligne each read to fasta file and obtain result in sam format
for INDEX in ${REF}/*.fna
do
        REF_N=$(basename ${INDEX%.*})
        ${BWA}/bwa mem -t 12 ${INDEX} ${OUT}/ART_SIMULATION/02_ART_TRIM/${NAME}_R1_filtered.fastq  > ${OUT}/ART_SIMULATION/03_ART_SAM_${NAME}/${REF_N}_${NAME}.sam
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

##make ISS_SAM directory
mkdir ${OUT}/ISS_SIMULATION/03_ISS_SAM_${NAME}

##Aligne each read to fasta file and obtain result in sam format
for INDEX in ${REF}/*.fna
do
        REF_N=$(basename ${INDEX%.*})
        ${BWA}/bwa mem -t 12 ${INDEX} ${OUT}/ISS_SIMULATION/02_ISS_TRIM/${NAME}_R1_filtered.fastq  > ${OUT}/ISS_SIMULATION/03_ISS_SAM_${NAME}/${REF_N}_${NAME}.sam
done
wait

############################################################################################################################################################
# FILTER SAM FILE
############################################################################################################################################################

# Remove supplementary alignments, not primary alignments and not mapped reads.

## Create a filter mapped samples directory
mkdir ${OUT}/ART_SIMULATION/04_ART_FILTERMAP_${NAME}

## Filter sam files
for SAM in ${OUT}/ART_SIMULATION/03_ART_SAM_${NAME}/*.sam
do
        FILTER=$(basename ${SAM})
        ${SAMTOOLS} view -@ 12 -S -F2308 ${SAM} -o ${OUT}/ART_SIMULATION/04_ART_FILTERMAP_${NAME}/${FILTER}
done

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

## Create a filter mapped samples directory
mkdir ${OUT}/ISS_SIMULATION/04_ISS_FILTERMAP_${NAME}

## Filter sam files
for SAM in ${OUT}/ISS_SIMULATION/03_ISS_SAM_${NAME}/*.sam
do
        FILTER=$(basename ${SAM})
        ${SAMTOOLS} view -@ 12 -S -F2308 ${SAM} -o ${OUT}/ISS_SIMULATION/04_ISS_FILTERMAP_${NAME}/${FILTER}
done

############################################################################################################################################################
# QUANTIFICATION AND IDENTIFICATION
############################################################################################################################################################


python3.6 ${GD}/g-d_algorithm_3.6.py -r ${RUN_ART} -g 0.99 -d 0.98 -m ${OUT}/ART_SIMULATION/04_ART_FILTERMAP_${NAME} -o ${OUT}/ART_SIMULATOR/art_${NAME}.csv -O ${OUT}/ART_SIMULATOR/genomic_regions_${NAME}.csv

#-----------------------------------------------------------------------------------------------------------------------------------------------------------

python3.6 ${GD}/g-d_algorithm_3.6.py -r ${RUN_ART} -g 0.99 -d 0.98 -m ${OUT}/ISS_SIMULATION/04_ISS_FILTERMAP_${NAME} -o ${OUT}/ISS_SIMULATOR/iss_${NAME}.csv -O ${OUT}/ISS_SIMULATOR/genomic_regions_${NAME}.csv

