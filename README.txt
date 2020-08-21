Prerrequisites:

Install python 3.6 -> sudo apt-get install python3.6
Install ART simulator -> https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm
Install InsilicoSeq simulator -> pip install InSilicoSeq https://insilicoseq.readthedocs.io/en/latest/
Install trimmomatic -> sudo apt-get install trimmomatic http://www.usadellab.org/cms/?page=trimmomatic
Install BWA -> http://bio-bwa.sourceforge.net/ 
Install samtools -> http://samtools.sourceforge.net/
Download g-d-algorithm -> https://github.com/LidiaGS/g-d_algorithm

1. Single_lirbaries_workflow.

It is an automatized pipeline of data simulation by ART and InsilicoSeq simulator tools and provides a 
quantitative estimate of species from samples.

Given a DNA sample, it simulates 150bp read length fastq libraries from ART and InSilicoSeq as Illumina 
NGS with a specific coverage and a set of reference genomes, it will trim data allign to ever reference 
genome and identify and quantify species in the library by gamma-delta algorithm. It will generate fastq 
files, sam files and csv file as output.

2. Mixed_libraries_workflow

This workflow is an automatized pipeline for mixed libraries generation and identification by g-d algorithm.

Given more than a DNA sample, it simulates 150bp read length fastq libraries from ART and InSilicoSeq as Illumina 
NGS with a specific coverage and a set of reference genomes, it will trim data allign to ever reference 
genome and identify and quantify species in the library by gamma-delta algorithm. It will generate fastq 
files, sam files and csv file as output.

3. multifasta_converter

Is a python script that convert multifasta files into a fasta with a single header. This make able to manage 
abundance files for InsilicoSeq simulation. If fasta files have only one header, this stepe is not neccessary.

4. meta_art

Is a python script that allows art_illumina simulation deal with more than one sample. It can be consider as
an extention of ART, such as it makes possible simulate metagenomics sampole with ART genomic simulator tool.

 
All reports and feedbacks are highly appreciate. Please report any suggestion on github or by email to
mariadelcarmen.martos@e-campus.uab.cat.


