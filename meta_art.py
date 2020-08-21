# -*- coding: utf-8 -*-
"""
Created on Thu May 28 18:02:28 2020
@author: María del Carmen Martos Contreras
meta_art as an extension for art_illumina. This allows as to generate mixed 
species libraries for metagenomics studies
"""

import argparse, os, glob

#Parameters definition
parser = argparse.ArgumentParser( prog='Metagenomics Simulator for ART', usage='meta_art [options]', description='meta_art: ART simulator extetion fro metagenomics libraries generating')
parser.add_argument("-g",   "--genomes",        type=str, nargs="+",  help="Input all reference files from species to create mixed specie libraries")
parser.add_argument("-n",  "--n_genomes",       type=int,   help="Number of species to generate, must be the same number of genomes introduced")
parser.add_argument("-ss",  "--seqSys",         choices=['MSv1', 'MSv3'], default='MSv1',   help="Reverse (R2) paired-end reads sample (FASTQ/FASTA)")
parser.add_argument("-p",   "--paired",  action='store_true',       help="Indicate a paired-end read simulation or to generate reads from both ends of amplicons")
parser.add_argument("-sam",   "--samout",   action='store_true',    help=" Indicate to generate SAM alignment file")
parser.add_argument("-mp",   "--matepair",    action='store_true',  help="Indicate a mate-pair read simulation")
parser.add_argument("-o",   "--out",            type=str,   help="The prefix of output filename")
parser.add_argument("-l",   "--len",            type=int,   help="The length of reads to be simulated")
parser.add_argument("-a",  "--abundance_file",  type=str,   help="Text file with abundances for each reference file")
parser.add_argument("-r",  "--nreads",          type=int,   help="Total number of reads (integer number)")
parser.add_argument("-m",   "--mflen",          type=int,   help="The mean size of DNA/RNA fragments for paired-end simulations")
parser.add_argument("-s",   "--sdev",           type=int,   help="The standard deviation of DNA/RNA fragment size for paired-end simulations")
parser.add_argument("-na",   "--noALN", action='store_true',        help="Do not output ALN alignment file")
args = parser.parse_args()

##Create a dictionary with the nucleotides of every library
def nucleotides(n_reads, read_len, a_file):
    nuc={}
    file = open(a_file, "rt")
    for line in file:
        line = line.replace("\n","")
        line = line.split(", ")
        nuc[line[0]]=int(float(line[1])*read_len*n_reads)
    file.close()
    return nuc

##Return the coverage for every specie
def coverage(n_genomes, genomes, nucleotide):
    cov={}
    for i in range(0,n_genomes):
        nuc=0
        file=open(genomes[i], "rt")
        for line in file:
            if not "^>" in line:
                line = line.replace("\n","")
                nuc+=len(line)
        file.close()
        cov[genomes[i]]=nucleotide[genomes[i]]/nuc
    return cov


x=nucleotides(args.nreads,args.len,args.abundance_file)
y=coverage(args.n_genomes, args.genomes, x)

output=os.path.dirname(args.out)+"/"

#Send art single especies libraries simulation for the abundance input
if args.paired:
    if args.samout and args.noALN:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -na -sam -i %s -p -l %d -f %s -m %d -s %d -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]], args.mflen, args.sdev, output+args.genomes[i]))
    elif args.samout and not args.noALN:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -sam -i %s -p -l %d -f %s -m %d -s %d -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]], args.mflen, args.sdev, output+args.genomes[i]))
    elif (not args.samout) and args.noALN:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -na -i %s -p -l %d -f %s -m %d -s %d -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]], args.mflen, args.sdev, output+args.genomes[i]))
    else:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -i %s -p -l %d -f %s -m %d -s %d -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]], args.mflen, args.sdev, output+args.genomes[i]))
elif args.matepair:
    if args.samout and args.noALN:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -na -sam -i %s -mp -l %d -f %s -m %d -s %d -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]], args.mflen, args.sdev,  output+args.genomes[i]))
    elif args.samout and not args.noALN:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -sam -i %s -mp -l %d -f %s -m %d -s %d -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]], args.mflen, args.sdev,  output+args.genomes[i]))
    elif (not args.samout) and args.noALN:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -na -i %s -mp -l %d -f %s -m %d -s %d -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]], args.mflen, args.sdev,  output+args.genomes[i]))
    else:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -i %s -mp -l %d -f %s -m %d -s %d -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]], args.mflen, args.sdev, output+args.genomes[i]))
else:
    if args.samout and args.noALN:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -na -sam -i %s -l %d -f %s -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]],  output+args.genomes[i]))
    elif args.samout and not args.noALN:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -sam -i %s -l %d -f %s -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]],  output+args.genomes[i]))
    elif (not args.samout) and args.noALN:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -na -i %s-l %d -f %s -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]],  output+args.genomes[i]))
    else:
        for i in range(0, args.n_genomes):
            os.system("art_illumina -ss %s -i %s -l %d -f %s -o %s" %(args.seqSys, args.genomes[i], args.len, y[args.genomes[i]],  output+args.genomes[i]))

#Abundance output file 
total_reads=0
n_reads={}
for file in glob.glob(output+"*.fa*1.fq"):
    reads=1
    fq=open(file, "rt")
    header=fq.readline()
    seq=fq.readline()
    while seq:
        seq=fq.readline()
        reads+=1
        sign=fq.readline()
        qualities=fq.readline()
    fq.close()
    if args.paired or args.matepair: 
        n_reads[file[:len(file)-4:]]=reads*2
        total_reads+=n_reads[file[:len(file)-4:]]
    else:
        n_reads[file[:len(file)-4:]]=reads
        total_reads+=n_reads[file[:len(file)-4:]]
    
ab_file=open(args.out+"_abundances.txt",'wt') 
for r in n_reads:
    abundance=n_reads[r]/total_reads
    ab_file.write("%s, %s\n"% (r,abundance)) 
ab_file.close()
         
#Merge differnet species libraries into one file

if args.paired or args.matepair:
    os.system("cat %s*.fa*1.fq > %s_R1.fq" % (output, args.out))
    os.system("cat %s*.fa*2.fq > %s_R2.fq" % (output, args.out))
    if not args.noALN:
        os.system("cat %s*.fa*1.aln > %s_R1.aln" % (output, args.out))
        os.system("cat %s*.fa*2.aln > %s_R2.aln" % (output, args.out))
else:
    os.system("cat %s*.fa*.fq > %s.fq" % (output, args.out))
    if not args.noALN:
        os.system("cat %s*.fa*.aln > %s.aln" % (output, args.out))
if args.samout:
    os.system("cat %s*.fa*.sam > %s.sam" % (output, args.out))

os.system("rm %s*.fa*" %(output))


