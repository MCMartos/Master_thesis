# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:15:39 2020
@author: María del Carmen Martos Contreras
Multifasta converter to fasta
"""

import argparse, os

#Parameters definition
parser = argparse.ArgumentParser( prog='Multifasta to fasta converter', usage='multifasta_converter [options]', description='multifasta_converter: Convert multifasta files to fasta and add Nx to unite sequences')
parser.add_argument("-N",   "--N_number",        type=int, help="Determine an specific number of N to add in sepuences space")
parser.add_argument("-f",  "--multifasta_file",  type=str,   help="Introduce the multifasta file name to convert in a single fasta file")
parser.add_argument("-o",   "--out",            type=str,   help="The prefix of output filename")
args = parser.parse_args()

seqs=[]
fasta = open(args.out,"wt")
with open(args.multifasta_file, "rt") as multifasta:
    prev_seq = []
    for line in multifasta:
        if line.startswith(">"):
            seqs.append("".join(prev_seq))
            prev_seq = []
        else:
            prev_seq.append(line.rstrip())

seqs.append("".join(prev_seq))
name=os.path.basename(args.out)
fasta.write(">"+name+"\n")
fasta.write(("N"*args.N_number).join(seqs)+"\n")
multifasta.close()
fasta.close()
