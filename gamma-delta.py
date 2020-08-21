# gamma-delta.py
# coding=utf-8
#
# Implementation of the gamma-delta algorithm
# This file is part of the https://github.com/LidiaGS/g-d_algorithm.
# Copyright (c) 2019 Universitat Autònoma de Barcelona
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import argparse, csv, operator, os
from decimal import Decimal

## ARGUMENTS DEFINITION
parser = argparse.ArgumentParser()
parser.add_argument('-r',   '--R',          type=str,  required = False, help='Single-end reads sample (FASTQ/FASTA format)')
parser.add_argument('-r1',  '--R1',         type=str,  required = False, help='Separeted paired-end reads sample, R1 (forward) file (FASTA/FASTQ format)')
parser.add_argument('-r2',  '--R2',         type=str,  required = False, help='Separeted paired-end reads sample, R2 (reverse) file (FASTA/FASTQ format)')
parser.add_argument('-m',   '--SAMfolder',  type=str,  required = True,  help='Folder containing the mapped reads files in SAM format (e.g: mapped_reads.sam) [MANDATORY OPTION]')
parser.add_argument('-g',   '--gamma',      type=str,  required = False,  help='Gamma threshold (per one unit) [0-1] (default 0.99)')
parser.add_argument('-d',   '--delta',      type=str,  required = False,  help='Delta threshold (per one unit) [0-1] (default 0.96)')
parser.add_argument('-o',   '--output',     type=str,  required = True,  help='Output file containing the summary information of the mapped reads (CSV format) [MANDATORY OPTION]')
parser.add_argument('-O',   '--outInfo',    type=str,  required = False, help='Output file containing the assignment information of each read (CSV format)')
options = parser.parse_args()

## ABREVIATIONS
# A1:  Highest mapping ratio
# A2:  Second highest mapping ratio

## FUNCTIONS

## Obtaining the mapping ratio (A) read against a reference by using CIGAR and NM from SAM file
def mapping_ratio(read_CIGAR, read_NM):
    num_cigar = 0
    match_mismatch = 0
    insertion = 0
    deletion = 0
    soft_clipped = 0
    hard_clipped = 0
    padding = 0
    skipped = 0
    mismatch = 0
    match = 0
    NM_mismatch = 0
    total_cigar = 0
    map_ratio = 0

    ## Trim CIGAR into parts and join same values:
    for cigar in read_CIGAR: # Looping through elements of CIGAR
        if cigar.isdigit(): # Save numeric values from CIGAR
            num_cigar = num_cigar * 10
            num_cigar = int(cigar) + num_cigar
        else: # Save numeric values from CIGAR to its cathegory (i.e.: Matching/mismatchig, indels, ...)
            if cigar == 'M': ## M: matches & mismatches
                match_mismatch = num_cigar + match_mismatch
                num_cigar = 0
            elif cigar == 'I': ## I: insertion
                insertion = num_cigar + insertion
                num_cigar = 0
            elif cigar == 'D': ## D: deletion
                deletion = num_cigar + deletion
                num_cigar = 0
            elif cigar == 'S': ## S: soft-clipped
                soft_clipped = num_cigar + soft_clipped
                num_cigar = 0
            elif cigar == 'H': ## H: hard-clipped
                hard_clipped = num_cigar + hard_clipped
                num_cigar = 0
            elif cigar == 'P': ## P: Padding
                padding = num_cigar + padding
                num_cigar = 0
            elif cigar == 'N': ## N: Skipped
                skipped = num_cigar + skipped
                num_cigar = 0
            elif cigar == 'X': ## X: Mismatch
                mismatch = num_cigar + mismatch
                num_cigar = 0
            elif cigar == '=': ## =: Match
                match = num_cigar + match
                num_cigar = 0
            elif cigar == '*':
                num_cigar = 0

    ## Length of CIGAR (also reads length)
    total_cigar = match_mismatch + insertion + deletion + soft_clipped + hard_clipped + padding + skipped + mismatch + match

    ## Mapping ratio
    map_ratio = (((total_cigar + deletion + skipped) - float(read_NM) - soft_clipped - hard_clipped - mismatch) / (total_cigar + deletion + skipped))
    return map_ratio

## When working with paried-end reads, add '/1' and '/2' at the end of the R1 and R2 names, respectively
def rename_paired_end(readName, flag):
    if int(flag) & 1 == 1: # # Check paired-end read; if true, this are paired-end reads
        if int(flag) & 128 == 128: # Check second in pair; if true, return R2
            return readName+str('/2')
        else: # R1
            return readName+str('/1')
    return readName

## Check that A1 is above gamma and A2 below delta
def check_distance(ReadDict):
    gamma = options.gamma if options.gamma is not None else 0.99
    delta = options.delta if options.delta is not None else 0.96

    for ReadId in ReadDict:
        if ReadDict[ReadId][0] == 0:
            ReadDict[ReadId][9] = 'Not-Map'
        elif ReadDict[ReadId][0] != 0 and float(ReadDict[ReadId][5]) >= float(delta):
            ReadDict[ReadId][9] = 'A2-above-delta'
        elif ReadDict[ReadId][0] != 0 and float(ReadDict[ReadId][0]) < float(gamma):
            ReadDict[ReadId][9] = 'A1-below-gamma'
        else:
            ReadDict[ReadId][9] = 'Assigned'
    return ReadDict

## Check whether the file is in FASTA or FASTQ format
def check_fastq_fastq(fileName):
    if fileName.endswith('fasta'): # FASTA file
        j = 2
    else: # FASTQ file
        j = 4
    return j

## Saving the mapping information for a given read and reference into the dictionary
def save_map_info(mapRow, mapRef, ReadDict):
    mapRow = mapRow.strip('\n')
    row = mapRow.split('\t')
    ReadId = str(rename_paired_end(row[0],row[1])) # RNAME

    # Create read in the dictionary if the key does not exist yet
    if ReadId not in ReadDict:
        ReadDict[ReadId] = list([0, '0', 'NA', 'NA', 'NA', 0, 'NA', 'NA', '0', 'NA']) # Create and empty tuple, so A1 and A2 are zero.

    # Obtain the mapping ratio of the read in the current row
    NM = [elem[5:] for elem in row if elem.startswith('NM')]
    newRatio = mapping_ratio(row[5], int(NM[0])) # CIGAR, NM are picked

    # New mapping ratio is higher or equal to A1
    if newRatio >= ReadDict[ReadId][0]:
        # Changing A2
        ReadDict[ReadId][5] = ReadDict[ReadId][0] # A2 (A2 to old A1)
        ReadDict[ReadId][6] = ReadDict[ReadId][1] # Reference id of A2 (A2 to old A1)
        ReadDict[ReadId][7] = ReadDict[ReadId][2] # RNAME A2 (A2 to old A1)
        ReadDict[ReadId][8] = ReadDict[ReadId][3] # POS A2 (A2 to old A1)
        # Changing A1
        ReadDict[ReadId][0] = newRatio # A1 (old A1 to new A1)
        ReadDict[ReadId][1] = mapRef # Reference id of A1 (old A1 to new A1)
        ReadDict[ReadId][2] = row[2] # RNAME A1 (old A1 to new A1)
        ReadDict[ReadId][3] = row[3] # POS A1 (old A1 to new A1)
        ReadDict[ReadId][4] = row[1] # FLAG A1 (old A1 to new A1)
    # New mapping ratio is higher than A2
    elif newRatio > ReadDict[ReadId][5]:
        ReadDict[ReadId][5] = newRatio # A2 (old A2 to new A1)
        ReadDict[ReadId][6] = mapRef # Reference id of A2 (old A2 to new A1)
        ReadDict[ReadId][7] = row[2] # RNAME A2 (old A2 to new A1)
        ReadDict[ReadId][8] = row[3] # POS A2 (old A2 to new A1)
    return ReadDict

## Read and load SAM files. Notice that SAM header is ignored.
def load_sam_files(samFiles):
    ReadDict = dict()

    print ('The ', options.SAMfolder,' folder contains ', len(samFiles), ' references')
    print('Loading mapping information...')

    for samFile in samFiles: # Loop through files in folder
        refName = samFile.split('.')[0] # Keep the file name as the reference code
        sam = open(samFile, 'r')
        for row in sam:  # Loop through lines in sam file
            if not row.startswith('@'): # Avoiding headers
                ReadDict = save_map_info(row, refName, ReadDict)
        sam.close()

    ## Check that the A1 and A2 fulfill the gamma and delta thresholds requirements (e.g. A1 above gamma and A2 below delta)
    ReadDict = check_distance(ReadDict)
    return ReadDict

## Create a dictionary for storing the number of reads assigned to each reference or discarded
def summary_dict(ReadDict):
    print('Summarizing data...')
    assignList = [[assignment[9], assignment[1]] for assignment in ReadDict.values()] # Take only the assignment information of each read
    # assignList contains the assignment decision of each read in assignList[0] and the reference of A1 in assignList[1]
    # We subsequently take the reference of the 'Assigned' reads and gamma-delta decision of the not assigned ones
    assignedReads = [assignedRead[1] if assignedRead[0] is 'Assigned' else assignedRead[0] for assignedRead in assignList ] # Take unique assignment elements

    SummDict = dict()
    for elem in list(set(assignedReads)): # Loop through unique assigned elements
        SummDict[elem] = assignedReads.count(elem)

    save_summary_dict_csv(SummDict) ## Save the number of reads assigned to each detected reference
    return

# Obtain header column name
def sample_header():
    if options.R is not None:
        return os.path.basename(options.R)
    elif options.R1 is not None and options.R2 is not None:
        return str(os.path.basename(options.R1))+' + '+str(os.path.basename(options.R2))
    elif options.SAMfolder is not None:
        return options.SAMfolder
    elif options.samFile is not None:
        return options.samFile
    else:
        import datetime
        return str('Sample:')+str(datetime.datetime.now().time())

# Count the number of reads in file (FASTQ/FASTA formats)
def count_reads_in_file(fileName):
    j = check_fastq_fastq(fileName)
    return int(len(open(fileName).readlines( ))/j)

# Return the number of NOT mapped reads
def count_totalReads():
    totalReads = 0
    if options.R is not None:
        totalReads = count_reads_in_file(options.R)
        print('> Initial reads: '+str(totalReads))
    elif options.R1 is not None and options.R2 is not None:
        totalReads = count_reads_in_file(options.R1)*2
        print('> Initial reads: '+str(totalReads))
    return totalReads

## Save on a CSV file the summary information of the detected species together with the total number of reads assigned and their relative species abundance
def save_summary_dict_csv(SummDict):
    totalMapp = sum(SummDict.values())
    totalReads = count_totalReads()
    totalNotMapp = 0 if totalReads is 0 else int(totalReads - totalMapp)
    print('> Mapped reads: '+str(totalMapp))

    # Pop top two elements from column
    rm_by_A1 = SummDict.pop('A1-below-gamma') if 'A1-below-gamma' in SummDict else 0
    rm_by_A2 = SummDict.pop('A2-above-delta') if 'A2-above-delta' in SummDict else 0

    totalAssig = sum(SummDict.values())
    print('> Assigned reads: '+str(totalAssig))

    ## Sort Dictionary by VALUES
    sorted_Tupla = sorted(SummDict.items(), key=operator.itemgetter(1), reverse=True)

    ## Checking if output file already exists, if so, adding a column to the existing file
    if os.path.isfile(options.output):
        # Store configuration file values
        csvinput = open(options.output,'r')
        reader = csv.reader(csvinput)

        reader_header = next(reader) # header row
        reader_NoMapp = next(reader) # reads that did not map
        reader_A1 = next(reader) # reads remove by A1 row
        reader_A2 = next(reader) # reads remove by A2 row

        all_items = []
        for line in reader:
            all_items.append(line)
        csvinput.close()

        csvoutput = open(options.output, 'w')
        writer = csv.writer(csvoutput)
        writer.writerow(reader_header+[str(sample_header())]) # Write sample name
        writer.writerow(reader_NoMapp+[str('Not-mapping-reads'+' ('+str(totalNotMapp)+')')]) # Write Not-Map
        writer.writerow(reader_A1+[str('A1-below-gamma'+' ('+str(rm_by_A1)+')')]) # Write 'A1-below-gamma'
        writer.writerow(reader_A2+[str('A2-above-delta'+' ('+str(rm_by_A2)+')')]) # Write 'A2-above-delta'

        for item in range(max(len(all_items), len(sorted_Tupla))):
            nitem = item + 1
            if nitem > len(all_items): # Input file has LESS rows than the new data
                writer.writerow(["0"]*len(all_items[0])+[str(sorted_Tupla[item][0]+' ('+str(sorted_Tupla[item][1])+'|'+str(round(Decimal(float(sorted_Tupla[item][1])/float(totalAssig)),5))+')')])
            elif nitem > len(sorted_Tupla): # Input file has MORE rows than the new data
                writer.writerow(all_items[item]+["0"])
            else: # Input file has THE SAME number of rows than the new data
                writer.writerow(all_items[item]+[str(sorted_Tupla[item][0]+' ('+str(sorted_Tupla[item][1])+'|'+str(round(Decimal(float(sorted_Tupla[item][1])/float(totalAssig)),5))+')')])
        csvoutput.close()

    # If output file doesn't exist yet, create it and save the mapping information
    else:
        csvoutput = open(options.output, 'w')
        writer = csv.writer(csvoutput)
        writer.writerow([str(sample_header())]) # Write sample name
        writer.writerow([str('Not-mapping-reads'+' ('+str(totalNotMapp)+')')]) # Write Not-Map
        writer.writerow([str('A1-below-gamma'+' ('+str(rm_by_A1)+')')]) # Write 'A1-below-gamma'
        writer.writerow([str('A2-above-delta'+' ('+str(rm_by_A2)+')')]) # Write 'A2-above-delta'

        for item in sorted_Tupla:
            writer.writerow([str(item[0]+' ('+str(item[1])+'|'+str(round(Decimal(float(item[1])/float(totalAssig)),5))+')')])

    csvoutput.close()
    return

# Creates a list containing all reads in the sample
def load_reads_in_sample(fileName):
    j = check_fastq_fastq(fileName)
    nLine = 0
    readList = list()
    sample = open(fileName, 'r')
    for row in sample:  # Loop through lines in sam files
        if nLine % j == 0: # This line has the read's name
            readList.append(row.split()[0][1:])
        nLine += 1
    return readList

# Add not mapped-reads into the dictionary
def add_not_mapped_reads_to_dict(ReadDict, readList):
    for ReadId in readList:
        if ReadId not in ReadDict:
            ReadDict[ReadId] = list([0, '0', 'NA', 'NA', 'NA', 0, 'NA', 'NA', '0', 'Not-map']) # Create and empty tuple, so A1 and A2 are zero.
    return ReadDict

# Save the assignment information of each mapping read in a CSV file
def save_dict_csv(ReadDict):
    if options.R is not None:
        readList = load_reads_in_sample(options.R)
        ReadDict = add_not_mapped_reads_to_dict(ReadDict, readList)
    elif options.R1 is not None or options.R2 is not None:
        readList = load_reads_in_sample(options.R1)
        readList = [ReadId + '/1' for ReadId in readList] # Adds '/1' at the end of each ReadId
        ReadDict = add_not_mapped_reads_to_dict(ReadDict, readList)
        readList = load_reads_in_sample(options.R2)
        readList = [ReadId + '/2' for ReadId in readList] # Adds '/2' at the end of each ReadId
        ReadDict = add_not_mapped_reads_to_dict(ReadDict, readList)

    csv_columns = ['RNAME','A1','REF_A1','SCAF_A1','POS_A1','FLAG_A1','A2','REF_A2','SCAF_A2','POS_A2','ASSIGNMENT']
    with open(options.outInfo, 'w') as csvoutput:
        writer = csv.writer(csvoutput, delimiter=';')
        writer.writerow(csv_columns)
        for key in ReadDict.keys():
            o = ReadDict[key]
            o.insert(0,key)
            writer.writerow(o)
    csvoutput.close()

## MAIN PROGRAMME
def main():

    # Save mapping information (SAM files)
    wd = os.getcwd() # Get working directory
    os.chdir(options.SAMfolder) # Change the working directory
    samFiles = [f for f in os.listdir(options.SAMfolder) if f.endswith('.sam')]

    # Add the mapping information of each read within a SAM in a dicitonary
    ReadDict = load_sam_files(samFiles)

    # Summary each read assignment
    os.chdir(wd) # Turn back to the original working directory
    summary_dict(ReadDict)

    # Save each read assingment information
    if (options.outInfo is not None):
        save_dict_csv(ReadDict)
    return

## MAIN
main()