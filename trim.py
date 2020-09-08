
# -*- coding: utf-8 -*-
"""given a FASTA or FASTQ file to preprocess all it sequences:
                    perform trim operation"""

import sys,re

#Check there are all the operations neede written 
minimun = ["--input", "--output", "--operation"]
for x in minimun:
    if not x in sys.argv:
        print ("Needed argument '%s'" % x)
        exit(1)
if "trim" in sys.argv and (not"--trim-right" in sys.argv or not "--trim-left" in sys.argv):
        print("The operation trim need both '--trim-left' and 'trim-right' elements")
        exit(1)  
        
#Define sys.argv from empty dictionary
arguments = {}
#CReate a dictionarie of arguments and values
i= 1
while i < len(sys.argv):
    #check if the argument is valid
    argument =  sys.argv[i]
    if "--" in argument and (argument != "--input" and argument != "--output" and argument != "--operation" and argument != "--trim-right" and argument != "--trim-left" and argument != "--adaptor"):
        print("Unknown element '%s'" % argument)
        exit(1) 
        break
    #if any previous conditional happens add arguments values
    else:
        arguments[sys.argv[i]]=sys.argv[1+i]
     #check missing values   
    if "--" in sys.argv[i+1]:
        print ("Argument '%s' need a value" % argument)
        exit(1)
    i+=2   

#1: validation input file
def validation_fast_file(input_file): 
    form = input_file
    if ".fasta" in form or "fastq" in form:
        #Open input-file as reading text
        input_file = open(input_file, "rt")   
        #read the first line
        first_line= input_file.readline()  
        #conditional for know the input-file format
        if first_line[0] == "@":#FASTQ format
            return form, input_file
        elif first_line[0] == ">":#FASTA format
            return form, input_file
    else:#WRONG format
        print("You have introduce a wrong input file")
        exit(1)     
 
#2: open output file
def output(output_file):
    form = output_file 
    output_file = open(output_file, "wt")
    return form, output_file

#3: preprocessing
##percentage function
count_n = {"reads":0,"bases":0,"A":0,"C":0,"G":0,"T":0, "N":0,}
def percentage(nucleotide_count, count_n):
    #upper nucleotides
    nucleotide_count = nucleotide_count.upper() 
    #Count total bases processed
    for nucleotide in "ACGTN":
            #Count each nocleotide processed
            count_n[nucleotide]= str(nucleotide_count).count(nucleotide)
 
    #sumatory of each nucleotide processed value
    count_n["bases"] = count_n["A"]+count_n["C"]+count_n["G"]+count_n["T"]+count_n["N"]
    
    #Percentage for sequences which are more than 0 baseselif "--" in sys.argv[i+1]:
    if count_n["bases"] > 0:   
        for N in "ACGTN":
            percent= int((count_n[N]/count_n["bases"])*100)
            count_n[N] = percent    
    #Percentage for sequence which are 0 bases
    else:
        for N in "ACGTN":     
            count_n[N]= 0 
    #Write the percent
    return count_n

#Print summary summary
def make_summary( adaptors, seq_dic, trim_dic):
    space=len(str(seq_dic["bases"]))-len(str(seq_dic["reads"]))
    print("Summary:")
    print(" "*(space-1), seq_dic["reads"], "reads processed")
    print(str(seq_dic["bases"]) + " bases processed ", "(" + str(seq_dic["A"]) + "% A, "+ str(seq_dic["C"])+ "% C, " + str(seq_dic["G"]) + "% G, "+ str(seq_dic["C"]) + "% T, "+ str(seq_dic["N"]) + "% N)")
    if arguments["--operation"] == "trim":
        space=len(str(seq_dic["bases"]))-len(str(trim_dic["bases"]))
        if space == 0:
            print(trim_dic["bases"], "bases trimmed   ", "(" + str(trim_dic["A"]) + "% A, "+ str(trim_dic["C"])+ "% C, " + str(trim_dic["G"]) + "% G, "+ str(trim_dic["C"]) + "% T, "+ str(trim_dic["N"]) + "% N)")
        else:
            print(" "*(space-1), trim_dic["bases"], "bases trimmed   ", "(" + str(trim_dic["A"]) + "% A, "+ str(trim_dic["C"])+ "% C, " + str(trim_dic["G"]) + "% G, "+ str(trim_dic["C"]) + "% T, "+ str(trim_dic["N"]) + "% N)")
    elif arguments["--operation"] == "adaptor-removal":    
        space=len(str(seq_dic["bases"]))-len(str(adaptors))
        print(" "*(space-1),adaptors,"adaptors-found")

###Create trim function
def trim(input_file, output_file, right, left, validation):
    total_count_seq = ""
    total_seq = ""
    count_n["reads"]=0
    trim_n = {"bases":0, "A":0,"C":0,"G":0,"T":0, "N":0,}
    input_file.seek(0)
    #Conditional to make sure the trims aren bigger than sequence or negative numbers
    ## Make sure that right and left trims are numbers
    match1=re.search (r'\d+',str(right))
    match2=re.search (r'\d+',str(left))
    
    #If match1 and match2 are numbers
    if match1 and match2:
        right =int(right)
        left = int(left)  
        #COnditional for negative trim numbers
        if right<0 or left<0:
            print("trims values are negative")
            exit(1)
            
        #while there is lines in the input file it keep reading
        while True:
                tag = input_file.readline()#Tag
                if not tag: break # End of file
                seq = input_file.readline()
                count_n["reads"]+=1
                
                #conditional for length od trims and sequence and if there are negative values
                if right+left < len(seq)-2 and (right > 0 or left > 0):               
                        #Give trim seq
                        trim_seq = seq[left-1:(len(seq)-1)-right:1]
                        output_file.write("%s%s\n" %(tag,trim_seq))
                       
                        #Conditional for FASTQ format
                        if "fastq" in validation:   
                            sign = input_file.readline()#+
                            qualities = input_file.readline()#qualities
                            #write trimed qualities
                            trim_q = qualities[left-1:(len(seq)-1)-right:1]
                            output_file.write("%s%s\n"% (sign,trim_q))
                            
                        #Create the sequence of trimmed bases
                        match = re.search(trim_seq,seq)
                        if match:
                            match = re.search(trim_seq,seq)
                            count_seq = seq[:match.start():1]
                            count_seq += seq[match.end()::1]
                        else:
                            count_seq = seq
                            
                #Conditional to write the trimmed sequence as non if ringht+left trim are longer than the sequence
                elif right+left > len(seq)-2:
                        count_seq =" "
                        #Conditional for FASTQ format
                        if "fastq" in validation:   
                            sign = input_file.readline()#+
                            qualities = input_file.readline()#qualities
                    
                total_count_seq += count_seq
                total_seq += seq           

    else:#match1 and match2 aren't numbers
        print("Make sure that trim-right and trim-left are numbers")
        exit(1)

    percentage(total_seq, count_n)
    percentage(total_count_seq, trim_n)       
    print("File '%s' has been successfully hard-trimmed ('%s')"%(form1, form2))
    make_summary(0, count_n, trim_n) 
    return output_file

##Validate files
form1, input_file = validation_fast_file(arguments["--input"])
form2, output_file = output(arguments["--output"])

## preprocesing databases
if arguments["--operation"] == "trim":
    trim(input_file, output_file, arguments["--trim-right"], arguments["--trim-left"],form1)
else:
    print("Operation '%s' is not recognized" % arguments["--operation"])
    exit(1)

input_file.close()
output_file.close()

