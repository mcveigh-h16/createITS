# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 08:39:05 2021

@author: mcveigh
"""

#
# processITS designed to validate ITS sequences and create a feature table that can be imported into gbench
# User must specify the input filename and the outputfile name. This version uses fasta as the input format

import pandas as pd
import Bio
import os
import sys
from datetime import datetime
import functools

startTime = datetime.now()
print("Start time is ", startTime) 

inputfile = sys.argv[1]
outputfile = sys.argv[2]
#print(inputfile)
#print(outputfile)

from Bio import SeqIO
sequences = [] 
sequenceLength = []
orgname = []
for seq_record in SeqIO.parse(inputfile, "fasta"):  
    #seq_record.description = seq_record.annotations["organism"]
    #seq_record.description = "contains"
    str_id = seq_record.id      
    sequences.append(seq_record)
    seqLength = '%s %i\n' %  (seq_record.id, len(seq_record))
    sequenceLength.append(seqLength)
    if seq_record.seq.count("NNNNN"):
        print(seq_record.id, "contains internal Ns, this may result in incorrect predictions")
    #orgname = seq_record.annotations["organism"]
    #if "Giardia" in str(orgname):
        #print(seq_record.id, "From Giarda sp. and will require special handling")  
SeqIO.write(sequences, "input.fsa", "fasta")  
seqlen_str = functools.reduce(lambda a,b : a + b, sequenceLength) 
f = open('my.seqlen', 'w')
f.write(seqlen_str)
f.close()
  
#Run CMscan on all sequences
#print('cmscan1')
os.system("cmscan --cpu 16 --mid -T 20 --verbose --tblout tblout.df.txt rrna.cm input.fsa > /dev/null")
#print('cmscan2')
os.system("cmscan --cpu 16 --mid -T 20 --verbose --anytrunc --tblout tblout.at.txt rrna.cm input.fsa > /dev/null")
cmscanTime = datetime.now()
print("CMscan time is ", cmscanTime) 
os.system("cat tblout.df.txt tblout.at.txt > tblout.both.txt")
os.system("perl cmsearch_tblout_deoverlap/cmsearch-deoverlap.pl --maxkeep -s --cmscan tblout.both.txt")
os.system("head -n2 tblout.both.txt > final.tblout")
os.system("cat tblout.both.txt.deoverlapped >> final.tblout")

#Add seq Length to cmscan output
#os.system("esl-seqstat -a ribo-out/ribo-out.ribodbmaker.final.fa | grep ^\\= | awk '{ printf(\"%s %s\\n\", $2, $3); }' > my.seqlen")
os.system("perl tblout-add.pl -t final.tblout 18 my.seqlen 3 > cmscan_final.tblout")

#Parse the final results of CMscan, sort and write fasta files
CMscan_output = (r'cmscan_final.tblout')
CMscan_df = pd.read_csv(CMscan_output,
                        sep='\t',
                        index_col=None,
                        low_memory=False,
                        usecols=[0,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17],
                        header=None,
                        names=["gene", "accession","model", "mdl_from",
                               "mdl_to", "seq_from", "seq_to", "strand",
                               "trunc", "pass", "gc", "bias", "score",
                               "E-value", "Inc", "Length"])
    
CMscan_df['mdl_from']=pd.to_numeric(CMscan_df['mdl_from'])
CMscan_df['mdl_to']=pd.to_numeric(CMscan_df['mdl_to'])
CMscan_df['seq_from']=pd.to_numeric(CMscan_df['seq_from']) 
CMscan_df['seq_to']=pd.to_numeric(CMscan_df['seq_to'])
CMscan_df['Length']=pd.to_numeric(CMscan_df['Length'])  
CMscan_df['score']=pd.to_numeric(CMscan_df['score'])
#CMscan_df = CMscan_df[1:]
#print(CMscan_df.head(10))

#Remove 5S model rows
CMscan_df2 = CMscan_df[CMscan_df['gene'] != "5S_rRNA"]
#print(CMscan_df2.head(20))

#Find sequences that do not pass CMscan tests
fail_test = CMscan_df2[CMscan_df2['Inc'] == "?"]
#print("Sequences that have a ? are\n ", fail_test)

SSU_RNA_df = CMscan_df2[CMscan_df2['gene'] == "SSU_rRNA_eukarya"]
Five_RNA_df = CMscan_df2[CMscan_df2['gene'] == "5_8S_rRNA"]
LSU_RNA_df = CMscan_df2[CMscan_df2['gene'] == "LSU_rRNA_eukarya"]
FiveComplete = Five_RNA_df[Five_RNA_df['trunc'] == "no"]
FiveCompleteStrongHit = FiveComplete[FiveComplete['score'] > 50]
FiveFivePartial = Five_RNA_df[Five_RNA_df['trunc'] == "5'"]
FiveThreePartial = Five_RNA_df[Five_RNA_df['trunc'] == "3'"]
FiveBothPartial = Five_RNA_df[Five_RNA_df['trunc'] == "5'&3'"]
LSUpartial = CMscan_df2.loc[(CMscan_df2['gene'] == "LSU_rRNA_eukarya") & (CMscan_df2['trunc'] != "no")]
LSUcomplete = CMscan_df2.loc[(CMscan_df2['gene'] == "LSU_rRNA_eukarya") & (CMscan_df2['trunc'] == "no")]
SSUpartial = CMscan_df2.loc[(CMscan_df2['gene'] == "SSU_rRNA_eukarya") & (CMscan_df2['trunc'] != "no")]
SSUcomplete = CMscan_df2.loc[(CMscan_df2['gene'] == "SSU_rRNA_eukarya") & (CMscan_df2['trunc'] == "no")]
LSUendmissing = LSUpartial[(LSUpartial['strand'] == "+") & (LSUpartial['seq_to'] != LSUpartial['Length'])]
SSUthreepartial = SSUpartial[SSUpartial['trunc'] == "3'"]
SSUfivefound = SSUpartial[SSUpartial['mdl_from'] == 1]
SSUfivepartial =  SSUpartial[SSUpartial['trunc'] == "5'"]
SSUbothpartial = SSUpartial[SSUpartial['trunc'] == "5'&3'"]
LSUendfound = LSUpartial[(LSUpartial['strand'] == "+") & (LSUpartial['mdl_to'] == 3401)]
LSUFivePartial = LSUpartial[LSUpartial['trunc'] == "5'"]
LSUThreePartial = LSUpartial[LSUpartial['trunc'] == "3'"]
LSUbothpartial = LSUpartial[LSUpartial['trunc'] == "5'&3'"]
five_minus_strand = FiveCompleteStrongHit[FiveCompleteStrongHit['strand'] != "+"]
five_minus_any = Five_RNA_df[Five_RNA_df['strand'] != "+"]
five_plus_any = Five_RNA_df[Five_RNA_df['strand'] == "+"]
LSUpartial_minus_strand = LSUpartial[LSUpartial['strand'] != "+"]
LSUcomplete_minus_strand = LSUcomplete[LSUcomplete['strand'] != "+"]
LSUpartial_plus_strand = LSUpartial[LSUpartial['strand'] == "+"]
LSUcomplete_plus_strand = LSUcomplete[LSUcomplete['strand'] == "+"]
LSUplusframe = [LSUpartial_plus_strand, LSUcomplete_plus_strand]
LSUframes = [LSUpartial_minus_strand, LSUcomplete_minus_strand]
LSUplus = pd.concat(LSUplusframe)
LSUminus = pd.concat(LSUframes)
SSUpartial_minus_strand = SSUpartial[SSUpartial['strand'] != "+"]
SSUcomplete_minus_strand = SSUcomplete[SSUcomplete['strand'] != "+"]
SSUpartial_plus_strand = SSUpartial[SSUpartial['strand'] == "+"]
SSUcomplete_plus_strand = SSUcomplete[SSUcomplete['strand'] == "+"]
SSUplusframe = [SSUpartial_plus_strand, SSUcomplete_plus_strand]
SSUplus = pd.concat(SSUplusframe)
SSUframes = [SSUpartial_minus_strand, SSUcomplete_minus_strand]
SSUminus = pd.concat(SSUframes)
Minusframes = [SSUminus, five_minus_any, LSUminus]
AnyMinus = pd.concat(Minusframes)
Plusframes = [SSUplus, five_plus_any, LSUplus]
AnyPlus = pd.concat(Plusframes)
#print("I found 5.8S sequences on the minus strand\n ", five_minus_strand)
#print("I found LSU sequences on the minus strand\n ", LSUminus)
#print("I found SSU sequences on the minus strand\n ", SSUminus)

#Sequences with extra sequence on the 5' the extends beyond position 1 of the SSU model
SSU_not5_end = SSU_RNA_df[(SSU_RNA_df['seq_from'] != 1) & (SSU_RNA_df['strand'] == "+")]
SSUextra = SSU_not5_end[SSU_not5_end['mdl_from'] == 1]
print("sequences with extra data on the 5' end of SSU\n", SSUextra)
#Sequnces with extra data on the 5' or 3' end beyond the end of LSU
LSU_RNA_df = CMscan_df2[(CMscan_df2['gene'] == "LSU_rRNA_eukarya") & (CMscan_df2['strand'] == "+")]
LSU_not5_end = LSU_RNA_df[(LSU_RNA_df['mdl_from'] != 1) & (LSU_RNA_df['strand'] == "+") & (LSU_RNA_df['seq_from'] != 1)]
print("sequences with extra data on the 5' end of LSU\n", LSU_not5_end)
LSUextra=LSU_RNA_df.loc[(LSU_RNA_df['seq_to'] != LSU_RNA_df['Length']) & (LSU_RNA_df['mdl_to'] == 3401) & (LSU_RNA_df['mdl_from'] == 1)]
print("sequences with extra data on the 3' end of LSU\n", LSUextra)


from Bio import SeqIO
sequences = []  
misassembled = []
removeacc = []
rna_not_found = []
rna_found = []
minus_strand_acc = []

for seq_record in SeqIO.parse("input.fsa", "fasta"): 
    s = seq_record
    SSUstart = []
    SSUend = []
    fivestart = []
    fiveend = []
    LSUstart = []
    LSUend = []
    #print("test1 ", seq_record.id, SSUstart, SSUend, fivestart, fiveend, LSUstart, LSUend)
#check for sequences with no rRNA gene
    if seq_record.id not in SSU_RNA_df['accession'].tolist():
        if seq_record.id not in Five_RNA_df['accession'].tolist():
            if seq_record.id not in LSU_RNA_df['accession'].tolist():
                print(seq_record.id, "No rRNA gene was found")    
                rna_not_found.append(s)
                removeacc.append(seq_record.id)
                SeqIO.write(rna_not_found, "rna_not_found_seqs", "fasta") 
                #print(rna_not_found)
    
#Check for Mixed Strand, Noncontiguous and Misassembled Sequences    
    if seq_record.id in AnyPlus['accession'].tolist():
        if seq_record.id in AnyMinus['accession'].tolist():
            print(seq_record.id, "Mixed strand sequence")
            misassembled.append(s)
            removeacc.append(seq_record.id)
        if seq_record.id in SSU_RNA_df['accession'].tolist():
            start_df = SSU_RNA_df[(SSU_RNA_df['accession'] == seq_record.id) & (SSU_RNA_df['gene'] == 'SSU_rRNA_eukarya')]
            SSUstart = int(start_df['seq_from'].iloc[0])
            SSUend = int(start_df['seq_to'].iloc[0])
            #print(seq_record.id, " SSU ", SSUstart, SSUend)
        if seq_record.id in Five_RNA_df['accession'].tolist():
            startfive_df = Five_RNA_df[(Five_RNA_df['accession'] == seq_record.id) & (Five_RNA_df['gene'] == '5_8S_rRNA')]
            fivestart = int(startfive_df['seq_from'].iloc[0])
            fiveend = int(startfive_df['seq_to'].iloc[0])
            #print(seq_record.id, " FIVE ", fivestart, fiveend)
        if seq_record.id in LSU_RNA_df['accession'].tolist(): 
            startLSU_df = LSU_RNA_df[(LSU_RNA_df['accession'] == seq_record.id) & (LSU_RNA_df['gene'] == 'LSU_rRNA_eukarya')]
            LSUstart = int(startLSU_df['seq_from'].iloc[0])
            LSUend = int(startLSU_df['seq_to'].iloc[0])
            #print(seq_record.id, " LSU ", LSUstart, LSUend)
    elif seq_record.id in AnyMinus['accession'].tolist():
        if seq_record.id in AnyPlus['accession'].tolist():
            print(seq_record.id, "Mixed strand sequence")
            misassembled.append(s)
            removeacc.append(seq_record.id)
        if seq_record.id in SSU_RNA_df['accession'].tolist():
            start_df = SSU_RNA_df[(SSU_RNA_df['accession'] == seq_record.id) & (SSU_RNA_df['gene'] == 'SSU_rRNA_eukarya')]
            SSUstart = int(start_df['seq_to'].iloc[0])
            SSUend = int(start_df['seq_from'].iloc[0])
        if seq_record.id in Five_RNA_df['accession'].tolist():
            startfive_df = Five_RNA_df[(Five_RNA_df['accession'] == seq_record.id) & (Five_RNA_df['gene'] == '5_8S_rRNA')]
            fivestart = startfive_df['seq_to'].iloc[0]
            fiveend = startfive_df['seq_from'].iloc[0]
        if seq_record.id in LSU_RNA_df['accession'].tolist():  
            startLSU_df = LSU_RNA_df[(LSU_RNA_df['accession'] == seq_record.id) & (LSU_RNA_df['gene'] == 'LSU_rRNA_eukarya')]
            LSUstart = startLSU_df['seq_to'].iloc[0]
            LSUend = startLSU_df['seq_from'].iloc[0]
    #print("test2 ", seq_record.id, SSUstart, SSUend, fivestart, fiveend, LSUstart, LSUend)
    
    if seq_record.id in SSU_RNA_df['accession'].tolist():
       if seq_record.id in Five_RNA_df['accession'].tolist():
            if SSUend > fivestart:
               print(seq_record.id, " ssu spans ", SSUstart, SSUend)
               print(seq_record.id, " five spans ", fivestart, fiveend)
               print(seq_record.id, "Misassembled sequence 1")
               misassembled.append(s)
               removeacc.append(seq_record.id)
    if seq_record.id in Five_RNA_df['accession'].tolist():
        if seq_record.id in LSU_RNA_df['accession'].tolist():
            if fiveend > LSUstart:
                print(seq_record.id, "Misassembled sequence 2")
                misassembled.append(s)
                removeacc.append(seq_record.id)

#noncontig seq check
    if seq_record.id in SSU_RNA_df['accession'].tolist():
        if seq_record.id in LSU_RNA_df['accession'].tolist():
            if seq_record.id not in FiveCompleteStrongHit['accession'].tolist():
                print(seq_record.id, "Noncontiguous sequence")
                misassembled.append(s)
                removeacc.append(seq_record.id)      
    SeqIO.write(misassembled, "misassembled_seqs", "fasta")  
    #reject_df = pd.DataFrame([removeacc])
    #reject_df.columns =['accession']
    #print("rejectdf accession", reject_df)
    if seq_record.id not in removeacc:
        rna_found.append(s)
        SeqIO.write(rna_found, "rna_found.seqs", "fasta")  

line = []
for seq_record in SeqIO.parse("rna_found.seqs", "fasta"):
    start = []
    stop = []
    #ITS1start = []
    #ITS1end = []
    #fivestart = []
    #fiveend = []
    #ITS2start = []
    #ITS2end = []
    LSUstart = []
    LSUend = []
    length = []
    mdlfrom = []
    mdlto = []
    SSUstart = []
    SSUend = []
    seq_record.description = "contains"
    
#parse CMscan output file and write five column feature table
#plus strand sequences
    if seq_record.id in AnyPlus['accession'].tolist(): 
        if seq_record.id in SSU_RNA_df['accession'].tolist(): 
            if seq_record.id in SSUcomplete['accession'].tolist(): 
                start_df = SSUcomplete[(SSUcomplete['accession'] == seq_record.id) & (SSUcomplete['gene'] == 'SSU_rRNA_eukarya')]
                start = start_df['seq_from'].iloc[0]
                stop = start_df['seq_to'].iloc[0]
                length = start_df['Length'].iloc[0]
                if (start == 1) & (stop == length):
                    seq_record.description = "small subunit ribosomal RNA" 
                if seq_record.id in Five_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + " small subunit ribosomal RNA, internal transcribed spacer 1"
                elif seq_record.id not in Five_RNA_df['accession'].tolist():
                    if stop != length:
                        seq_record.description = seq_record.description + " small subunit ribosomal RNA and internal transcribed spacer 1"
                        stop = ">" + str(length)       
            elif seq_record.id in SSUpartial['accession'].tolist():  
                start_df = SSUpartial[(SSUpartial['accession'] == seq_record.id) & (SSUpartial['gene'] == 'SSU_rRNA_eukarya')]
                start = start_df['seq_from'].iloc[0]
                stop = start_df['seq_to'].iloc[0]
                length = start_df['Length'].iloc[0]
                mdlfrom = start_df['mdl_from'].iloc[0]
                mdlto = start_df['mdl_to'].iloc[0]
                if seq_record.id in SSUfivefound['accession'].tolist():
                    if (start == 1) & (stop == length):
                        seq_record.description = "small subunit ribosomal RNA" 
                        stop = ">" + str(stop)
                    elif (start > 1):
                        seq_record.description = "small subunit ribosomal RNA"
                if seq_record.id in SSUfivepartial['accession'].tolist():
                    start = "<" + str(start)
                    if stop == length:
                        seq_record.description = "small subunit ribosomal RNA"                
                    else:
                        stop = ">" + str(length)  
                        seq_record.description = "contains small subunit ribosomal RNA and internal transcribed spacer 1"
                if seq_record.id in Five_RNA_df['accession'].tolist():
                    seq_record.description = "contains small subunit ribosomal RNA, internal transcribed spacer 1"
                if seq_record.id in SSUbothpartial['accession'].tolist():
                    if (start == 1) & (stop == length):
                        seq_record.description = "small subunit ribosomal RNA"
                        start = "<" + str(start)
                        stop = ">" + str(stop)
                    else:
                        seq_record.description = "small subunit ribosomal RNA"
                        start = "<" + str(start)
                        stop = ">" + str(stop)
                        print(seq_record.id, " is incomplete on both ends but the RNA prediction does not extend to the end")
        if seq_record.id in Five_RNA_df['accession'].tolist(): 
            if seq_record.id in FiveCompleteStrongHit['accession'].tolist():  
                if seq_record.id not in SSU_RNA_df['accession'].tolist():
                    start_df = FiveCompleteStrongHit[(FiveCompleteStrongHit['accession'] == seq_record.id) & (FiveCompleteStrongHit['gene'] == '5_8S_rRNA')]
                    start = start_df['seq_from'].iloc[0]
                    if start > 1:
                        seq_record.description = seq_record.description + " internal transcribed spacer 1, 5.8S ribosomal RNA"
                        start = "<1"
                if seq_record.id in SSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", 5.8S ribosomal RNA" 
                if seq_record.id not in LSU_RNA_df['accession'].tolist():
                    stop_df = FiveCompleteStrongHit[(FiveCompleteStrongHit['accession'] == seq_record.id) & (FiveCompleteStrongHit['gene'] == '5_8S_rRNA')]
                    length_df = FiveCompleteStrongHit[(FiveCompleteStrongHit['accession'] == seq_record.id) & (FiveCompleteStrongHit['gene'] == '5_8S_rRNA')]
                    stop = stop_df['seq_to'].iloc[0]
                    length = length_df['Length'].iloc[0]
                    if stop < length:
                        seq_record.description = seq_record.description + ", and internal transcribed spacer 2"
                        stop = ">" + str(length)
                    if stop == length:
                        seq_record.description = "5.8S ribosomal RNA"
            elif seq_record.id in FiveFivePartial['accession'].tolist():
                start = "<1"
                stop_df = FiveFivePartial[(FiveFivePartial['accession'] == seq_record.id) & (FiveFivePartial['gene'] == '5_8S_rRNA')]
                stop = stop_df['Length'].iloc[0]
                if seq_record.id in LSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + " 5.8S ribosomal RNA, internal transcribed spacer 2"
                elif seq_record.id not in LSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + " 5.8S ribosomal RNA and internal transcribed spacer 2"
                    stop = ">" + str(stop)
            elif seq_record.id in FiveThreePartial['accession'].tolist():
                if seq_record.id not in SSU_RNA_df['accession'].tolist():
                    start_df = FiveThreePartial[(FiveThreePartial['accession'] == seq_record.id) & (FiveThreePartial['gene'] == '5_8S_rRNA')]
                    start = start_df['seq_from'].iloc[0]
                    if start == 1:
                        seq_record.description = "5.8S ribosomal RNA"
                    else:
                        start = "<1"
                        seq_record.description = seq_record.description + " internal transcribed spacer 1 and 5.8S ribosomal RNA"
                elif seq_record.id in SSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", 5.8S ribosomal RNA"
                stop_df = FiveThreePartial[(FiveThreePartial['accession'] == seq_record.id) & (FiveThreePartial['gene'] == '5_8S_rRNA')]
                stop = stop_df['seq_to'].iloc[0]
                stop = ">" + str(stop)
            elif seq_record.id in FiveBothPartial['accession'].tolist():
                seq_record.description = "5.8S ribosomal RNA"
                start = "<1" 
                stop_df = FiveBothPartial[(FiveBothPartial['accession'] == seq_record.id) & (FiveBothPartial['gene'] == '5_8S_rRNA')]
                stop = stop_df['seq_to'].iloc[0] 
                stop = ">" + str(stop) 
            else:
                if seq_record.id not in SSU_RNA_df['accession'].tolist():
                    start_df = Five_RNA_df[(Five_RNA_df['accession'] == seq_record.id) & (Five_RNA_df['gene'] == '5_8S_rRNA')]
                    start = start_df['seq_from'].iloc[0]
                    length_df = Five_RNA_df[(Five_RNA_df['accession'] == seq_record.id) & (Five_RNA_df['gene'] == '5_8S_rRNA')]
                    length = length_df['Length'].iloc[0]
                    if start > 1:
                        seq_record.description = seq_record.description + " internal transcribed spacer 1, 5.8S ribosomal RNA"
                        start = "<1"
                if seq_record.id in SSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", 5.8S ribosomal RNA" 
                stop_df = Five_RNA_df[(Five_RNA_df['accession'] == seq_record.id) & (Five_RNA_df['gene'] == '5_8S_rRNA')]
                stop = stop_df['seq_to'].iloc[0]
                #if seq_record.id in LSU_RNA_df['accession'].tolist():
                    #seq_record.description = seq_record.description + ", internal transcribed spacer 2"      
                if seq_record.id not in LSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", and internal transcribed spacer 2"
                    stop = ">" + str(length)
                print(seq_record.id, "Low scoring complete 5.8S gene found this is likely a bad sequence")
        if seq_record.id in LSU_RNA_df['accession'].tolist():  
            if seq_record.id in LSUcomplete['accession'].tolist():
                stop_df = LSUcomplete[(LSUcomplete['accession'] == seq_record.id) & (LSUcomplete['gene'] == 'LSU_rRNA_eukarya')]
                stop = stop_df['seq_to'].iloc[0] 
                LSUstart = stop_df['seq_from'].iloc[0]
                if seq_record.id not in Five_RNA_df['accession'].tolist(): 
                    if LSUstart == 1:
                        start = "1"
                        seq_record.description = "large subunit ribosomal RNA"
                    else:
                        seq_record.description = seq_record.description + " internal transcribed spacer 2 and large subunit ribosomal RNA"
                        start = "<1"
                elif seq_record.id in Five_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", internal transcribed spacer 2 and large subunit ribosomal RNA"
            elif seq_record.id in LSUpartial['accession'].tolist():
                stop_df = LSUpartial[(LSUpartial['accession'] == seq_record.id) & (LSUpartial['gene'] == 'LSU_rRNA_eukarya')]
                LSUstart = stop_df['seq_from'].iloc[0]
                stop = stop_df['seq_to'].iloc[0]
                mdlfrom = stop_df['mdl_from'].iloc[0]
                stop = ">" + str(stop)
                if seq_record.id not in Five_RNA_df['accession'].tolist():
                    if seq_record.id in LSUbothpartial['accession'].tolist():
                        start = "<1"
                        seq_record.description = "large subunit ribosomal RNA"
                    if seq_record.id in LSUThreePartial['accession'].tolist():               
                        if LSUstart == "1":
                            seq_record.description = "large subunit ribosomal RNA"
                        else: 
                            seq_record.description = seq_record.description + " internal transcribed spacer 2 and large subunit ribosomal RNA"
                            start = "<1"
                    elif seq_record.id in LSUFivePartial['accession'].tolist():
                        seq_record.description = "large subunit ribosomal RNA"
                        start = "<1"
                elif seq_record.id in Five_RNA_df['accession'].tolist():
                    if seq_record.id in FiveFivePartial['accession'].tolist():
                        seq_record.description = seq_record.description + ", and large subunit ribosomal RNA"
                    else:
                        seq_record.description = seq_record.description + ", internal transcribed spacer 2 and large subunit ribosomal RNA"
        sequences.append(s)
#Minus strand sequences        
    if seq_record.id in AnyMinus['accession'].tolist():
        minus_strand_acc.append(seq_record.id)
        if seq_record.id in SSU_RNA_df['accession'].tolist(): 
            if seq_record.id in SSUcomplete['accession'].tolist():  
                seq_record.description = seq_record.description + " small subunit ribosomal RNA" 
                start_df = SSUcomplete[(SSUcomplete['accession'] == seq_record.id) & (SSUcomplete['gene'] == 'SSU_rRNA_eukarya')]
                start = start_df['seq_from'].iloc[0]           
            #if seq_record.id in SSUpartial['accession'].tolist():  
            else:
                seq_record.description = seq_record.description + " small subunit ribosomal RNA"
                start_df = SSUpartial[(SSUpartial['accession'] == seq_record.id) & (SSUpartial['gene'] == 'SSU_rRNA_eukarya')]
                start = start_df['seq_from'].iloc[0]
                start = "<" + str(start)
            if seq_record.id in Five_RNA_df['accession'].tolist():
                seq_record.description = seq_record.description + ", internal transcribed spacer 1"
            elif seq_record.id not in Five_RNA_df['accession'].tolist():
                stop_df = SSU_RNA_df[(SSU_RNA_df['accession'] == seq_record.id) & (SSU_RNA_df['gene'] == 'SSU_rRNA_eukarya')]
                length_df = SSU_RNA_df[(SSU_RNA_df['accession'] == seq_record.id) & (SSU_RNA_df['gene'] == 'SSU_rRNA_eukarya')]
                stop = stop_df['seq_to'].iloc[0]
                length = length_df['Length'].iloc[0]
                if stop == 1:
                    if seq_record.id in SSUcomplete['accession'].tolist():
                        seq_record.description = "small subunit ribosomal RNA"
                    elif seq_record.id in SSUpartial['accession'].tolist():
                        stop = ">" + str(stop)
                        seq_record.description = "small subunit ribosomal RNA"
                else:
                    seq_record.description = seq_record.description + "and internal transcribed spacer 1"
                    stop = "<1"
        if seq_record.id in Five_RNA_df['accession'].tolist(): 
            if seq_record.id in FiveCompleteStrongHit['accession'].tolist():  
                if seq_record.id not in SSU_RNA_df['accession'].tolist():
                    start_df = FiveCompleteStrongHit[(FiveCompleteStrongHit['accession'] == seq_record.id) & (FiveCompleteStrongHit['gene'] == '5_8S_rRNA')]
                    length = start_df['Length'].iloc[0]
                    start = start_df['seq_from'].iloc[0]
                    if start < length:
                        seq_record.description = seq_record.description + " internal transcribed spacer 1, 5.8S ribosomal RNA"
                        start = "<" + str(length)
                if seq_record.id in SSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", 5.8S ribosomal RNA" 
                if seq_record.id not in LSUminus['accession'].tolist():
                    stop_df = FiveCompleteStrongHit[(FiveCompleteStrongHit['accession'] == seq_record.id) & (FiveCompleteStrongHit['gene'] == '5_8S_rRNA')]
                    length_df = FiveCompleteStrongHit[(FiveCompleteStrongHit['accession'] == seq_record.id) & (FiveCompleteStrongHit['gene'] == '5_8S_rRNA')]
                    stop = stop_df['seq_to'].iloc[0]
                    length = length_df['Length'].iloc[0]
                    if stop > 1:
                        seq_record.description = seq_record.description + ", and internal transcribed spacer 2"
                        stop = ">1"
                    if stop == 1:
                        seq_record.description = "5.8S ribosomal RNA"
            elif seq_record.id in FiveFivePartial['accession'].tolist():
                length_df = FiveFivePartial[(FiveFivePartial['accession'] == seq_record.id) & (FiveFivePartial['gene'] == '5_8S_rRNA')]
                length = length_df['Length'].iloc[0]
                start = ">" + str(length)
                stop_df = FiveFivePartial[(FiveFivePartial['accession'] == seq_record.id) & (FiveFivePartial['gene'] == '5_8S_rRNA')]
                stop = stop_df['seq_to'].iloc[0]
                if seq_record.id in LSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + " 5.8S ribosomal RNA, internal transcribed spacer 2"
                elif seq_record.id not in LSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + " 5.8S ribosomal RNA and internal transcribed spacer 2"
                    stop = "<" + str(stop)
            elif seq_record.id in FiveThreePartial['accession'].tolist():
                if seq_record.id not in SSU_RNA_df['accession'].tolist():
                    start_df = FiveThreePartial[(FiveThreePartial['accession'] == seq_record.id) & (FiveThreePartial['gene'] == '5_8S_rRNA')]
                    start = start_df['seq_from'].iloc[0]
                    length_df = FiveThreePartial[(FiveThreePartial['accession'] == seq_record.id) & (FiveThreePartial['gene'] == '5_8S_rRNA')]
                    length = length_df['Length'].iloc[0]
                    if start == length:
                        seq_record.description = "5.8S ribosomal RNA"
                    else: 
                        start = ">" + str(length)
                        seq_record.description = seq_record.description + " internal transcribed spacer 1 and 5.8S ribosomal RNA"
                elif seq_record.id in SSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", 5.8S ribosomal RNA"
                stop_df = FiveThreePartial[(FiveThreePartial['accession'] == seq_record.id) & (FiveThreePartial['gene'] == '5_8S_rRNA')]
                stop = stop_df['seq_to'].iloc[0]
                stop = ">" + str(stop)
            elif seq_record.id in FiveBothPartial['accession'].tolist():
                seq_record.description = "5.8S ribosomal RNA"
                start = "<" + str(length) 
                stop_df = FiveBothPartial[(FiveBothPartial['accession'] == seq_record.id) & (FiveBothPartial['gene'] == '5_8S_rRNA')]
                stop = stop_df['seq_to'].iloc[0] 
                stop = ">" + str(stop)    
            else:
                if seq_record.id not in SSU_RNA_df['accession'].tolist():
                    start_df = Five_RNA_df[(Five_RNA_df['accession'] == seq_record.id) & (Five_RNA_df['gene'] == '5_8S_rRNA')]
                    start = start_df['seq_from'].iloc[0]
                    length_df = Five_RNA_df[(Five_RNA_df['accession'] == seq_record.id) & (Five_RNA_df['gene'] == '5_8S_rRNA')]
                    length = length_df['Length'].iloc[0]
                    if start < length:
                        seq_record.description = seq_record.description + " internal transcribed spacer 1, 5.8S ribosomal RNA"
                        start = "<" + str(start)
                if seq_record.id in SSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", 5.8S ribosomal RNA" 
                stop_df = Five_RNA_df[(Five_RNA_df['accession'] == seq_record.id) & (Five_RNA_df['gene'] == '5_8S_rRNA')]
                stop = stop_df['seq_to'].iloc[0]
                if seq_record.id in LSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", internal transcribed spacer 2"      
                elif seq_record.id not in LSU_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", and internal transcribed spacer 2"
                    stop = "<1" 
                print(seq_record.id, "Low scoring complete 5.8S gene found this is likely a bad sequence")    
        if seq_record.id in LSUminus['accession'].tolist():             
            if seq_record.id in LSUcomplete['accession'].tolist():
                stop_df = LSUcomplete[(LSUcomplete['accession'] == seq_record.id) & (LSUcomplete['gene'] == 'LSU_rRNA_eukarya')]
                stop = stop_df['seq_to'].iloc[0] 
                LSUstart = stop_df['seq_from'].iloc[0]
                length = stop_df['length'].iloc[0]
                if seq_record.id not in Five_RNA_df['accession'].tolist(): 
                    if LSUstart == length:
                        start = LSUstart
                        seq_record.description = "large subunit ribosomal RNA"
                    else:
                        seq_record.description = seq_record.description + " internal transcribed spacer 2 and large subunit ribosomal RNA"
                        start = "<" + str(LSUstart)
                elif seq_record.id in Five_RNA_df['accession'].tolist():
                    seq_record.description = seq_record.description + ", internal transcribed spacer 2 and large subunit ribosomal RNA"
            elif seq_record.id in LSUpartial['accession'].tolist():
                stop_df = LSUpartial[(LSUpartial['accession'] == seq_record.id) & (LSUpartial['gene'] == 'LSU_rRNA_eukarya')]
                LSUstart = stop_df['seq_from'].iloc[0]
                stop = stop_df['seq_to'].iloc[0]
                mdlfrom = stop_df['mdl_from'].iloc[0]
                length = stop_df['Length'].iloc[0]
                stop = ">" + str(stop)
                if seq_record.id not in Five_RNA_df['accession'].tolist():
                    if LSUstart == length:
                        seq_record.description = "large subunit ribosomal RNA"
                        if mdlfrom == 1:
                            start = LSUstart
                        else:
                            start = "<" + str(LSUstart)
                    else:
                        seq_record.description = seq_record.description + " internal transcribed spacer 2 and large subunit ribosomal RNA"
                        start = "<" + str(length)
                elif seq_record.id in Five_RNA_df['accession'].tolist():
                    if seq_record.id in FiveFivePartial['accession'].tolist():
                        seq_record.description = seq_record.description + ", and large subunit ribosomal RNA"
                    else:
                        seq_record.description = seq_record.description + ", internal transcribed spacer 2 and large subunit ribosomal RNA"
        sequences.append(s)
        
    #print(seq_record.description, start, stop)
    if seq_record.description == str("small subunit ribosomal RNA"):
        line.append(">Feature " + seq_record.id + "\n" + str(start) + "\t" + str(stop) + "\trRNA\n\t\t\t\tproduct\t" + seq_record.description + "\n")
    if seq_record.description == str("large subunit ribosomal RNA"):
        line.append(">Feature " + seq_record.id + "\n" + str(start) + "\t" + str(stop) + "\trRNA\n\t\t\t\tproduct\t" + seq_record.description + "\n")
    if seq_record.description == str("5.8S ribosomal RNA"):
        line.append(">Feature " + seq_record.id + "\n" + str(start) + "\t" + str(stop) + "\trRNA\n\t\t\t\tproduct\t" + seq_record.description + "\n")
    if "internal" in str(seq_record.description):
        line.append(">Feature " + seq_record.id + "\n" + str(start) + "\t" + str(stop) + "\tmisc_RNA\n\t\t\t\tnote\t" + seq_record.description + "\n")
        
#print("sequences on the minus strand are: \n", minus_strand_acc) 
stdout_fileno = sys.stdout 
for ip in minus_strand_acc:
    stdout_fileno.write(ip + ' is on the minus strand\n')      
#print(removeacc, "removed accessions")
#SeqIO.write(sequences, outputfile, "fasta")     
#print(line)
file = open(outputfile, "w")
file.writelines(line)
file.close()
scriptTime = datetime.now()
print("script time is ", scriptTime)             
            