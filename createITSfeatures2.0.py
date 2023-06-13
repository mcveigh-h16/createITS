# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:50:11 2023

@author: mcveigh

    processITS designed to validate ITS sequences and create a feature table that can be imported into gbench
    User must specify the input filename and the outputfile name. This version uses fasta as the input format.
    Minus sequences are reverse complemented and CMscan is rerun. This puts all sequences on the plus strand
    prior to analysis simplifying the code. 
"""

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
#inputfile = (r'test17.fsa')
#outputfile = (r'test17.out')

from Bio import SeqIO
sequences = [] 
sequenceLength = []
orgname = []
for seq_record in SeqIO.parse(inputfile, "fasta"):  
    str_id = seq_record.id      
    sequences.append(seq_record)
    seqLength = '%s %i\n' %  (seq_record.id, len(seq_record))
    sequenceLength.append(seqLength)
    if seq_record.seq.count("NNNNN"):
        print(seq_record.id, "contains internal Ns, this may result in incorrect predictions")
SeqIO.write(sequences, "input.fsa", "fasta")  
seqlen_str = functools.reduce(lambda a,b : a + b, sequenceLength) 
f = open('my.seqlen', 'w')
f.write(seqlen_str)
f.close()

"""  
Run CMscan on all sequences
"""

os.system("cmscan --cpu 16 --mid -T 20 --verbose --tblout tblout.df.txt rrna.cm input.fsa > /dev/null")
os.system("cmscan --cpu 16 --mid -T 20 --verbose --anytrunc --tblout tblout.at.txt rrna.cm input.fsa > /dev/null")
cmscanTime = datetime.now()
print("CMscan time is ", cmscanTime) 
os.system("cat tblout.df.txt tblout.at.txt > tblout.both.txt")
os.system("perl cmsearch_tblout_deoverlap/cmsearch-deoverlap.pl --maxkeep -s --cmscan tblout.both.txt")
os.system("head -n2 tblout.both.txt > final.tblout")
os.system("cat tblout.both.txt.deoverlapped >> final.tblout")

"""
Add lengths to cmscan out files
"""
os.system("perl tblout-add.pl -t final.tblout 18 my.seqlen 3 > cmscan_final.tblout")

"""
Parse the results of CMscan. If any rRNA gene is found on minus strand these
are reverse complemented and rerun through CMscan. 
"""
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

#FIND ALL POSSIBLE MINUS STRAND GENES
FiveCompleteMinus = CMscan_df.loc[(CMscan_df['gene'] == "5_8S_rRNA") & (CMscan_df['trunc'] == "no") & (CMscan_df['score'] > 50) & (CMscan_df['strand'] != "+")]
FiveMinus = CMscan_df.loc[(CMscan_df['gene'] == "5_8S_rRNA") & (CMscan_df['strand'] != "+")]
SSUminus = CMscan_df.loc[(CMscan_df['gene'] == "SSU_rRNA_eukarya") & (CMscan_df['strand'] != "+")]
LSUminus = CMscan_df.loc[(CMscan_df['gene'] == "LSU_rRNA_eukarya") & (CMscan_df['strand'] != "+")]
#print('5.8S on minus strand \n', FiveCompleteMinus)
#print('SSU on minus strand \n', SSUminus)
#print('LSU on minus strand \n', LSUminus)

reverse_seq = []
reverse_acc = []
forward_seq = []

for seq_record in SeqIO.parse("input.fsa", "fasta"): 
    s = seq_record
    if seq_record.id in FiveCompleteMinus['accession'].tolist():
        print("I reverse complemented a full length 5.8S gene", seq_record.id)
        s.seq = s.seq.reverse_complement()
        reverse_seq.append(s)
        reverse_acc.append(s.id)
    elif seq_record.id in FiveMinus['accession'].tolist():
        print("I reverse complemented a partial 5.8S gene", seq_record.id)
        s.seq = s.seq.reverse_complement()
        reverse_seq.append(s)
        reverse_acc.append(s.id)
    elif seq_record.id in SSUminus['accession'].tolist():
        print("I reverse complemented a SSU gene", seq_record.id)
        s.seq = s.seq.reverse_complement()
        reverse_seq.append(s)
        reverse_acc.append(s.id)
    elif seq_record.id in LSUminus['accession'].tolist():
        print("I reverse complemented a LSU gene", seq_record.id)
        s.seq = s.seq.reverse_complement()
        reverse_seq.append(s)
        reverse_acc.append(s.id)
    else:
        forward_seq.append(s)
    
SeqIO.write(reverse_seq, "reverse_seqs", "fasta")
SeqIO.write(forward_seq, "forward_seqs", "fasta")
os.system("cat forward_seqs reverse_seqs > plus_strand_seqs")
print('reversed accesssions are: ', reverse_acc)

def repeat_cmscan():
    """
    repeat cmscan for only sequences that have been reverse complemented and
    combines output original. 
    """  

    CMscan_df.drop(CMscan_df[CMscan_df.accession.isin(reverse_acc)].index.tolist())
    Plus_df = CMscan_df[~CMscan_df.accession.isin(reverse_acc)]
    

    """
    Run CMscan on the reverse complemented sequences
    """
    print('cmscan3')
    os.system("cmscan --cpu 16 --mid -T 20 --verbose --tblout repeat.df.txt rrna.cm reverse_seqs > /dev/null")
    print('cmscan4')
    os.system("cmscan --cpu 16 --mid -T 20 --verbose --anytrunc --tblout repeat.at.txt rrna.cm reverse_seqs > /dev/null")
    repeatcmscanTime = datetime.now()
    print("CMscan time is ", repeatcmscanTime) 
    os.system("cat repeat.df.txt repeat.at.txt > repeat.both.txt")
    os.system("perl cmsearch_tblout_deoverlap/cmsearch-deoverlap.pl --maxkeep -s --cmscan repeat.both.txt")
    os.system("head -n2 repeat.both.txt > repeat.tblout")
    os.system("cat repeat.both.txt.deoverlapped >> repeat.tblout")
    os.system("perl tblout-add.pl -t repeat.tblout 18 my.seqlen 3 > repeat_cmscan_final.tblout")

    """
    Collect reverse complemented CMscan data and merging this with the data for
    sequences on the plus strand. All seqs should now be plus strand 5.8S rRNA
    """
    repeat_CMscan_output = (r'repeat_cmscan_final.tblout')
    repeat_CMscan_df = pd.read_csv(repeat_CMscan_output,
                        sep='\t',
                        index_col=None,
                        low_memory=False,
                        usecols=[0,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17],
                        header=None,
                        names=["gene", "accession","model", "mdl_from",
                               "mdl_to", "seq_from", "seq_to", "strand",
                               "trunc", "pass", "gc", "bias", "score",
                               "E-value", "Inc", "Length"])

    combined = [Plus_df, repeat_CMscan_df]
    combined_df = pd.concat(combined)
    print(combined_df)
    print('i finished the function')
    return combined_df
        
if len(reverse_acc) != 0:
    #print('Minus strand sequences are found')
    print('reversed accesssions are: ', reverse_acc)
    repeat_cmscan()
    combined_df = repeat_cmscan()
else:
    combined_df = CMscan_df


"""
Remove 5S model rows from data frame
Find and report sequences that have truncated models
Find sequences that do not pass CMscan tests
"""
combined_df = combined_df[combined_df['gene'] != "5S_rRNA"]
#print(combined_df)

truncated=combined_df.loc[(combined_df['gene'] != "5_8S_rRNA") & (combined_df['trunc'] == "5'&3'")]
print("I found truncated models suggesting the presence of an intron\n ", truncated)

fail_test = combined_df[combined_df['Inc'] == "?"]
print("Sequences that have a ? are\n ", fail_test)


"""
Build Dataframes
"""

Fivecomplete = combined_df.loc[(combined_df['gene'] == "5_8S_rRNA") & (combined_df['trunc'] == "no") & (combined_df['score'] > 50)]
Fivepartial = combined_df.loc[(combined_df['gene'] == "5_8S_rRNA") & (combined_df['trunc'] != "no")]
SSU_RNA_df = combined_df.loc[(combined_df['gene'] == "SSU_rRNA_eukarya")]
LSU_RNA_df = combined_df.loc[(combined_df['gene'] == "LSU_rRNA_eukarya")]  
Five_RNA_df = combined_df.loc[(combined_df['gene'] == "5_8S_rRNA") & (combined_df['strand'] == "+")]
LSUpartial = combined_df.loc[(combined_df['gene'] == "LSU_rRNA_eukarya") & (combined_df['trunc'] != "no")]
LSUcomplete = combined_df.loc[(combined_df['gene'] == "LSU_rRNA_eukarya") & (combined_df['trunc'] == "no")]
SSUpartial = combined_df.loc[(combined_df['gene'] == "SSU_rRNA_eukarya") & (combined_df['trunc'] != "no")]
SSUcomplete = combined_df.loc[(combined_df['gene'] == "SSU_rRNA_eukarya") & (combined_df['trunc'] == "no")]                         
#LSUendmissing = LSUpartial[(LSUpartial['strand'] == "+") & (LSUpartial['seq_to'] != LSUpartial['Length'])]
#SSUendmissing = SSUpartial[(SSUpartial['strand'] == "+") & (SSUpartial['seq_from'] > 1)]
#SSUendfound = SSUpartial[(SSUpartial['strand'] == "+") & (SSUpartial['mdl_from'] == 1)]
#LSUendfound = LSUpartial[(LSUpartial['strand'] == "+") & (LSUpartial['mdl_to'] == 3401)]
SSUminus = combined_df.loc[(combined_df['gene'] == "SSU_rRNA_eukarya") & (combined_df['strand'] != "+")]
LSUminus = combined_df.loc[(combined_df['gene'] == "LSU_rRNA_eukarya") & (combined_df['strand'] != "+")]
SSUextra = SSUcomplete.loc[(SSUcomplete['seq_from'] != 1) & SSUcomplete['mdl_from'] == 1]
LSUextra=LSUcomplete.loc[(LSUcomplete['seq_to'] != LSUcomplete['Length']) & (LSUcomplete['mdl_to'] == 3401) & (LSUcomplete['mdl_from'] == 1)]
print("sequences with extra data on the 5' end \n", SSUextra)
print("sequences with extra data on the 3' end \n", LSUextra)

misassembled = []
removeacc = []
rna_not_found = []
rna_found = []
line = []

"""
Iterate through the sequences and evaluate each
Check for sequences with errors
"""
for seq_record in SeqIO.parse("plus_strand_seqs", "fasta"): 
    s = seq_record
    SSUstart = []
    SSUend = []
    fivestart = []
    fiveend = []
    LSUstart = []
    LSUend = []
    start = []
    stop = []
    length = []
    mdlfrom = []
    mdlto = []
    seq_record.description = "contains"

    """
    check for sequences with no rRNA gene
    """
    if seq_record.id not in SSU_RNA_df['accession'].tolist():
        if seq_record.id not in Five_RNA_df['accession'].tolist():
            if seq_record.id not in LSU_RNA_df['accession'].tolist():
                print(seq_record.id, "No rRNA gene was found")    
                rna_not_found.append(s)
                removeacc.append(seq_record.id)
                SeqIO.write(rna_not_found, "rna_not_found_seqs", "fasta") 
                #print(rna_not_found)
    
    """            
    Check for Mixed Strand, Noncontiguous and Misassembled Sequences 
    """
    if seq_record.id in Fivecomplete['accession'].tolist():
        if seq_record.id in SSUminus['accession'].tolist():
            print(seq_record.id, "Mixed strand sequence 1")
            misassembled.append(s)
            removeacc.append(seq_record.id)
        if seq_record.id in LSUminus['accession'].tolist():
            print(seq_record.id, "Mixed strand sequence 2")
            misassembled.append(s)
            removeacc.append(seq_record.id)
        if seq_record.id in SSU_RNA_df['accession'].tolist():
            start_df = SSU_RNA_df[(SSU_RNA_df['accession'] == seq_record.id) & (SSU_RNA_df['gene'] == 'SSU_rRNA_eukarya')]
            SSUstart = int(start_df['seq_from'].iloc[0])
            SSUend = int(start_df['seq_to'].iloc[0])
            #print(seq_record.id, " SSU ", SSUstart, SSUend)
        if seq_record.id in Fivecomplete['accession'].tolist():
            startfive_df = Fivecomplete[(Fivecomplete['accession'] == seq_record.id) & (Fivecomplete['gene'] == '5_8S_rRNA')]
            fivestart = int(startfive_df['seq_from'].iloc[0])
            fiveend = int(startfive_df['seq_to'].iloc[0])
            #print(seq_record.id, " FIVE ", fivestart, fiveend)
        if seq_record.id in LSU_RNA_df['accession'].tolist(): 
            startLSU_df = LSU_RNA_df[(LSU_RNA_df['accession'] == seq_record.id) & (LSU_RNA_df['gene'] == 'LSU_rRNA_eukarya')]
            LSUstart = int(startLSU_df['seq_from'].iloc[0])
            LSUend = int(startLSU_df['seq_to'].iloc[0])
            #print(seq_record.id, " LSU ", LSUstart, LSUend)
            
        if seq_record.id in SSU_RNA_df['accession'].tolist():
            if seq_record.id in Fivecomplete['accession'].tolist():
                if SSUend > fivestart:
                    print(seq_record.id, " ssu spans ", SSUstart, SSUend)
                    print(seq_record.id, " five spans ", fivestart, fiveend)
                    print(seq_record.id, "Misassembled sequence 1")
                    misassembled.append(s)
                    removeacc.append(seq_record.id)
        if seq_record.id in Fivecomplete['accession'].tolist():
            if seq_record.id in LSU_RNA_df['accession'].tolist():
                if fiveend > LSUstart:
                    print(seq_record.id, "Misassembled sequence 2")
                    print(seq_record.id, " lsu spans ", LSUstart, LSUend)
                    misassembled.append(s)
                    removeacc.append(seq_record.id)   
    SeqIO.write(misassembled, "misassembled_seqs", "fasta") 
    if seq_record.id not in removeacc:
        rna_found.append(s)
        SeqIO.write(rna_found, "rna_found.seqs", "fasta")   
        
#noncontig seq check
    if seq_record.id in SSU_RNA_df['accession'].tolist():
        if seq_record.id in LSU_RNA_df['accession'].tolist():
            if seq_record.id not in Fivecomplete['accession'].tolist():
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
        length = combined_df['Length'].iloc[0]
        
        """
        Generate feature table
        """   
        if seq_record.id in SSUcomplete['accession'].tolist(): 
            start_df = SSUcomplete[(SSUcomplete['accession'] == seq_record.id) & (SSUcomplete['gene'] == 'SSU_rRNA_eukarya')]
            start = start_df['seq_from'].iloc[0]
            stop = start_df['seq_to'].iloc[0]
            if (start == 1) & (stop == length):
                seq_record.description = "small subunit ribosomal RNA" 
            if seq_record.id in Five_RNA_df['accession'].tolist():
                seq_record.description = seq_record.description + " small subunit ribosomal RNA, internal transcribed spacer 1"
            elif seq_record.id not in Five_RNA_df['accession'].tolist():
                if stop != length:
                    seq_record.description = seq_record.description + " small subunit ribosomal RNA and internal transcribed spacer 1"
                    stop = ">" + str(length)       
        elif seq_record.id in SSUpartial['accession'].tolist():  
            #start_df = SSUpartial[(SSUpartial['accession'] == seq_record.id) & (SSUpartial['gene'] == 'SSU_rRNA_eukarya')]
            #start = start_df['seq_from'].iloc[0]
            #stop = start_df['seq_to'].iloc[0]
            #mdlfrom = SSUpartial['mdl_from'].iloc[0]
            #mdlto = SSUpartial['mdl_to'].iloc[0]
            start = SSUpartial['seq_from'].iloc[0]
            stop = SSUpartial['seq_to'].iloc[0]
                      
          
            if (start == 1) & (stop == length):
                seq_record.description = "small subunit ribosomal RNA" 
                start = "<" + str(start)
                stop = ">" + str(stop)  
                #add more scenarios for 5' and 3' partial
            if seq_record.id in Five_RNA_df['accession'].tolist():
                seq_record.description = seq_record.description + " small subunit ribosomal RNA, internal transcribed spacer 1"
                start = "<" + str(start)
            elif seq_record.id not in Five_RNA_df['accession'].tolist():
                start = "<" + str(start)
                if stop != length:
                    seq_record.description = seq_record.description + " small subunit ribosomal RNA and internal transcribed spacer 1"
                    stop = ">" + str(length)  

        if seq_record.id in Fivecomplete['accession'].tolist():  
            if seq_record.id not in SSU_RNA_df['accession'].tolist():
                start_df = Fivecomplete[(Fivecomplete['accession'] == seq_record.id) & (Fivecomplete['gene'] == '5_8S_rRNA')]
                start = start_df['seq_from'].iloc[0]
                if start > 1:
                    seq_record.description = seq_record.description + " internal transcribed spacer 1, 5.8S ribosomal RNA"
                    start = "<1"
                elif start == 1 & stop == length:
                    seq_record.description = "5.8S ribosomal RNA"
                else:
                    seq_record.description = seq_record.description + " 5.8S ribosomal RNA"

            if seq_record.id in SSU_RNA_df['accession'].tolist():
                seq_record.description = seq_record.description + ", 5.8S ribosomal RNA" 

            if seq_record.id not in LSU_RNA_df['accession'].tolist():
                stop_df = Fivecomplete[(Fivecomplete['accession'] == seq_record.id) & (Fivecomplete['gene'] == '5_8S_rRNA')]
                length_df = Fivecomplete[(Fivecomplete['accession'] == seq_record.id) & (Fivecomplete['gene'] == '5_8S_rRNA')]
                stop = stop_df['seq_to'].iloc[0]
                length = length_df['Length'].iloc[0]
                if stop < length:
                    seq_record.description = seq_record.description + ", and internal transcribed spacer 2"
                    stop = ">" + str(length)

        elif seq_record.id in Fivepartial['accession'].tolist():
            start_df = Fivepartial[(Fivepartial['accession'] == seq_record.id) & (Fivepartial['gene'] == '5_8S_rRNA')]
            start = start_df['seq_from'].iloc[0]
            #start = Fivepartial['seq_from'].iloc[0]
            stop = Fivepartial['Length'].iloc[0]
            mdlfrom = Fivepartial['mdl_from'].iloc[0]
            mdlto = Fivepartial['mdl_to'].iloc[0]
            if seq_record.id not in SSU_RNA_df['accession'].tolist():
                if start == 1:
                      seq_record.description = "5.8S ribosomal RNA"
                else:
                    seq_record.description = seq_record.description + " internal transcribed spacer 1 and 5.8S ribosomal RNA"            
            if seq_record.id in LSU_RNA_df['accession'].tolist():
                seq_record.description = seq_record.description + " 5.8S ribosomal RNA, internal transcribed spacer 2"
            elif seq_record.id not in LSU_RNA_df['accession'].tolist():
                seq_record.description = seq_record.description + " 5.8S ribosomal RNA and internal transcribed spacer 2"
            if mdlfrom != 1:
                start = "<" + str(start)
            if stop != mdlto:
                stop = ">" + str(stop)   
        elif seq_record.id in Five_RNA_df['accession'].tolist():
            print(seq_record.id, "Low scoring complete 5.8S gene found this is likely a bad sequence")
 
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
                if seq_record.id in Five_RNA_df['accession'].tolist():
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
        

stdout_fileno = sys.stdout   
#print(removeacc, "removed accessions")
file = open(outputfile, "w")
file.writelines(line)
file.close()
scriptTime = datetime.now()
print("script time is ", scriptTime)             
            