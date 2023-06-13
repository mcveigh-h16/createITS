# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 08:39:05 2020

@author: mcveigh
"""

#
# processITS designed to validate ITS sequences. Sequences that pass all tests are saved to a outfile file in fasta format
# Selectively sort GenBank flatfiles removing selected seqeunces and saving them to a new fasta file
# User must specify the input filename and the outputfile name

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

#Read in the reject list and find the accession
reject_file_name = (r'ITS_reject_seqs3.txt') 
#rejectlist_df = pd.read_csv(reject_file_name, sep='\t', index_col=None, low_memory=False, header=None, names=["accession", "type", "reason"])
rejectlist_df = pd.read_csv(reject_file_name, sep='\t', index_col=None, low_memory=False, header=None, names=["accession"])
#rejectlist = rejectlist_df['accession']
reject_list = set(rejectlist_df['accession'].tolist())
   

#Parse the GenBank file and remove any sequences found on the reject list
from Bio import SeqIO
sequences = [] 
found = []
sequencelength = []
#missingRNA = []

rejectTime = datetime.now()
print("Reject time is ", rejectTime) 

for seq_record in SeqIO.parse(inputfile, "genbank"):    
    str_id = seq_record.id
    #print(seq_record.id)
    #if str_id.find('.') != -1:        
    #    str_id = str_id[:str_id.find('.')]        
    if seq_record.name not in reject_list:                       
        seq_record.description = seq_record.annotations["organism"]
        sequences.append(seq_record)
        seqlength = '%s %i\n' %  (seq_record.id, len(seq_record))
        sequencelength.append(seqlength)
    else:
        found.append(seq_record)
        print("I found this accession on the reject list and wrote the sequence to found.fsa: ", seq_record.id)  
SeqIO.write(sequences, "stripped.fsa", "fasta")  
SeqIO.write(found, "found.fsa", "fasta")  
#delimiter = ','
seqlen_str = functools.reduce(lambda a,b : a + b, sequencelength) 
f = open('my.seqlen', 'w')
f.write(seqlen_str)
f.close()
  
#Run short version of ribodbmaker.pl    
os.system("ribodbmaker -f --skipfribo1 --skipfribo2 --skipfmspan --skipingrup --skipclustr --skiplistms --skipmstbl --skipfblast --taxin /panfs/pan1/dnaorg/rrna/git-ncbi-rrna-project/taxonomy-files/ncbi-taxonomy-tree.ribodbmaker.txt stripped.fsa ribo-out")
#os.system("ribodbmaker -f --skipfribo1 --skipfribo2 --skipfmspan --skipingrup --skipclustr --skiplistms --skipmstbl --skipfblast --skipftaxid  stripped.fsa ribo-out")
#os.system("ribodbmaker -f --skipfribo1 --skipfribo2 --skipfmspan --skipingrup --skipclustr --skiplistms --skipmstbl --skipfblast stripped.fsa ribo-out")
riboTime = datetime.now()
print("ribodbmaker time is ", riboTime) 

#Run CMscan using the output of ribodbmaker.pl
print('cmscan1')
os.system("cmscan --cpu 16 --mid -T 20 --verbose --tblout tblout.df.txt rrna.cm ribo-out/ribo-out.ribodbmaker.full.fa > /dev/null")
print('cmscan2')
os.system("cmscan --cpu 16 --mid -T 20 --verbose --anytrunc --tblout tblout.at.txt rrna.cm ribo-out/ribo-out.ribodbmaker.full.fa > /dev/null")
cmscanTime = datetime.now()
print("CMscan time is ", cmscanTime) 
os.system("cat tblout.df.txt tblout.at.txt > tblout.both.txt")
os.system("perl cmsearch_tblout_deoverlap/cmsearch-deoverlap.pl --maxkeep -s --cmscan tblout.both.txt")
os.system("head -n2 tblout.both.txt > final.tblout")
os.system("cat tblout.both.txt.deoverlapped >> final.tblout")

#Add seq length to cmscan output
#os.system("esl-seqstat -a ribo-out/ribo-out.ribodbmaker.full.fa | grep ^\\= | awk '{ printf(\"%s %s\\n\", $2, $3); }' > my.seqlen")
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
#CMscan_df = CMscan_df[1:]
#print(CMscan_df.head(10))

#Remove 5S model rows
CMscan_df2 = CMscan_df[CMscan_df['gene'] != "5S_rRNA"]
#print(CMscan_df2.head(20))
#Find sequences on the minus strand

#Find sequences that have truncated models
#truncatedSSU = CMscan_df2[CMscan_df2['trunc'] == "5'&3'"]
truncated=CMscan_df2.loc[(CMscan_df2['gene'] != "5_8S_rRNA") & (CMscan_df2['trunc'] == "5'&3'")]
print("I found truncated models suggesting the presence of an intron\n ", truncated)
#truncatedLSU = CMscan_df2[CMscan_df2['trunc'] == "5'&3'"]
#print("I found truncated LSU models suggesting the presence of an intron\n ", truncatedLSU)

#Find sequences that do not pass CMscan tests
fail_test = CMscan_df2[CMscan_df2['Inc'] == "?"]
print("Sequences that have a ? are\n ", fail_test)

#show rows containing SSU or LSU
SSU_RNA_df = CMscan_df2[CMscan_df2['gene'] == "SSU_rRNA_eukarya"]
#print("these sequences have SSU\n ", SSU_RNA_df)

Five_RNA_df = CMscan_df2[CMscan_df2['gene'] == "5_8S_rRNA"]
FiveComplete = Five_RNA_df[Five_RNA_df['trunc'] == "no"]
FiveComplete['score'] = pd.to_numeric(FiveComplete['score'])
FiveCompleteStrongHit = FiveComplete[FiveComplete['score'] > 50]
#FiveFivePartial = Five_RNA_df[Five_RNA_df['trunc'] == "5'"]
#FiveThreePartial = Five_RNA_df[Five_RNA_df['trunc'] == "3'"]
#FiveBothPartial = Five_RNA_df[Five_RNA_df['trunc'] == "5'&3'"]
LSUpartial = CMscan_df2.loc[(CMscan_df2['gene'] == "LSU_rRNA_eukarya") & (CMscan_df2['trunc'] != "no")]
LSUcomplete = CMscan_df2.loc[(CMscan_df2['gene'] == "LSU_rRNA_eukarya") & (CMscan_df2['trunc'] == "no")]
SSUpartial = CMscan_df2.loc[(CMscan_df2['gene'] == "SSU_rRNA_eukarya") & (CMscan_df2['trunc'] != "no")]
SSUcomplete = CMscan_df2.loc[(CMscan_df2['gene'] == "SSU_rRNA_eukarya") & (CMscan_df2['trunc'] == "no")]
LSUendmissing = LSUpartial[(LSUpartial['strand'] == "+") & (LSUpartial['seq_to'] != LSUpartial['Length'])]
SSUendmissing = SSUpartial[(SSUpartial['strand'] == "+") & (SSUpartial['seq_from'] > 1)]
SSUendfound = SSUpartial[(SSUpartial['strand'] == "+") & (SSUpartial['mdl_from'] == 1)]
LSUendfound = LSUpartial[(LSUpartial['strand'] == "+") & (LSUpartial['mdl_to'] == 3401)]

five_minus_strand = FiveCompleteStrongHit[FiveCompleteStrongHit['strand'] != "+"]
five_minus_any = Five_RNA_df[Five_RNA_df['strand'] != "+"]
five_plus_any = Five_RNA_df[Five_RNA_df['strand'] == "+"]
LSUpartial_minus_strand = LSUpartial[LSUpartial['strand'] != "+"]
LSUcomplete_minus_strand = LSUcomplete[LSUcomplete['strand'] != "+"]
LSUpartial_plus_strand = LSUpartial[LSUpartial['strand'] == "+"]
LSUcomplete_plus_strand = LSUpartial[LSUpartial['strand'] == "+"]
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
print("I found 5.8S sequences on the minus strand\n ", five_minus_strand)
print("I found LSU sequences on the minus strand\n ", LSUminus)
print("I found SSU sequences on the minus strand\n ", SSUminus)

#print("LSU partial\n", LSUcomplete)
#print("SSU complete\n", SSUcomplete)
#print("LSU complete\n", LSUcomplete)
#print("These 5.8S rRNA are 5'&3'\n", FiveBothPartial)

#Sequences with extra sequence on the 5' the extends beyond position 1 of the SSU model
SSU_not5_end = SSU_RNA_df[(SSU_RNA_df['seq_from'] != 1) & (SSU_RNA_df['strand'] == "+")]
SSUextra = SSU_not5_end[SSU_not5_end['mdl_from'] == 1]
print("sequences with extra data on the 5' end \n", SSUextra)
#Sequnces with extra data on the 3' end beyond the end of LSU
LSU_RNA_df = CMscan_df2[(CMscan_df2['gene'] == "LSU_rRNA_eukarya") & (CMscan_df2['strand'] == "+")]
LSUextra=LSU_RNA_df.loc[(LSU_RNA_df['seq_to'] != LSU_RNA_df['Length']) & (LSU_RNA_df['mdl_to'] == 3401) & (LSU_RNA_df['mdl_from'] == 1)]
print("sequences with extra data on the 3' end \n", LSUextra)


#Trim the sequences with extra sequence flanking the SSU and LSU then rewrite the editted fasta to a new file
from Bio import SeqIO
sequences = []  
misassembled = []
removeacc = []
rna_not_found = []
rna_found = []
minus_strand_acc = []

for seq_record in SeqIO.parse("ribo-out/ribo-out.ribodbmaker.full.fa", "fasta"): 
    s = seq_record
    start = []
    to = []
    SSUend = []
    fivestart = []
    LSUstart = []
    fiveend = []
    if seq_record.id not in Five_RNA_df['accession'].tolist(): 
        #sequences.remove(seq_record)
        print("No 5.8S rRNA was found in ", seq_record.id)
        #missingRNA.append(seq_record)
    #elif seq_record.id in FiveBothPartial['accession'].tolist():
        #print("Truncated 5.8S rRNA was found in ", seq_record.id)
    #elif seq_record.id in FiveFivePartial['accession'].tolist():
        #print("5' partial 5.8S rRNA was found in ", seq_record.id)   
    #elif seq_record.id in FiveThreePartial['accession'].tolist():
        #print("3' partial 5.8S rRNA was found in ", seq_record.id)
    if seq_record.id in FiveCompleteStrongHit['accession'].tolist():
        if seq_record.id in LSUextra['accession'].tolist():           
            #start_df = CMscan_df2[(CMscan_df2['accession'] == seq_record.id) & (CMscan_df2['gene'] == 'SSU_rRNA_eukarya')]
            #start = start_df['seq_from'].iloc[0]
            to_df = CMscan_df2[(CMscan_df2['accession'] == seq_record.id) & (CMscan_df2['gene'] == 'LSU_rRNA_eukarya')]
            to = to_df['seq_to'].iloc[0] 
            s.seq = s.seq[0:int(to)]
            #seq_record.description = seq_record.description + " small subunit ribosomal RNA, internal transcribed spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2 and large subunit ribosomal RNA, complete sequence"
            #seq_record.description = seq_record.description + " COMPLETE LSU"
            #print(seq_record.id,seq_record.description)
        if seq_record.id in SSUextra['accession'].tolist():
            start_df = SSUextra[(SSUextra['accession'] == seq_record.id) & (SSUextra['gene'] == 'SSU_rRNA_eukarya')]
            start = start_df['seq_from'].iloc[0]
            to = start_df['Length'].iloc[0]
            s.seq = s.seq[start-1:to]
#trim overhanging unconfirmed sequences from partial predictions    
        if seq_record.id in LSUendmissing['accession'].tolist(): 
            if seq_record.id not in LSUendfound['accession'].tolist(): 
                to_df = CMscan_df2[(CMscan_df2['accession'] == seq_record.id) & (CMscan_df2['gene'] == 'LSU_rRNA_eukarya')]
                to = to_df['seq_to'].iloc[0]
                s.seq = s.seq[0:int(to)]
            
        if seq_record.id in SSUendmissing['accession'].tolist():  
            if seq_record.id not in SSUendfound['accession'].tolist(): 
                if seq_record.id not in SSUextra['accession'].tolist():
                    start_df = CMscan_df2[(CMscan_df2['accession'] == seq_record.id) & (CMscan_df2['gene'] == 'SSU_rRNA_eukarya')]
                    start = start_df['seq_from'].iloc[0]
                    to = start_df['Length'].iloc[0]
                    s.seq = s.seq[start-1:to]


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
                    print(seq_record.id, " lsu spans ", LSUstart, LSUend)
                    misassembled.append(s)
                    removeacc.append(seq_record.id)
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
                if SSUstart < fiveend:
                    print(seq_record.id, " ssu spans ", SSUstart, SSUend)
                    print(seq_record.id, " five spans ", fivestart, fiveend)
                    print(seq_record.id, "Misassembled sequence 1")
                    misassembled.append(s)
                    removeacc.append(seq_record.id)
        if seq_record.id in Five_RNA_df['accession'].tolist():
            if seq_record.id in LSU_RNA_df['accession'].tolist():
                if fiveend < LSUend:
                    print(seq_record.id, "Misassembled sequence 2")
                    print(seq_record.id, " lsu spans ", LSUstart, LSUend)
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

    
#Definition Line Generator               
        if seq_record.id not in removeacc:
            if seq_record.id in FiveCompleteStrongHit['accession'].tolist():
                if seq_record.id in five_minus_strand['accession'].tolist():
                    print("I reverse complemented ", seq_record.id)
                    s.seq = s.seq.reverse_complement()
                if seq_record.id in SSUcomplete['accession'].tolist():  
                    seq_record.description = seq_record.description + " SSU"
                elif seq_record.id in SSUendfound['accession'].tolist():   
                    seq_record.description = seq_record.description + " <SSU"
                elif seq_record.id in SSUpartial['accession'].tolist(): 
                    if seq_record.id not in SSUendfound['accession'].tolist():
                        seq_record.description = seq_record.description + " <SSU"
            #if seq_record.id in FiveComplete['accession'].tolist(): 
                if seq_record.id in SSU_RNA_df['accession'].tolist(): 
                    seq_record.description = seq_record.description + " ITS1 5.8S ITS2"
                else: 
                    seq_record.description = seq_record.description + " <ITS1 5.8S ITS2"       
                if seq_record.id in LSUcomplete['accession'].tolist():  
                    seq_record.description = seq_record.description + " LSU"
                elif seq_record.id in LSUpartial['accession'].tolist():  
                    seq_record.description = seq_record.description + " LSU>"
                else:
                    seq_record.description = seq_record.description + ">"
                sequences.append(s)
SeqIO.write(sequences, outputfile, "fasta")  
#SeqIO.write(missingRNA, "missing5_8.fsa", "fasta")

#Print seq_record.id and sequence length from output file
from Bio import SeqIO
for record in SeqIO.parse(outputfile, "fasta"):  
    tmp = LSUextra[LSUextra['accession'] == record.id]
    tmp2 = SSUextra[SSUextra['accession'] == record.id]
    tmp3 = LSUendmissing[LSUendmissing['accession'] == record.id]
    tmp4 = SSUendmissing[SSUendmissing['accession'] == record.id]
    tmpframe = [tmp, tmp2, tmp3, tmp4]
    extra = pd.concat(tmpframe)
    if not extra.empty:
        length = len(record)
        print("Sequences that were trimmed are: ", record.id, length)
        
#Count the number of sequences in the final output file
fh = open(outputfile)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
fh.close()
print("The number of sequences in the final output file is: ", n)


print ("processITS script is done and output is saved in:", outputfile)