from Bio import SeqIO
import functools
import pandas as pd
configfile: "config.yaml"
gb_input = config['gb_input']
final_tblout = config['final_tblout']
split_size = config['split_size']
feature_table = config['feature_table']
rejects = config['rejects']

rule all:
    input: 
        'stripped.fsa',
        'found.fsa',
        final_tblout,
        'final.ribomaker.fa',
        'blast_fasta.fsa',

rule parse_input_gb:
    input:
        gb_in = gb_input,
        rej = rejects,
    output:
        stripped_fsa = 'stripped.fsa',
        found_fsa = 'found.fsa',
        seqlen = 'my.seqlen',
    run:
        ## create a list of reject accs 
        reject_accs = set()
        with open(input.rej, 'rt') as f:
            for line in f:
                reject_accs.add(line.strip().strip('\n'))

        ## parse genbank flatfile
        sequences = [] 
        found = []
        sequencelength = []

        for seq_record in SeqIO.parse(input.gb_in, "genbank"):    
            str_id = seq_record.id
            if seq_record.name not in reject_accs:                       
                seq_record.description = seq_record.annotations["organism"]
                sequences.append(seq_record)
                seqlength = '%s %i\n' %  (seq_record.id, len(seq_record))
                sequencelength.append(seqlength)
            else:
                found.append(seq_record)
                print("I found this accession on the reject list and wrote the sequence to found.fsa: ", seq_record.id)  
        SeqIO.write(sequences, output.stripped_fsa, "fasta")  
        SeqIO.write(found, output.found_fsa, "fasta") 
        #delimiter = ','
        seqlen_str = functools.reduce(lambda a,b : a + b, sequencelength) 
        f = open(output.seqlen, 'w')
        f.write(seqlen_str)
        f.close()  


checkpoint split_fasta:
    input:
        stripped_fsa = rules.parse_input_gb.output.stripped_fsa,
    output: directory('split_fasta_files')
    shell:
        '''
        /usr/local/seqkit/0.11.0/bin/seqkit split2 {input.stripped_fsa} -s {split_size} -O {output} -f
        '''
## output files are 
## stripped.part_001.fsa  stripped.part_002.fsa

rule run_ribomaker:
    input: 'split_fasta_files/stripped.part_{fapart}.fsa',
    output: 'ribomaker_out/ribomaker_{fapart}/ribomaker_{fapart}.ribodbmaker.final.fa',
    log: 'log_files/ribomaker_{fapart}.log'
    params:
        outprefix = 'ribomaker_out/ribomaker_{fapart}'
    shell:
        '''
        export RIBODIR="/panfs/pan1/dnaorg/ssudetection/code/ribovore-install/ribovore"
        export RIBOINFERNALDIR="/usr/local/infernal/1.1.2/bin"
        export RIBOEASELDIR="/usr/local/infernal/1.1.2/bin"
        export RIBOTIMEDIR="/usr/bin"
        export SENSORDIR="/panfs/pan1/dnaorg/ssudetection/code/ribovore-install/rRNA_sensor"
        export EPNOPTDIR="/panfs/pan1/dnaorg/ssudetection/code/ribovore-install/epn-options"
        export EPNOFILEDIR="/panfs/pan1/dnaorg/ssudetection/code/ribovore-install/epn-ofile"
        export EPNTESTDIR="/panfs/pan1/dnaorg/ssudetection/code/ribovore-install/epn-test"
        export PERL5LIB="$RIBODIR":"$EPNOPTDIR":"$EPNOFILEDIR":"$EPNTESTDIR"
        export PATH="$RIBODIR":"$SENSORDIR":"$PATH"
        export BLASTDB="$SENSORDIR"
        export RIBOBLASTDIR="/panfs/pan1/dnaorg/ssudetection/code/ribovore-install/ncbi-blast-2.8.1+/bin"
        export VECPLUSDIR="/panfs/pan1/dnaorg/ssudetection/code/ribovore-install/vecscreen_plus_taxonomy"
        export PERL5LIB="$PERL5LIB":"/usr/local/perl/5.16.3/lib/perl5"
        export BLASTDB="/panfs/pan1/dnaorg/ssudetection/code/ribovore-install/vecscreen_plus_taxonomy/univec-files/":"$SENSORDIR"


        ## create a directory for ribomaker output
        mkdir -p 'ribomaker_out'
        mkdir -p {params.outprefix}

        ## run ribomaker.pl script
        /panfs/pan1.be-md.ncbi.nlm.nih.gov/dnaorg/ssudetection/code/ribotyper-v1/ribodbmaker.pl -f \
          --skipfribo1 --skipfribo2 --skipfmspan \
          --skipingrup --skipclustr --skiplistms \
          --skipmstbl --skipfblast --skipftaxid \
          {input} {params.outprefix} > {log} 2>&1

        ## touch output file (workaround because ribodbmaker output is a folder)
        touch {output}
        '''

## output file we are interested in:
## ribomaker_out/ribomaker_001/ribomaker_001.ribodbmaker.final.fa


rule run_cmscan_df:
    input: 
        rm_fa = rules.run_ribomaker.output[0],
    output: 'cm_out/cm_tblout_{fapart}.df.txt'
    params: 
        cmdb = 'rrna.cm'
    threads: 4
    log: 'log_files/cmscan_df_{fapart}.log'
    shell:
        '''
        mkdir -p 'cm_out'

        /usr/local/infernal/1.1.2/bin/cmscan --cpu {threads} \
          --mid -T 20 --verbose \
          --tblout {output} \
          {params.cmdb} \
          {input.rm_fa} > /dev/null 
        ''' 

rule run_cmscan_at:
    input: 
        rm_fa = rules.run_ribomaker.output[0],
    output: 'cm_out/cm_tblout_{fapart}.at.txt'
    params: 
        cmdb = 'rrna.cm'
    threads: 4
    log: 'log_files/cmscan_at_{fapart}.log'
    shell:
        '''
        mkdir -p 'cm_out'
        
        /usr/local/infernal/1.1.2/bin/cmscan --cpu {threads} \
          --mid -T 20 --verbose --anytrunc \
          --tblout {output} \
          {params.cmdb} \
          {input.rm_fa} > /dev/null 
        ''' 

rule cm_deoverlap:
    input: 
        [ 'cm_out/cm_tblout_{fapart}.df.txt', 
          'cm_out/cm_tblout_{fapart}.at.txt' ]
    output: 
        combined_tblout = 'cm_out/cm_tblout_{fapart}.both.txt',
        deoverlap_tbl = touch('cm_out/cm_tblout_{fapart}.both.txt.deoverlapped'),
    shell:
        '''
        cat {input} > {output.combined_tblout}
        perl cmsearch_tblout_deoverlap/cmsearch-deoverlap.pl \
          --maxkeep -s --cmscan {output.combined_tblout}
        '''

def create_file_list(wildcards):
    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    faparts = glob_wildcards('split_fasta_files/stripped.part_{fapart}.fsa').fapart
    file_list = expand('cm_out/cm_tblout_{fapart}.both.txt.deoverlapped', fapart=faparts)
    print('final cm_tblout files: ', file_list)
    return file_list

rule aggregate_tblout_files:
    input: create_file_list
    output: 'agg_final.tblout' 
    shell:
        '''
        cat {input} > {output} 
        '''

rule tblout_add:
    input: 
        seqlen = 'my.seqlen',  
        finaltbl = 'agg_final.tblout'
    output: final_tblout
    shell:
        '''
        perl tblout-add.pl -t {input.finaltbl} 18 {input.seqlen} 3 > {output}
        '''       

def create_ribodb_file_list(wildcards):
    checkpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
    faparts = glob_wildcards('split_fasta_files/stripped.part_{fapart}.fsa').fapart
    file_list = expand('ribomaker_out/ribomaker_{fapart}/ribomaker_{fapart}.ribodbmaker.final.fa', fapart=faparts)
    print('final ribomaker_tblout files: ', file_list)
    return file_list
    
rule aggregate_ribodbmaker_files:
    input: create_ribodb_file_list
    output: 'final.ribomaker.fa' 
    shell:
        '''
        cat {input} > {output}        
        '''   

rule runParser:
    input:
        final_tblout,
        'final.ribomaker.fa',
    output:
         config['blast_fsa']
    shell:
        '''
        /home/mcveigh/master/bin/python3 ParseCMscan1.67.py {output}
        '''