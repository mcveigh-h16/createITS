from Bio import SeqIO
import functools
import pandas as pd
configfile: "config.yaml"
gb_input = config['gb_input']
final_tblout = config['final_tblout']
split_size = config['split_size']
feature_table = config['feature_table']
#outdir = config['outdir']

rule all:
    input: 
        'input.fsa',
        final_tblout,
        #split_size,
        feature_table,
    #output:
        #outdir,

rule parse_input_gb:
    input:
        gb_in = gb_input,
    output:
        input_fsa = 'input.fsa',
        seqlen = 'my.seqlen',
    run:
        ## parse genbank flatfile
        sequences = [] 
        orgname = []
        sequencelength = []
        for seq_record in SeqIO.parse(input.gb_in, "fasta"):    
            str_id = seq_record.id
            sequences.append(seq_record)
            seqlen = '%s %i\n' %  (seq_record.id, len(seq_record))
            sequencelength.append(seqlen)
            if seq_record.seq.count("NNNNN"):
                print(seq_record.id, "contains internal Ns, this may result in incorrect predictions")
                with open('log_files/internalNs.log', 'wt') as f:
                    print(seq_record.id, 'contains internal Ns', file=f)
            #orgname = seq_record.annotations["organism"]
            #if "Giardia" in str(orgname):
                #print(seq_record.id, "From Giarda sp. and will require special handling")  
                #with open('log_files/Giardia.log', 'wt') as f:
                    #print(seq_record.id, ' Giardia sp.', file=f)
        SeqIO.write(sequences, output.input_fsa, "fasta")  
        seqlen_str = functools.reduce(lambda a,b : a + b, sequencelength) 
        f = open(output.seqlen, 'w')
        f.write(seqlen_str)
        f.close()  


checkpoint split_fasta:
    input:
        input_fsa = rules.parse_input_gb.output.input_fsa,
    output: directory('split_fasta_files')
    shell:
        '''
        /usr/local/seqkit/0.11.0/bin/seqkit split2 {input.input_fsa} -s {split_size} -O {output} -f
        '''
## output files are 
## input.part_001.fsa  input.part_002.fsa

rule run_cmscan_df:
    input: 
        rm_fa = rules.parse_input_gb.output[0],
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
        rm_fa = rules.parse_input_gb.output[0],
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
    faparts = glob_wildcards('split_fasta_files/input.part_{fapart}.fsa').fapart
    file_list = expand('cm_out/cm_tblout_{fapart}.both.txt.deoverlapped', fapart=faparts)
    print('final cm_tblout files: ', file_list)
    return file_list

rule aggregate_tblout_files:
    input: create_file_list
    output: 'agg_final.tblout' 
    log: 'log_files/aggregate.log'
    shell:
        '''
        cat {input} > {output} 2> {log}
        '''

rule tblout_add:
    input: 
        seqlen = 'my.seqlen',  
        finaltbl = 'agg_final.tblout'
    output: final_tblout
    log: 'log_files/tblout.log'
    shell:
        '''
        perl tblout-add.pl -t {input.finaltbl} 18 {input.seqlen} 3 > {output} 2> {log}
        '''       

rule runParser:
    input:
        final_tblout,
        'input.fsa',
    output:
         config['feature_table']
    log: 'log_files/parser.log'
    shell:
        '''
        /home/mcveigh/master/bin/python3 CMscantoITSfeatures.py {output} 2> {log}
        '''