# createITS
CreateITSfeatures1.0.py
Python application to analyze a gbk or fasta file containing eukaryotic nuclear ribosomal RNA and/or internal transcribed spacers (ITS1/ITS2) and generate a five column feature table that can be imported into gbench. 
Is it designed to work with any of these scenarios:
SSU
SSU-ITS1
SSU-ITS1-5.8
SSU-ITS1-5.8-ITS2
SSU-ITS1-5.8-ITS2-LSU
ITS1-5.8
ITS1-5.8-ITS2
ITS1-5.8-ITS2-LSU
5.8
5.8-ITS2
5.8-ITS2-LSU
ITS2-LSU
LSU
Plus or minus strand and partial or complete for every possible combination. 
Usage is
python createITSfeatures1.65.py inputfile.fasta outfile.tbl
there are additional files created by the script and saved in your directory. Some are intermediate files the script needs others are there for your reference if you need to investigate anything further. These include
cmscan_final.tblout	cmscan results file
rna_found.seqs		output fasta for any sequences where RNA was found
rna_not_found_seqs	output fasta for any sequences where an rRNA was not found
misassembled_seqs	output fasta for any sequences found to be misassembled
final.tblout		intermediate file
input.fsa  		intermediate file
my.seqlen		intermediate file
tblout.at.txt		intermediate file
tblout.both.txt		intermediate file
tblout.both.txt.deoverlapped		intermediate file
tblout.both.txt.sort	intermediate file
tblout.df.txt		intermediate file


SET UP
Configure your linux environment for python
This is a linux command line application and it requires some set up to install. 
Running Python at NCBI. 
https://confluence.ncbi.nlm.nih.gov/display/PY/Python+at+NCBI
Follow the instructions for checking your .ncbi_hints file and Creating a virtualenv
You must create a virtual environment and activate it daily when you want to use python.
Activation is something like this   > source myenv/bin/activate.csh
Your linux cursor will now change to look like this
((myenv) iebdev11[998]
Create the virtual environment in your home directory and it will be available in all your directories. If you create it in a sub-directory then it will only work there. 
Python Dependencies that must be installed before initial use. 
Pip install the following modules
Pandas
biopython
Datetime
For example >pip install pandas

Configure your linux environment to run CMscan (skip this if was previously done)
It is easiest if you run CMscan and the script in the same directory all the time. Recommend you set up a new directory for this. 
CMscan instructions and dependencies
http://eddylab.org/infernal/Userguide.pdf
Set up your linux environment to use CMscan and test that it works
Add infernal to the facilities line of your .ncbi_hints file
Download the Rfam models
If so, download and gunzip Rfam.cm.gz from here: ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz

Then create a text file rrna.list.txt with 4 lines:
RF00001
RF00002
RF01960
RF02543

And do these commands:
> cmfetch --index Rfam.cm
> cmfetch -f Rfam.cm rrna.list.txt > rrna.cm
> cmpress rrna.cm

Execute a cmscan test job in your working directory. Do both these commands to be sure they both work. 
> cmscan rrna.cm <fasta file>
>cmscan --anytrunc --tblout myseqs.tblout rrna.cm <fasta file>

copy these files and install them in the directory where you are going to do run the script. Be sure they are executable.
~mcveigh/createITS
cmsearch-deoverlap.pl
cmsearch_tblout_deoverlap
cmsearch_tblout_deoverlap-master
createITSfeatures.1.0.py

The script will write a number of intermediate files and some log files to your directory. The file names are hard coded and any existing files with the same names (i.e. previous run) will be overwritten. You need to be sure you do not run out of space in your linux directory or something rude will happen. 
Files needed or written
Input file	starting input
Output file	output feature table
cmscan_final.tblout	complete cmscan table
misassembled_seqs	report of any misassembled sequences
rna_found.seqs		report of all sequences where an RNA gene was found
rna_not_found_seqs	report of all sequences where no RNA gene was found (i.e. ITS1 or ITS2 only)

intermediate files the script reads and writes but the user doesn?t need to look at.
my.seqlen
final.tblout
input.fsa
tblout.at.txt
tblout.both.txt
tblout.both.txt.deoverlapped
tblout.both.txt.sort
tblout.df.txt

NOTES:
All Seqids must be unique! Any seqid is fine just be sure it is unique.
The feature table created as a misc_RNA with a contains note or a rRNA for single rRNA genes. 
CMscan is not a very fast tool. Processing time is dependent on the size of the input file. If you have a larger file but only need to evaluate a few sequences I suggest you just pull out those few. The script will run a file of any size but you may get a nastygram from systems if you try to run a large file on a shared machine such as iebdev11. Processing time on a shared machine will also vary with the load on the machine. Large files should be run on the farm and I can help you set that up if you need to. A multithreaded version of the script is planned and you will want to use this for really giant files. 
CMscan needs approximately 20 bp of the RNA gene to detect it. Similar to ITSx.
The script can only detect eukaryotic nuclear ribosomal RNAs (SSU, 5.8S and LSU). ITS only sequences and sequences without a RNA gene are reported and saved in rna_not_found_seqs with the original definition line. 
Known issues. If the sequence contains internal Ns, the script will detect this and report a warning message. It will still run but the results for any sequence with internal Ns maybe incorrect depending on the location and size of the internal Ns. 
If the RNA gene has introns it will work and create the appropriate misc_RNA. 
Giardia sp. if found are reported. Giardia have unusual rRNA gene structure and require special handling. The script will sometimes make a mistake in reporting whether the 5.8S rRNA is complete as this gene is smaller than normal in Giardia. Future version may include handling for Giardia but they are very rare so we added a warning message for now. 


