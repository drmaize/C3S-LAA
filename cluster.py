### C3S-LAA
### Code fo clustering ccs read IDs based on the presence of barcode and primer pairs
### Version 0.7.0: 05/09/2017
### Author: Felix Francis (felixfrancier@gmail.com), Randall J Wisser (rjw@udel.edu)

############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()

############################################################
#### IMPORT FUNCTIONS
############################################################

from Bio import SeqIO
import pandas as pd
import os
from parameters import *

############################################################
#### FUNCTIONS
############################################################

### function to get reverse complement of an input sequence
def rev_complement(seq):
    seq = seq.upper()
    basecomplement = {'A':'T', 
                      'C':'G', 
                      'G':'C', 
                      'T':'A', 
                      '-':'-', 
                      'N':'N'}
    letters = list(seq)
    letters = [basecomplement[base] for base in letters]
    complement = (''.join(letters))
    return complement[::-1]
	
### function for clustering reads
def primer_based_cluster(primer_info_file, ccs, node, processors, fofn, primer_search_space, max_barcode_length, max_padding_length):	
	if barcode_subset == 0:
		os.makedirs(consensus_output)
		df = pd.read_csv(primer_info_file, sep='\t', skiprows=0, header=0)
		adj_primer_search_space = primer_search_space + max_barcode_length + max_padding_length              
		with open("consensus_calling.sh", "w") as shell_script_output: 
			cwd = os.getcwd()
			shell_script_output.write("#!/bin/sh" + "\n")
			shell_script_output.write("#PBS -N consensus_calling" +"\n")
			shell_script_output.write("#PBS -r n"+ "\n")
			shell_script_output.write("#PBS -l walltime="+ str(walltime) +":00:00" + "\n")
			shell_script_output.write("#PBS -l nodes=" + node + ":ppn=" + str(processors) + "\n")
			shell_script_output.write("#PBS -d " + cwd +  "\n"+  "\n")
			shell_script_output.write("export SMRT=/opt/smrtanalysis" +  "\n")
			shell_script_output.write("export INPUT=" + fofn +  "\n")
			shell_script_output.write("alias smrtwrap='$SMRT/smrtcmds/bin/smrtwrap'" + "\n" + "\n")
			output_dir = "./output/"
			for index, row in df.iterrows():
				f_primer_seq, r_primer_seq =  str(row['f_primer_sequence']), str(row['r_primer_sequence'])
				f_primer_name, r_primer_name =  str(row['f_primer_name']), str(row['r_primer_name'])
				chr_no, amplicon_start, amplicon_stop =  f_primer_name.split("_")[1], f_primer_name.split("_")[2], r_primer_name.split("_")[2]
				rcf_primer_seq, rcr_primer_seq = rev_complement(f_primer_seq), rev_complement(r_primer_seq)  
				amplicon_name = str(f_primer_name) + "_" + str(r_primer_name)
				reads_file_name = output_dir + amplicon_name + "_"  "reads.txt"
				with open(reads_file_name, "w") as output: 
					shell_script_output.write("### laa error corrected consensus calling "+ amplicon_name + "_ reads" +  "\n")
					shell_script_output.write("smrtwrap ConsensusTools.sh AmpliconAnalysis $INPUT --noPhasing --noClustering -n " + str(processors) + " --whiteList " + reads_file_name + " -o " + output_dir + "amplicon_" + str(chr_no) + "_" + str(amplicon_start) + "_" + str(amplicon_stop) + "/" + "\n" + "\n")
					fasta_sequences = SeqIO.parse(open(ccs),'fasta')
					for fasta in fasta_sequences:
							read, header = str(fasta.seq), fasta.id
							if min_read_len_filter == 1:
								if len(read) >= min_read_length:
									if primer_search_space > 20:
										if f_primer_seq in read[:adj_primer_search_space] and rcr_primer_seq in read[-adj_primer_search_space:] :         
											output.write(header[:-4] + "\n")
										elif r_primer_seq in read[:adj_primer_search_space] and rcf_primer_seq in read[-adj_primer_search_space:] :       
											output.write(header[:-4] + "\n")
									else:
										print "Primer search space is too small (should be > 20 bp)"
										break
							
							
	elif barcode_subset == 1:	 
		os.makedirs(consensus_output)
		df_barcodes = pd.read_csv(barcode_list, sep='\t', skiprows=0, header=0)
		df_barcodes['rcf_barcode_sequence'] = df_barcodes.apply(lambda row: rev_complement(row['f_barcode_sequence']), axis=1)
		df_barcodes['rcr_barcode_sequence'] = df_barcodes.apply(lambda row: rev_complement(row['r_barcode_sequence']), axis=1)
		for index, row_b in df_barcodes.iterrows():
			f_barcode_name, r_barcode_name = str(row_b['f_barcode_name']), str(row_b['r_barcode_name'])
			barcode_dir = "./output/" + str(f_barcode_name) + "_" + str(r_barcode_name)+ "/"
			os.makedirs(barcode_dir)
		df = pd.read_csv(primer_info_file, sep='\t', skiprows=0, header=0)
		adj_primer_search_space = primer_search_space + max_barcode_length + max_padding_length              
		with open("consensus_calling.sh", "w") as shell_script_output: 
			cwd = os.getcwd()
			shell_script_output.write("#!/bin/sh" + "\n")
			shell_script_output.write("#PBS -N consensus_calling"+ "\n")
			shell_script_output.write("#PBS -r n"+ "\n")
			shell_script_output.write("#PBS -l walltime="+ str(walltime) +":00:00" + "\n")
			shell_script_output.write("#PBS -l nodes=" + node + ":ppn=" + str(processors) + "\n")
			shell_script_output.write("#PBS -d " + cwd +  "\n"+  "\n")
			shell_script_output.write("export SMRT=/opt/smrtanalysis" +  "\n")
			shell_script_output.write("export INPUT=" + fofn +  "\n")
			shell_script_output.write("alias smrtwrap='$SMRT/smrtcmds/bin/smrtwrap'" + "\n" + "\n")
			for index_b, row_b in df_barcodes.iterrows():
				f_barcode_seq, r_barcode_seq = str(row_b['f_barcode_sequence']), str(row_b['r_barcode_sequence'])
				rcf_barcode_seq, rcr_barcode_seq = str(row_b['rcf_barcode_sequence']), str(row_b['rcr_barcode_sequence'])
				f_barcode_name, r_barcode_name = str(row_b['f_barcode_name']), str(row_b['r_barcode_name'])
				barcode_dir = "./output/" + str(f_barcode_name) + "_" + str(r_barcode_name)+ "/"
				for index, row in df.iterrows():
					f_primer_seq, r_primer_seq =  str(row['f_primer_sequence']), str(row['r_primer_sequence'])
					f_primer_name, r_primer_name =  str(row['f_primer_name']), str(row['r_primer_name'])
					chr_no, amplicon_start, amplicon_stop =  f_primer_name.split("_")[1], f_primer_name.split("_")[2], r_primer_name.split("_")[2]
					rcf_primer_seq, rcr_primer_seq = rev_complement(f_primer_seq), rev_complement(r_primer_seq)  
					barcode_amplicon_name = str(f_barcode_name) + "_" + str(r_barcode_name) + "_" + str(f_primer_name) + "_" + str(r_primer_name)
					reads_file_name = barcode_dir + barcode_amplicon_name + "_"  "reads.txt"
					with open(reads_file_name, "w") as output: 
						shell_script_output.write("### laa error corrected consensus calling "+ barcode_amplicon_name + "_ reads" +  "\n")
						shell_script_output.write("smrtwrap ConsensusTools.sh AmpliconAnalysis $INPUT --noPhasing --noClustering -n " + str(processors) + " --whiteList " + reads_file_name + " -o " + barcode_dir + "amplicon_" + str(chr_no) + "_" + str(amplicon_start) + "_" + str(amplicon_stop) + "/" + "\n" + "\n")
						fasta_sequences = SeqIO.parse(open(ccs),'fasta')
						for fasta in fasta_sequences:
							read, header = str(fasta.seq), fasta.id
							if min_read_len_filter == 1:
								if len(read) >= min_read_length:
									if primer_search_space > 20:
										if f_primer_seq in read[:adj_primer_search_space] and rcr_primer_seq in read[-adj_primer_search_space:] :         
											if f_barcode_seq in read[:adj_primer_search_space] and rcr_barcode_seq in read[-adj_primer_search_space:] :        
												output.write(header[:-4] + "\n")
										elif r_primer_seq in read[:adj_primer_search_space] and rcf_primer_seq in read[-adj_primer_search_space:] :      
											if r_barcode_seq in read[:adj_primer_search_space] and rcf_barcode_seq in read[-adj_primer_search_space:] :  
												output.write(header[:-4] + "\n")
									else:
										print "Primer, barcode search space is too small (should be > 20 bp)"
										break


############################################################
#### CODE
############################################################

if __name__ == '__main__':
	primer_based_cluster(primer_info_file, ccs, node, processors, fofn, primer_search_space, max_barcode_length, max_padding_length)

############################################################
#Time to run the code: end timer
############################################################
t1 = time.time()
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))
print "total time to run = ", total, " seconds"



### Version log (SemVer format)
### 0.1.0: read mapping based
### 0.2.0: based on presence of primer sequence at the end of reads     
### 0.3.0: based on presence of primer sequence anywhere in the reads (useful for barcoded and/or padded reads, where the primers are not at the ends) 
### 0.4.0: Automated writing of shell commands for error correction and consensus calling for reads from each amplicon
### 0.5.0: Search for primer sequence within n bases from the read terminal. n is a user defined parameter.
### 0.6.0: use external parameter file
### 0.7.0: with or without barcode option