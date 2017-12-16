### C3S-LAA
### Code for demultiplexing and assembly of error corrected PacBio read clusters. 
### Version 0.2.2: 05/09/2017
### Author: Felix Francis (felixfrancier@gmail.com), Randall J Wisser (rjw@udel.edu)

############################################################
#Time to run the code: start timer
############################################################
import time
t0 = time.time()

############################################################
#### IMPORT FUNCTIONS
############################################################

import multiprocessing
from Bio import SeqIO
import pandas as pd
import datetime
import re
import os
from pydna import Assembly, Dseqrecord
import subprocess as sp  
import os.path
import time
import shutil
from parameters import *

############################################################
#### FUNCTIONS
############################################################

### convert a given fasta file containing sequence reads to afg format and carry out minimus assembly.
def minimus_assembly(barcode_subset, amos_path, consensus_output):
	### convert fasta to afg format
	if barcode_subset == 1:
		reads = consensus_output + barcode_name + "/" + barcode_name + "merged_reads"
	else:
		reads = consensus_output + "merged_reads"
	p = sp.Popen(["%stoAmos" %amos_path,"-s","%s" %reads+".fasta", "-o", "%s" %reads+"_assembly_temp.afg"], stdout=sp.PIPE)
	out, err = p.communicate()
	### assemble afg reads
	p = sp.Popen(["%sminimus" %amos_path, "%s" %reads+"_assembly_temp.afg"], stdout=sp.PIPE)
	out, err = p.communicate()
	shutil.rmtree(reads+ "_assembly_temp.bnk")
	os.remove(reads+ "_assembly_temp.afg.runAmos.log")
	os.remove(reads+ "_assembly_temp.afg") 
	os.remove(reads+ "_assembly_temp.contig") 
	renamed_assembly = open(consensus_output +"c3slaa_assembly" +".fasta", "w")
	fasta_assembly = SeqIO.parse(open(consensus_output + "merged_reads_assembly_temp.fasta"), 'fasta') 
	for fasta in fasta_assembly:
		header, sequence = fasta.id.split('_'), str(fasta.seq)
		renamed_assembly.write(">C3SLAA_" +str(header[0]) + "\n")
		renamed_assembly.write(sequence + "\n")
	renamed_assembly.close()
	os.remove(consensus_output +"merged_reads_assembly_temp.fasta") 
	
### convert a given fasta file containing sequence reads to afg format and carry out minimus assembly (barcoded)
def minimus_assembly_barcoded(barcode_name, barcode_subset, amos_path, consensus_output):
	### convert fasta to afg format
	if barcode_subset == 1:
		reads = consensus_output + barcode_name + "/" + barcode_name + "merged_reads"
	else:
		reads = consensus_output + "merged_reads"
	p = sp.Popen(["%stoAmos" %amos_path,"-s","%s" %reads+".fasta", "-o", "%s" %reads+"_assembly_temp.afg"], stdout=sp.PIPE)
	out, err = p.communicate()
	### assemble afg reads
	p = sp.Popen(["%sminimus" %amos_path, "%s" %reads+"_assembly_temp.afg"], stdout=sp.PIPE)
	out, err = p.communicate()
	shutil.rmtree(reads+ "_assembly_temp.bnk")
	os.remove(reads+ "_assembly_temp.afg.runAmos.log")
	os.remove(reads+ "_assembly_temp.afg") 
	os.remove(reads+ "_assembly_temp.contig") 
	renamed_assembly = open(consensus_output + barcode_name+"c3slaa_assembly" +".fasta", "w")
	fasta_assembly = SeqIO.parse(open(consensus_output + barcode_name+ "merged_reads_assembly_temp.fasta"), 'fasta') 
	for fasta in fasta_assembly:
		header, sequence = fasta.id.split('_'), str(fasta.seq)
		renamed_assembly.write(">C3SLAA_" +str(header[0]) + "\n")
		renamed_assembly.write(sequence + "\n")
	renamed_assembly.close()
	os.remove(consensus_output + barcode_name +"merged_reads_assembly_temp.fasta") 
	
	
	
	
	
### pick reads from each barcoded sample and assemble (barcoded)
def pooled_locus_read_assembly_barcoded(primer_info_file, barcode_subset, amos_path,consensus_output):
	df_barcodes = pd.read_csv(barcode_list, sep='\t', skiprows=0, header=0)
	for index, row in df_barcodes.iterrows():
		barcode_name = str(row['f_barcode_name']) + "_" + str(row['r_barcode_name'])
		merged_fasta = open(str(consensus_output) + barcode_name + "/"+ barcode_name + "merged_reads" +".fasta", "w")
		df_primers = pd.read_csv(primer_info_file, sep='\t', skiprows=0, header=0)  ### read primer info
		for index, row in df_primers.iterrows():
			f_primer_name, r_primer_name =  str(row['f_primer_name']), str(row['r_primer_name'])
			primer_chr_no, amplicon_start, amplicon_stop = int(f_primer_name.split("_")[1]), int(f_primer_name.split("_")[2]), int(r_primer_name.split("_")[2])
			seq_path = str(consensus_output) + barcode_name + "/amplicon_"+ str(primer_chr_no) + "_" +  str(amplicon_start) + "_" + str(amplicon_stop) + "/"
			fasta_sequences = SeqIO.parse(open(seq_path + "amplicon_analysis.fasta"), 'fasta') 
			for fasta in fasta_sequences:
				header, sequence = fasta.id.split('_'), str(fasta.seq)
				no_reads = int(header[3][8:])
				if no_reads >= no_reads_threshold:
					merged_fasta.write(">" +str(barcode_name) + "_" +str(primer_chr_no) + "_" +str(amplicon_start)+ "_" +str(amplicon_stop) + "\n")
					merged_fasta.write(sequence[trim_bp:-trim_bp] + "\n")
		merged_fasta.close()
		minimus_assembly_barcoded(barcode_name, barcode_subset, amos_path, consensus_output) 

### pick reads from each sample and assemble (no barcodes)
def pooled_locus_read_assembly(primer_info_file, barcode_subset, amos_path,consensus_output):
	merged_fasta = open(str(consensus_output) + "merged_reads" +".fasta", "w")
	df_primers = pd.read_csv(primer_info_file, sep='\t', skiprows=0, header=0)  ### read primer info
	for index, row in df_primers.iterrows():
		f_primer_name, r_primer_name =  str(row['f_primer_name']), str(row['r_primer_name'])
		primer_chr_no, amplicon_start, amplicon_stop = int(f_primer_name.split("_")[1]), int(f_primer_name.split("_")[2]), int(r_primer_name.split("_")[2])
		seq_path = str(consensus_output) + "/amplicon_"+ str(primer_chr_no) + "_" +  str(amplicon_start) + "_" + str(amplicon_stop) + "/"
		fasta_sequences = SeqIO.parse(open(seq_path + "amplicon_analysis.fasta"), 'fasta') 
		for fasta in fasta_sequences:
			header, sequence = fasta.id.split('_'), str(fasta.seq)
			no_reads = int(header[3][8:])
			if no_reads >= no_reads_threshold:
				merged_fasta.write(">" +str(primer_chr_no) + "_" +str(amplicon_start)+ "_" +str(amplicon_stop) + "\n")
				merged_fasta.write(sequence[trim_bp:-trim_bp] + "\n")
	merged_fasta.close()
	minimus_assembly(barcode_subset, amos_path, consensus_output) 

############################################################
#### CODE
############################################################

if __name__ == '__main__':
	if barcode_subset == 1:
		pooled_locus_read_assembly_barcoded(primer_info_file, barcode_subset, amos_path,consensus_output)
	elif barcode_subset ==0:
		pooled_locus_read_assembly(primer_info_file, barcode_subset, amos_path,consensus_output)

#############################################
#Time to run the code: end timer
############################################################
t1 = time.time()
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))
print "total time to run = ", total, " seconds"

