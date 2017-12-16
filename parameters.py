################################################################
### User input parameters for C3S-LAA
################################################################

### Required dependencies and input files

# Path for the AMOS package that contains minimus assembler
amos_path					= "/usr/local/amos/bin/"	

# C3S-LAA input files	
primer_info_file				= "primer_pairs_info.txt"
barcode_info_file			= "barcode_pairs_info.txt"
fofn								= "/mnt/data27/ffrancis/PacBio_sequence_files/EqPCR_raw/F03_1/Analysis_Results/m160901_060459_42157_c101086112550000001823264003091775_s1_p0.bas.h5"
ccs								= "/mnt/data27/ffrancis/PacBio_sequence_files/old/primer_pair_based_grouping/Eq_wisser_PCR-ccs-opt-smrtanalysis-userdata-jobs-020-020256-data-reads_of_insert.fasta"

# Output path
consensus_output			= "./output/"


### C3S-LAA parameters
trim_bp							= 21				# number of bases corresponding to padding + barcode that need to be trimmed from the amplicon consensus
barcode_subset				= 0				# (1: yes; 0: no) 
min_read_length			= 0				# reads >= "min_read_length" will be searched for the presence of primer sequences
min_read_len_filter		= 1             	# (1: filter; 0: no filter)
primer_search_space		= 100           # searches for the primer sequence within n bases from the read terminals
max_barcode_length		= 0     			# barcode seq length
max_padding_length		= 5     			# padding seq length


### torque script settings
walltime						= 190			# walltime for consensus calling
node								= "1"			# node no. for consensus calling
processors					= 12				# no. processors for consensus calling
no_reads_threshold		= 100			# consens sequences generated from >= "no_reads_threshold" will be used for assembly
