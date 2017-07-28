# C3S-LAA
Clustering of Circular Consensus Sequences (C3S) for Long Amplicon Analysis of PacBio Data


Overview of C3S-LAA 
================================================
To improve and extend the functionality of the long amplicon analysis (LAA) module from PacBio, we restructured the clustering approach for LAA based on the high quality circular consensus sequence data, grouping these reads based on the primer sequences used to amplify the DNA and the molecular barcodes used to track individuals. This directs error correction and consensus sequence analysis to be performed on sequences that should be the same, from a given sample, leading to improved accuracy of the outputted data. In addition, integration of Minimus (Sommer et al. BMC Bioinformatics, 2007 8:64) for automated assembly of any overlapping amplicon consensus sequences allows for efficient processing of tiled amplicon resequence data.

Usage
================================================

Running C3S-LAA involves five steps. These are outlined below.

###  _1) Generate Circular Consensus Sequence (CCS) Reads_

Run the reads of insert protocol in SMRT Portal to generate circular concensus sequence (CCS) reads. 

###  _2) Set Up The Parameters File, Primer Information File, and Barcode Information File_

The required dependencies and input files along with the output directory need to be specified. The parameters file is used by the cluster.py script (step 3) and includes information about the location of the Minimus assembler, the location of output from step 1 (fofn and ccs files), the location of a barcode (optional) and a primer pairs file (required), specific run parameters, and settings information for automated generation of a torque script.

Example of parameters.py file:

    #########################################
    ### User input parameters for C3S-LAA ###
    #########################################
    
    ### Required dependencies and input files
      # path for the AMOS package that contains minimus assembler
    amos_path = "/usr/local/amos/bin/"
      # path to the primer pair info file
    primer_info_file = "primer_pairs_info.txt"
       # path to the barcode info file
    barcode_info_file = "barcode_pairs_info.txt"
      # path to PacBio file of file names (fofn) directing to the raw reads
    fofn = "/path/m160901_060459_42157_c101086112550000001823264003091775_s1_p0.bas.h5"
      # path to ccs reads
    ccs = "/path/reads_of_insert.fasta"
      # directory where the consensus files will be saved
    consensus_output = "./output/"
    
    ### C3S-LAA parameters
      # number of bases corresponding to padding + barcode that need to be trimmed from the amplicon consensus
    trim_bp = 21
      # 1: yes; 0: no
    barcode_subset = 0
      # reads >= "min_read_length" will be searched for the presence of primer sequences
    min_read_length = 0
      # 1: filter; 0: no filter
    min_read_len_filter = 1
      # searches for the primer sequence within n bases from the read terminals
    primer_search_space = 100
      # Maximum barcode seq length
    max_barcode_length = 0
      # Maximum padding seq length
    max_padding_length = 5

    ### torque script settings
      # walltime for consensus calling
    walltime = 190
      # node no./name for consensus calling
    node = "1"
      # no. of processors for consensus calling
    processors = 12
      # consensus sequences generated from >= "no_reads_threshold" will be used for assembly
    no_reads_threshold = 100
    
Example of barcode_info_file (tab delimited text file): 
    
    f_barcode_name	f_barcode_sequence	r_barcode_name	r_barcode_sequence
    Sample1_BC.F1	CTATACATGACTCTGC	Sample1_BC.R1	GCAGAGTCATGTATAG
    Sample2_BC.F1	CTATACATGACTCTGC	Sample2_BC.R2	CATGTACTGATACACA
    Sample3_BC.F1	CTATACATGACTCTGC	Sample3_BC.R3	GAGAGACGATCACATA
    Sample4_BC.F1	CTATACATGACTCTGC	Sample4_BC.R4	CTGATATGTAGTCGTA

Example of primers_info_file (tab delimited text file): 

    f_primer_name           r_primer_name           f_primer_sequence           r_primer_sequence
    TA_1_25390617_27_F	TA_1_25395472_24_R	AAACATTGGTGTGGAAAGCAACTGAAG	AGGGTCACAGCACAGGACAGATTC
    TA_1_25391952_24_F	TA_1_25396540_27_R	AGGGACAACGTAGGGAGCCTTTGG	CGTCGACCACCGAATCAAGCAAGCATG
    TA_2_37562840_25_F	TA_2_37567441_24_R	GGGTGTTGTTCGGTCACCTCCTTTG	ATCCTTTGAGTGACTGAGGGTGTG
    TA_2_37564580_25_F	TA_2_37569533_24_R	TACGAGGTGTTTGGTTTGGTGAACG	CATGCATGCACACCTTCCAAGCTC
    TA_3_33503980_28_F	TA_3_33507309_27_R	TGTCCACCGACGACATCATTGGAGAGTG	TGCACAGTGCACATATGCTGCTTGGTG
    TA_3_33506037_24_F	TA_3_33509162_25_R	CAAACCGCAGAGGATAGAGATCGC	GGGTTCTCATCAACATTTGGACCTC
    TA_6_7045710_25_F	TA_6_7050495_28_R	TAGGGAGAGGTGGGAATATAATGGG	CCATCAAGTACAACAACGCATGATCATC
    TA_6_7047707_27_F	TA_6_7052049_28_R	CAGCATGCGTATAAAGAAGGCGAGCTC	CCCGATGTGCGACGCCGTAACAAATCTC

** NOTE: Primers must be named using the following convention: Name_Chromosome_StartPos_PrimerLength_Direction

###  _3) Run Command 1_

        python cluster.py
    
The CCS reads are used to cluster the data based on the presence of barcodes (optional) and both the forward and reverse primer sequences for each amplicon (the pipeline considers the sense and antisense primer sequences). This will generate a list of CCS identifiers belonging to each primer pair cluster. The list is used to link the corresponding raw reads, using the whitelist option in LAA, to carry out Quiver-based consensus calling using only the raw reads belonging to an amplicon-specific cluster. Executing cluster.py will compile the required information on barcode and amplicon-specific read clusters and output a submission script for running Quiver using the TORQUE resource management system (https://en.wikipedia.org/wiki/TORQUE). If you don't use TORQUE, you can execute cluster.py and edit the TORQUE realted elements of the submission script for your specific platform.

###  _4) Run Command 2_

        qsub consensus_calling.sh

The shell script (consensus_calling.sh) will be generated by running "cluster.py" (from step 3). For users running this shell script as a TORQUE submission, the standard qsub parameters are indicated at the beginning of this shell script (first 6 lines as given below), based on the user input parameters. This may be removed/modified.
        
        #!/bin/sh
        #PBS -N consensus_calling
        #PBS -r n
        #PBS -l walltime=190:00:00
        #PBS -l nodes=1:ppn=12
        #PBS -d /path/
        
        <remainder of script [not shown] will include information for running Quiver on barcode and amplicon-specific clusters>

###  _5) Run Command 3_

        python consensus_assembly.py

This script first generates a multi-fasta file (merged_reads.fasta) containing the error corrected consensus sequences from all the amplicon for a given sample. Using this multi-fasta file, the script carries out assembly to generate an output file "merged_reads_assembly.fasta" that contains the assembly results. The assembled contigs are named as >1, >2, >3 etc. within this file.


###  _When multiple lanes of sequence data are available_
A custom file of file names (fofn) as a .txt file, containing the absolute path and file names of the raw reads should be made available. This should be depicted in the "fofn_pacbio_raw_reads" option in the parameters.py file.

Here is an example of the contents of the fofn.txt file, where xx_s1_p0.x.bax.h5 files are from lane 1 and xx_s2_p0.x.bax.h5 files are from lane 2:

    /absolute_path/m170410_233007_42157_c101187522550000001823244205011702_s1_p0.1.bax.h5
    /absolute_path/m170410_233007_42157_c101187522550000001823244205011702_s1_p0.2.bax.h5
    /absolute_path/m170410_233007_42157_c101187522550000001823244205011702_s1_p0.3.bax.h5
    /absolute_path/m170410_233007_42157_c101187522550000001823244205011702_s2_p0.1.bax.h5
    /absolute_path/m170410_233007_42157_c101187522550000001823244205011702_s2_p0.2.bax.h5
    /absolute_path/m170410_233007_42157_c101187522550000001823244205011702_s2_p0.3.bax.h5


Dependencies for C3S-LAA
================================================
C3S-LAA can be executed via Python. This requires installation of the following components. Other versions of these components have not been tested.

* Linux/Unix

* <a href="https://github.com/PacificBiosciences/SMRT-Analysis">SMRT-Analysis v2.3.0</a>
* <a href="http://python.org/">Python 2.7</a>
* <a href="http://www.numpy.org/">NumPy 1.9.2</a>
* <a href="http://pandas.pydata.org/">pandas 0.18.1</a>
* <a href="http://biopython.org/wiki/Download">Biopython 1.69</a>


Licensing and Availability
================================================
C3S-LAA is released under an MIT open source license.
The source code is available on GitHub: https://github.com/drmaize/C3S-LAA

