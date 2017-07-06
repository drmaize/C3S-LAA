# C3S-LAA
Clustering of Circular Consensus Sequences (C3S) for Long Amplicon Analysis of PacBio Data


##########################
### Overview of C3S-LAA ###
##########################
To improve and extend the functionality of the long amplicon analysis (LAA) module from PacBio, we restructured the clustering approach for LAA based on the high quality circular consensus sequence data, grouping these reads based on the primer sequences used to amplify the DNA and the molecular barcodes used to track individuals. This directs error correction and consensus sequence analysis to be performed on sequences that should be the same, from a given sample, leading to improved accuracy of the outputted data. In addition, integration of Minimus (Sommer et al. BMC Bioinformatics, 2007, 8:64) for automated assembly of any overlapping amplicons allows for efficient processing of tiled amplicon resequence data.

