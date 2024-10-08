This folder contains materials and code pertaining to the Rloop/HUV integration sites paper.

Most of the code is built upon previous work by Maja Kuzman and Dunja Glavaš, and can be found here;
https://github.com/MaKuzman/IS_mapping
https://zenodo.org/records/7529755
Simulating random genomic fragments was done as described by Greg Bedwell; 
https://gbedwell.github.io/random-genome-fragments/

Folders raw_integrations_1 and 2 contain the raw sequencing data for the HIV integration sites.
The respective metadata sheets for the samples are in NGS_Metadata_Sheet_1 and 2.
The folder "integration_analysis_scripts" contains the trimming/fltering/mapping bash scripts.

The "code" folder contains;

Random_insertions.R - code for making the random integration sites to compare with actual integration site data along with the code for making the comparison graphs.
MapInts.R - code for mapping the integration site data we get from the raw sequencing files onto the genome and comparing that data to previous data for rloop positions
