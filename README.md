# CRISPR_Array_Reformat

**What it does?**

The Perl script CRISPRFileToGFF.pl converts CRISPR arrays predicted by CRT, PILER-CR, CRISPRFinder and CRISPRCasFinder to CRISPRDetect-GFF format. CRISPRDetect reports predicted CRISPR-arrays in GFF format, with lengths displayed in the score column (i.e. the CRISPRDetect-GFF format).

The GFF file can then be used for display, procesing or as input to other software, such as CRISPRStudio.

**Motivation**

The basic structure for a CRISPR array is a number of direct repeats (DR) separated by unique spacers, functional arrays will be associated with Cas genes on the same genome, often nearby.

Early software for CRISPR detection grew from programs that described repeats. More recently software has tried to capture more information e.g. related to CRISPR system or strand. This has been extended to considering evolutionary information, e.g. mutations in DR or sequences comparisons. Examples include CRISPRCasFinder, CRT, CRISPRFinder, and CRISRDetect.

There are a number of different systems available to predict CRISPR arrays and systems in genomes. Currently prediction programs have a range of outputs. A common standard would facilitate the downstream analysis of these arrays e.g. prediction of  targets, genome comparisons. Databases are available of CRISPRs but these also use different formats.

Here we define a draft standard that will facilitate interchange of data, and provide a conversion program. The extensible standard is the GFF format introduced by CRISPRDetect.

**Usage**

The script requires the Perl library "JSON", as CRISPRCasFinder reports predicted CRISPR arrays in JSON format. 

The script will detect file types automatically based on the contents of file headers, for example, PILER-CR files has the word "pilercr" in the header, whereas CRISPRCasFinder files always start with a "{" as part of the JSON syntax. Users do not need to specify the file type. 

At the moment, only CRT, pilercr, CRISPRCasFinder and CRISPRfinder formats are supported.

Example commands:

        perl CRISPRFileToGFF.pl -in sample_crt.txt -out crt_out.gff
        perl CRISPRFileToGFF.pl -in sample_ccf.json -out ccf_out.gff
        perl CRISPRFileToGFF.pl -in sample_cf.txt -out cf_out.gff
        perl CRISPRFileToGFF.pl -in sample_pilercr.txt -out pc_out.gff
        
**Notes on the different input formats**

For CRISPRCasFinder, all information are available, so displayed them as a GFF.

For CRT, PILER-CR and CRISPRFinder, strand and score are not available. So score = NA, and strand = "+" since data are extracted from the original sequence provided.

For PILER_CR, consensus repeat is given, but individual repeat sequences are not given. However, the coordinates of individual repeats are given. So the start-end coordinates were taken as given for repeats, but the sequences were taken as the consensus.

For CRISPRFinder, consensus repeat is given, but neither individual repeat sequences nor their start-end coordinates are given. So calculate the start-end coordinates of repeats using the start-end coordinates of the entire array and spacers, and the sequences were taken as the consensus.

For CRT, individual repeat sequences are given along with their start-end coordinates, but the consensus repeat is not given​. So the consensus as the most-common repeat sequence (by making a tally using a hash map data structure).


