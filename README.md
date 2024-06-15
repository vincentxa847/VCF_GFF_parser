# Parsing variant data and generating summary information

### This program demonstrates fundamental bioinformatics skills by manipulating variant and genomic data from VCF and GFF files in Python and generating plots using matplotlib

## Running
The script accepts VCF, GFF, and FASTA files from the command line, with the option to specify the output file name.
```
python Manipulating_VCF_and_GFF files.py <filename>.vcf.gz <filename>.gff <filename>.fasta --output <outputFileName>
```

## Result
The script will iterate through variants stored in VCF file, if the SNP Quality of variant is greater than 20, then the coordinate of the variant will be extracted to index the CDS that the variant located in. 

This script handles VCF file without binary tabix-indexed, therefore VCF file with extension vcf.gz should be input not vcf.gz.tbi. 

CDS that the variant located is indexed using GFF file, the parent (mRNA) of the CDS will also be indexed. 
The DNA sequence of the CDSs under same parent will be extracted from the FASTA file, which will be utilized to get the DNA codons that the variant located in. 
In the DNA codon, the REF of the variant will be replaced by ALT of the variant. 
DNA codon before and after replaced will be extracted and translated to amino acid, which will be stored as Ref AA and Alt AA. 
If the replacement results in change in amino acid (Ref AA is not equal to Alt AA), then the variant will be marked as *non-synonymous*. 
Otherwise, the variant will be marked as *synonymous*. For the variant that SNP Quality is greater than 20 but the CDS of variant cannot be indexed (variant not in CDS), the variant will be marked as *non-coding*.

The output of the script is composed of three files, which are a log file (txt format), a bar plot (png format) and a tab-separated table (tsv format). 

**Log file** will record the filenames given at the command line, the number of the variants where SNP quality are smaller than 20 and the location of the output files (current directory). 
If the script run incompletely, then the filenames of the input files and the error will be recorded in log file. 
**Bar plot** will plot the proportion of variants with non-synonymous, synonymous and non-coding for variants where SNP quality are greater than 20. 
**Tab-separated table** will store the information for variants with SNP quality > 20, including *Chrom, Pos, Ref and Alt* from VCF file. 
*Type, Transcript, Protein Location, Ref AA and Alt AA* that calculated by script will also be stored in tab-separated table, but the Transcript, Protein Location, Ref AA and Alt AA will not be calculated for variants with type non-coding.

| CHROM | POS | REF | ALT | TYPE | TRANSCRIPT | Protein Location | Ref AA | Alt AA |
|----------|----------|----------|----------|----------|----------|----------|----------|----------|
| Pf3D7_10_v3	1468717 | 1468717  | C | A | Non-synonymous   | PF3D7_1037000.1   | 1948   | T   | N   |
| Pf3D7_10_v3	1509783 | 1509783  | G | T | Non-synonymous   | PF3D7_1038100.1   | 248    | P   | T   |

*Example of the output table*

![image](https://github.com/vincentxa847/VCF_GFF_parser/assets/118545004/d34f9e35-f1e0-4321-ad6d-07e00527a6e3)\
*Example of the output plot*
