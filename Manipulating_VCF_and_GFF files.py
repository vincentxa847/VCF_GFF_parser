import vcf
import os
import gffutils
import math
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('data' , help="Take in VCF,GFF and fasta file", nargs='*') # input more than one files
args = parser.parse_args()
VCF = ""
GFF = ""
FASTA = ""
for InputData in args.data:
    if str(InputData).find(".vcf") != -1 :
        VCF = str(InputData)
    elif str(InputData).find(".gff") != -1 :
        GFF = str(InputData)
    elif str(InputData).find(".fasta") != -1 :
        FASTA = str(InputData)
    else:
        error = "invalid data was input"
        with open("2802815logfile.txt", "w") as LogFile:
            LogFile.write("file given at the command line {} \nerror happened when running script, {}".format(str(args.data).replace("[","").replace("]",""),error))
        raise SystemExit(1)

if VCF == "" or GFF == "" or FASTA == "" :
    error = "incomplete file was input"
    with open("2802815logfile.txt", "w") as LogFile:
        LogFile.write("file given at the command line {} \nerror happened when running script, {}".format(str(args.data).replace("[", "").replace("]", ""), error))
    raise SystemExit(1)

# Handle gff file
if not os.path.isfile("2802815.db"):
    db = gffutils.create_db(GFF, dbfn="2802815.db")
else:
    db = gffutils.FeatureDB("2802815.db", keep_order=True)
# Read in vcf file
try:
    vcf_Reader1 = vcf.Reader(filename=VCF)
except UnicodeDecodeError:
    print("This script is not suit for tabix-indexed vcf file, please use vcf file without binary tabix-indexed")
    raise SystemExit(1)

VariantQualityGreaterThan20ID = []
# not all variant has CDS, targetVariant store the data of variant that has cds and mRNA
targetVariant = [] # DATA FROM VCF AND GFF
VariantsQualitySmallerThan20 = 0
for record in vcf_Reader1:
    if record.QUAL > 20:
        VariantQualityGreaterThan20ID.append(str(record.CHROM) + str(record.POS)) # using CHROM and POS as variant id
        for target in db.region(seqid=record.CHROM, start=record.POS, end=record.POS,featuretype="CDS"):
            for parent in db.parents(target, featuretype="mRNA"): # take the coordinate of mRNA to retrieve the DNA sequence
                transcript_id = parent.id
                Location = 0
                seq = ""
                if parent.strand == "+":
                    for cds in db.children(parent, featuretype="CDS", order_by="start"): # store the sequence for CDS under parent
                        seq = seq + cds.sequence(FASTA, use_strand="True")

                    for cds in db.children(parent, featuretype="CDS", order_by="start"): # store the coordinate of variant
                        if cds == target:
                            Location = Location + (int(record.POS) - cds.start + 1)
                            break
                        else:
                            cdsLength = cds.end - cds.start + 1
                            Location = Location + cdsLength
                elif parent.strand == "-":
                    for cds in db.children(parent, featuretype="CDS", order_by="start",reverse= True): # store the sequence for CDS under parent
                        seq = seq + cds.sequence(FASTA, use_strand="True")

                    for cds in db.children(parent, featuretype="CDS", order_by="start",reverse= True): # store the coordinate of variant
                        if cds == target:
                            Location = Location + (cds.end -int(record.POS) + 1)
                            break
                        else:
                            cdsLength = cds.end - cds.start + 1
                            Location = Location + cdsLength
                else:
                    error = "GFF file in bad format"
                    with open("2802815logfile.txt", "w") as LogFile:
                        LogFile.write("file given at the command line {} \nerror happened when running script, {}".format(str(args.data).replace("[", "").replace("]", ""), error))
                        raise SystemExit(1)
                ProteinLocation = math.ceil(Location / 3)
                if Location % 3  == 0: # variant in the third codon
                    REFDNACodon = seq[(Location-3):Location]
                    if parent.strand == "+":
                        ALTDNACodon = REFDNACodon[0] + REFDNACodon[1] + str(record.ALT)[1] # string plus string
                    elif parent.strand == "-":
                        ALTDNACodon = REFDNACodon[0] + REFDNACodon[1] + str(Seq(str(record.ALT)[1]).reverse_complement())
                elif Location % 3  == 2: # variant in the second codon
                    REFDNACodon = seq[(Location-2):(Location+1)]
                    if parent.strand == "+":
                        ALTDNACodon = REFDNACodon[0] + str(record.ALT)[1] + REFDNACodon[2]
                    elif parent.strand == "-":
                        ALTDNACodon = REFDNACodon[0] + str(Seq(str(record.ALT)[1]).reverse_complement()) + REFDNACodon[2]
                elif Location % 3 == 1: # variant in the first codon
                    REFDNACodon = seq[(Location-1):(Location+2)]
                    if parent.strand == "+":
                        ALTDNACodon = str(record.ALT)[1] + REFDNACodon[1] + REFDNACodon[2]
                    elif parent.strand == "-":
                        ALTDNACodon = str(Seq(str(record.ALT)[1]).reverse_complement()) + REFDNACodon[1] + REFDNACodon[2]
                RefAA = str(Seq(REFDNACodon).translate())
                AltAA = str(Seq(ALTDNACodon).translate())
                if RefAA == AltAA:
                    targetVariant.append([record.CHROM,str(record.POS),record.REF,str(record.ALT)[1],"Synonymous",transcript_id,str(ProteinLocation),RefAA,"NA"])
                elif RefAA != AltAA:
                    targetVariant.append([record.CHROM,str(record.POS),record.REF,str(record.ALT)[1],"Non-synonymous",transcript_id,str(ProteinLocation),RefAA,AltAA])
    else: VariantsQualitySmallerThan20 += 1

# find out the variant that not has cds
try:
    vcf_Reader2 = vcf.Reader(filename=VCF) # first vcf is exhausted, so need to open it again
except UnicodeDecodeError:
    print("This script is not suit for tabix-indexed vcf file,please use vcf file without binary tabix-indexed")
    raise SystemExit(1)


NonCodingVariant = [] # 2559
for k in targetVariant:  # remove variant that has cds
    for i in VariantQualityGreaterThan20ID:
        if i == str(k[0]) + str(k[1]):
            VariantQualityGreaterThan20ID.remove(i)
for record2 in vcf_Reader2:
    for NonCoding in VariantQualityGreaterThan20ID:
        if NonCoding == str(record2.CHROM) + str(record2.POS):
            NonCodingVariant.append([record2.CHROM, str(record2.POS), record2.REF, str(record2.ALT)[1], "Non-coding", "NA", "NA", "NA", "NA"])

# Output data
WholeData = targetVariant + NonCodingVariant
# For bar plot
SynonymousNumber = 0
NonSynonymousNumber = 0
NonCodingNumber = 0
for data in WholeData:
    if data[4] == "Synonymous":
        SynonymousNumber += 1
    elif data[4] == "Non-synonymous":
        NonSynonymousNumber += 1
    elif data[4] == "Non-coding":
        NonCodingNumber += 1
fig = plt.figure()
variantType = ["Non-synonymous","Synonymous","Non-coding"]
variantNumber = [NonSynonymousNumber,SynonymousNumber,NonCodingNumber]
plt.bar(variantType,variantNumber)
plt.title("2802815 bar plot")
fig.savefig("2802815.png")
plt.clf()

# For tab-seperated table
with open("2802815.tsv","w") as tsvfile:
    tsvfile.write("CHROM  POS REF ALT Type    Transcript  Protein Location  Ref AA  Alt AA \n") # header
    for variant in WholeData:
        row = "\t".join(variant)
        tsvfile.write(row+"\n")

# For log file (txt format) # situation without error
# output files will be store at current directory , so using getcwd to represent the location of thr output files
with open("2802815logfile.txt", "w") as LogFile:
    LogFile.write("Files given at the command line are {} \nThe count of variant where QUAL <= 20 is {} \nThe location of the output files is {} \n".format(str(args.data).replace("[","").replace("]",""),VariantsQualitySmallerThan20,os.getcwd()))