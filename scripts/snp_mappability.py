#############################################################################
# Annotates the snp file from compare_vcfs to include mappability of each 
# site from genmap
#############################################################################

# genmap commands:

# time -p genmap index -F Mus_musculus.GRCm39.dna.primary_assembly.chromes.fa -I Mus_musculus.GRCm39.dna.primary_assembly.chromes-genmap-index
# time -p genmap map -K 150 -E 0 -I Mus_musculus.GRCm39.dna.primary_assembly.chromes-genmap-index/ -t -w -bg -T 12 -O Mus_musculus.GRCm39.dna.primary_assembly.chromes-genmap

#############################################################################

import sys
import os
import gzip
import datetime

#############################################################################

def getDateTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m.%d.%Y | %H:%M:%S");

####################

def PWS(o_line, o_stream=False, std_stream=True, newline=True):
# Function to print a string AND write it to the file.
    if std_stream:
        print(o_line);
    if o_stream:
        o_stream.write(o_line);
        if newline:
            o_stream.write("\n");

####################

def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a string to make it a given length
    spaces = sep * (totlen - len(string));
    return string + spaces;

####################

def runTime(msg=False, writeout=False, printout=True):
# Print out a message at runtime with some info
    if msg:
        if not msg.startswith("#"):
            msg = "# " + msg;
        PWS(msg, writeout, printout);

    PWS("# PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])), writeout, printout)
    PWS("# Script call:    " + " ".join(sys.argv), writeout, printout)
    PWS("# Runtime:        " + datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S"), writeout, printout);
    PWS("# ----------------", writeout, printout);

####################

def fastaReadSeqs(filename, header_sep=False):
# Read a FASTA formatted sequence file
# Great iterator and groupby code from: https://www.biostars.org/p/710/ 
# Returns dictionary with the key:value format as title:sequence.

    from itertools import groupby

    file_stream = open(filename); 
    fa_iter = (x[1] for x in groupby(file_stream, lambda line: line[0] == ">"));
    readstr = lambda s : s.strip();
    # Read the lines of the file
    # file_stream opens the file as an iterable
    # groupby takes an iterable (file_stream) and a function that indicates the key of the group. It iterates over
    # the iterable and when it encounters a key, it groups all following items to it (perfect for FASTA files).
    # fa_iter is a generator object holding all the grouped iterators.
    # readstr is a function that changes depending on compression level 

    seqdict = {};
    # A dictionary of sequences:
    # <sequence id/header> : <sequence>

    for header_obj in fa_iter:
        #print(header_obj)
        header = readstr(header_obj.__next__());
        # The header object is an iterator. This gets the string.

        curkey = header[1:];
        # This removes the ">" character from the header string to act as the key in seqdict

        if header_sep:
            curkey = curkey.split(header_sep)[0];

        seq = "".join(readstr(s) for s in fa_iter.__next__());
        # The current header should correspond to the current iterator in fa_iter. This gets all those
        # lines and combines them as a string.

        #print(header, len(seq));

        seqdict[curkey] = seq;
        # Save the sequence in seqdict

    return seqdict;

#############################################################################

snpfile = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19/summary-files/mm39-30X-vcf-snps.csv.gz";
genmap_txt = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/reference-genomes/mm39/Mus_musculus.GRCm39.dna.primary_assembly.chromes-genmap.txt";
outfile_name = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19/summary-files/mm39-30X-vcf-snps-genmap.csv";
# File names

print("# " + getDateTime() + " Counting lines in SNP file...");
with gzip.open(snpfile, 'rb') as f:
    for i, l in enumerate(f):
        pass
num_lines = str(i + 1);
# Count the number of lines in the SNP file, just for progress updates... only takes a minute...

pad = 20;
with open(outfile_name, "w") as outfile:
    runTime("# Add mappability from genmap to SNPs", outfile);
    PWS(spacedOut("# SNP file:", pad) + snpfile, outfile);
    PWS(spacedOut("# SNP file lines:", pad) + num_lines, outfile);
    PWS(spacedOut("# genmap file:", pad) + genmap_txt, outfile);
    PWS(spacedOut("# Output file:", pad) + outfile_name, outfile);
    PWS("# ----------------", outfile);
    # Some basic info for the summary file

    PWS("# " + getDateTime() + " Reading genmap mappability file...", outfile);
    map_seqs = fastaReadSeqs(genmap_txt);
    for seq in map_seqs:
        map_seqs[seq] = map_seqs[seq].split(" ");
    # Reading the genmap annotation file in FASTA format and splitting the scores into a list

    # print("!", map_seqs["1"][3049999]); # 0
    # print("@", map_seqs["1"][3050000]); # 1
    # print("#", map_seqs["1"][3050001]); # 1
    # print("$", map_seqs["1"][3050445]); # 1
    # print("%", map_seqs["1"][3050446]); # 0.333333
    # print("^", map_seqs["1"][3050447]); # 0.333333
    # print("&", map_seqs["1"][3109465]); # 1
    # print("*", map_seqs["1"][3109466]); # 0.00200803
    # print("(", map_seqs["1"][3109467]); # 0.0020202
    # Some checks

    PWS("# " + getDateTime() + " Annotating SNPs...", outfile);
    first = True;
    i = 0;
    for line in gzip.open(snpfile):        
        line = line.decode();
        # Decode the gzipped line

        if first:
            line = line.replace("\n", ",mappability\n")
            if line.startswith("# "):
                line = line.replace("# ", "");
            outfile.write(line);
            first = False;
            continue;
        # For the first line, add the mappability column to the headers, and remove
        # the # if present

        if line[0] == "#":
            continue;
        # Skip other lines with headers (since this file is catted from multiple)

        line = line.strip().split(",");
        chrome, pos = line[5], int(line[6]);
        line.append(map_seqs[chrome][pos]);
        outfile.write(",".join(line) + "\n");
        # Lookup the mappability of the current SNP and add it to the line and write it

        i += 1;
        if i == 1 or i % 10000 == 0:
            print("# " + getDateTime() + " Annotating SNPs -> " + str(i) + " / " + num_lines);
        # Progress update
    ## End SNP line loop
## Close output file

