#############################################################################
# Compares variants called between two VCF files
#############################################################################

import sys
import os
import gzip
import datetime
from collections import Counter

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

def readVCF(vcffile, detailed=False):
# Reads the SNPs in a vcf file

    variants, non_variants = [], [];
    # The basic info about the SNPs (chr, pos, ref, alt)

    details = {};
    # Details about the SNPs (dp, ad, gq, pl)

    for line in gzip.open(vcffile):
    # Open the file with gzip and read each line
        line = line.decode().strip().split("\t");
        # Read the current line

        if line[0][0] == "#":
            continue;
        # Skip comment lines

        variant = ":".join([ line[0], line[1], line[3], line[4] ]);
        # Parse the basic infor of the SNP and join into a string for set comparisons and use as dict key later

        if detailed:
        # For the query VCF, we get some more detailed info to output per SNP
            fmt = line[9].split(":");
            info = line[7].split(";");

            if info == ['.']:
                cur_details = { "dp" : "NA", "mq": "NA", "alt.dp" : "NA", "gq" : "NA", "pl" : "NA" };

            else:
                dp, mq = "NA", "NA";

                for entry in info:
                    if "DP=" in entry:
                        dp = entry.replace("DP=", "");
                    if "MQ=" in entry:
                        mq = entry.replace("MQ=", "");

                if line[4] != ".":
                    #print(line);
                    pl_ind = 4
                    if "PGT" in line[9]:
                        pl_ind = 6;
                    cur_details = { "dp" : dp, "mq": mq, "alt.dp" : fmt[1].split(",")[1], "gq" : fmt[3], "pl" : fmt[pl_ind] };
                else:
                    cur_details = { "dp" : dp, "mq": mq, "alt.dp" : "NA", "gq" : "NA", "pl" : "NA" };
                details[variant] = cur_details;


            # Save the details for the current SNP in the details dict with the SNP string as the key
        ## End detailed block

        if line[4] == ".":
            non_variants.append(variant);
        else:
            variants.append(variant);
    ## End file loop

    return variants, non_variants, details;
    
#############################################################################

coverage, divergence, golden_vcf_file, query_vcf_file, summary_outfilename, snp_outfilename = sys.argv[1:];
# Inputs

#coverage = "20";
#divergence = "0.02";
#golden_vcf_file = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39/simulated-reads/20X/0.02/regions/mm39-19_golden.vcf.gz";
#query_vcf_file = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39/called-variants/gatk/20X/0.02/regions/mm39-19-20X-0.02.vcf.gz";
#summary_outfilename = "test-summary.csv";
#snp_outfilename = "test-snps.csv";
# Test inputs

region = os.path.basename(query_vcf_file).split("-")[1];
# The region as a string from the VCF file name

####################

pad = 20;
with open(summary_outfilename, "w") as sumfile, open(snp_outfilename, "w") as snpfile:
# Open summary file and output file for writing

    runTime("# Compare VCFs", sumfile);
    PWS(spacedOut("# Golden VCF:", pad) + golden_vcf_file, sumfile);
    PWS(spacedOut("# Query VCF:", pad) + query_vcf_file, sumfile);
    PWS(spacedOut("# Summary file:", pad) + summary_outfilename, sumfile);
    PWS(spacedOut("# Output file:", pad) + snp_outfilename, sumfile);
    PWS("# ----------------", sumfile);
    # Some basic info for the summary file

    ####################

    PWS("# Writing headers for SNP file...", sumfile);
    snp_headers = ["divergence", "coverage", "type", "chr", "pos", "ref", "alt", "dp", "mq", "alt.dp", "gq", "pl"];
    snpfile.write(",".join(snp_headers) + "\n");
    # Write the headers for the detailed SNP file

    ####################

    PWS("# Reading variants from golden file...", sumfile);
    golden_variants, golden_non_variants, golden_details = readVCF(golden_vcf_file);
    # Read the SNPs from the golden VCF
    # golden_details is an empty dict since it doesn't have any extra info

    golden_variants = set(golden_variants);
    num_golden = str(len(golden_variants));
    PWS("# " + num_golden + " variants read", sumfile);
    # Convert the list of SNPs to a set and count

    ####################

    PWS("# Reading variants from query file...", sumfile);
    query_variants, query_non_variants, query_details = readVCF(query_vcf_file, detailed=True);
    print(len(query_details));
    # Read the SNPs from the query VCF
    # query_details is retrieved since details is set to True
    
    query_variants = set(query_variants);
    print(bool("19:49023442:T:C" in query_variants));
    print(bool("19:49023442:T:C" in golden_variants));
    num_query = str(len(query_variants));
    PWS("# " + num_query + " variants read", sumfile);
    # Convert the list of SNPs to a set and count

    ####################

    PWS("# Counting true positives...", sumfile);
    tp = set.intersection(golden_variants, query_variants);
    num_tp = str(len(tp));
    #print(num_tp);
    #print(list(tp)[:100])
    # True positives are the intersect between the golden and query SNP sets

    for variant in tp:
        #print(variant);
        #print(query_details[variant]);
        outline = [coverage, divergence, "tp"] + variant.split(":") + [ query_details[variant][val] for val in snp_headers[7:] ];
        snpfile.write(",".join(outline) + "\n");
    # Get the details of the query SNPs classified as true positives and write them to the SNP file

    ####################

    PWS("# Counting false positives...", sumfile);
    fp = query_variants - golden_variants
    num_fp = str(len(fp));
    # False positives are SNPs are those found in the query but not in the golden set

    for variant in fp:
        outline = [coverage, divergence, "fp"] + variant.split(":") + [ query_details[variant][val] for val in snp_headers[7:] ];
        snpfile.write(",".join(outline) + "\n");
    # Get the details of the query SNPs classified as false positives and write them to the SNP file

    fp_pos = list(fp);
    fp_pos = [ ":".join(v.split(":")[:1]) for v in fp_pos ];
    fp_pos = set(fp_pos);
    # This creates a list of false positive SNPs with only their location in the genome
    # to check for SNPs that are called erroneously, but at a correct position

    ####################

    PWS("# Counting false negatives...", sumfile);
    fn = golden_variants - query_variants;
    num_fn = str(len(fn));
    #print(num_fn);
    # False negatives are SNPs are those found in the golden but not in the query set

    for variant in fn:
        outline = [coverage, divergence, "fn"] + variant.split(":") + [ query_details[variant][val] for val in snp_headers[7:] ];
        snpfile.write(",".join(outline) + "\n");
    # Get the details of the query SNPs classified as false positives and write them to the SNP file

    fn_pos = list(fn);
    fn_pos = [ ":".join(v.split(":")[:1]) for v in fn_pos ];
    fn_pos = set(fn_pos);
    # This creates a list of false negative SNPs with only their location in the genome
    # to check for SNPs that are called erroneously, but at a correct position

    ####################

    tp_pos = fp_pos & fn_pos;
    num_tp_pos = str(len(tp_pos));
    #print(num_tp_pos);
    # SNPs that are in both the false positive and false negative sets BY POSITION are those where a SNP was called
    # at a correct location, but with the wrong alternate allele, e.g. 19:1:A>T vs. 19:1:A>G
    # Get the union of the position sets to count these true positives by position

    ####################

    headers = ["region", "coverage", "divergence", "golden variants", "called variants", "tp", "fp", "fn", "shared pos"];
    PWS(",".join(headers), sumfile);
    # Write the headers for the summary file

    outline = [region, coverage, divergence, num_golden, num_query, num_tp, num_fp, num_fn, num_tp_pos];
    PWS(",".join(outline), sumfile, newline=False);
    # Write the final counts to the summary file
## End file block and close files

