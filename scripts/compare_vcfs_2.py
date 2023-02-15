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

def readVCF(vcffile):
# Reads the SNPs in a vcf file

    variant_pos, variants = [], [];
    # The basic info about the SNPs (chr, pos, ref, alt)

    confusing_pos = [];
    # A list of positions that are confusing because NEAT simulated multiple variants at them
    # According to https://github.com/samtools/bcftools/issues/600, only the first variant should
    # be inserted, but I'll just skip them for now
    # Rough count of 4792 out of 1764233 sites

    prev_pos = "";
    x = 1;

    for line in gzip.open(vcffile):
    # Open the file with gzip and read each line
        line = line.decode().strip().split("\t");
        # Read the current line

        if line[0][0] == "#":
            continue;
        # Skip comment lines

        pos = ":".join([ line[0], line[1] ]);
        phase = line[7].replace("WP=", "");
        variant = ":".join([ line[0], line[1], line[3], line[4].replace(",", ";") ]);
        # Parse the basic infor of the SNP and join into a string for set comparisons and use as dict key later

        if pos == prev_pos:
            confusing_pos.append(pos);
            #print(variant);
            x += 1;
            # if x == 100:
            #     sys.exit();
            continue;

        variant_pos.append(pos);
        variants.append(variant);
        prev_pos = pos;
    ## End file loop

    #print(x);

    return variants, set(confusing_pos);
    
#############################################################################

coverage, divergence, heterozygosity, iteration, golden_vcf_file, query_vcf_file, outdir, summary_outfilename, snp_outfilename = sys.argv[1:];
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

if not os.path.isdir(outdir):
    os.makedirs(outdir);

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

    PWS("#" + getDateTime() + " Writing headers for SNP file...", sumfile);
    snp_headers = ["# coverage", "divergence", "heterozygosity", "iteration", "type", "chr", "pos", "ref", "alt", "dp", "mq", "alt.dp", "gq", "pl"];
    snpfile.write(",".join(snp_headers) + "\n");
    # Write the headers for the detailed SNP file

    ####################

    PWS("#" + getDateTime() + " Reading variants from golden file...", sumfile);
    golden_variants, golden_dup_pos = readVCF(golden_vcf_file);
    # Read the SNPs from the golden VCF
    # golden_details is an empty dict since it doesn't have any extra info

    golden_variants = set(golden_variants);
    num_golden = str(len(golden_variants));
    golden_variants_list = list(golden_variants);
    golden_pos = set([ ":".join(v.split(":")[:2]) for v in golden_variants_list ]);
    golden_pos_list = list(golden_pos);
    golden_variants_dict = { golden_pos_list[i] : golden_variants_list[i] for i in range(len(golden_pos_list)) };
    num_golden_pos = str(len(golden_pos));
    PWS("#" + getDateTime() + " " + num_golden + " variants read at " + num_golden_pos + " positions.", sumfile);
    #print(len(golden_variants_dict));
    # Convert the list of SNPs to a set and count
    # Some variants at same site but other haplotype

    ####################

    PWS("#" + getDateTime() + " Reading variants from query file and checking overlaps...", sumfile);
    num_query, num_golden, num_golden_dup_pos, num_non_golden, num_tp, num_fp, num_fn, num_tn, no_info = 0,0,0,0,0,0,0,0,0;
    num_het_in_golden = 0;
    num_het_not_in_golden = 0;
    # The counts of each type of site category

    for line in gzip.open(query_vcf_file):
    # Open the file with gzip and read each line

        line = line.decode().strip().split("\t");
        # Read the current line

        if line[0][0] == "#":
            continue;
        # Skip comment lines

        if line[1] in golden_dup_pos:
            num_golden_dup_pos += 1;
            continue;
        # Skip the positions that had more than 1 variant simulated, for now

        query_variant = False;
        query_gt = False;
        golden_variant = False;
        # Flags for the current line to determine site overlaps

        cur_details = { "dp" : "NA", "mq": "NA", "alt.dp" : "NA", "gq" : "NA", "pl" : "NA" };
        # Initial values for details if there is a SNP at either of the golden or query sites

        if line[4] != ".":
        # If the current query site is a variant
            query_variant = True;
            num_query += 1;
            # Set the query_variant flag to True

            variant = ":".join([ line[0], line[1], line[3], line[4].replace(",", ";") ]);
            query_gt = line[9].split(":")[0];        
            # Parse the basic info of the SNP and join into a string for set comparisons and use as dict key later
        ## End query variant block

        pos = ":".join([ line[0], line[1] ]);
        if pos in golden_pos:
            golden_variant = True;
            num_golden += 1;
        else:
            num_non_golden += 1;
        # Parse the current position and set golden_variant to True if the position is in golden_variants

        type_str = "NA";
        # Categories of variant overlaps: tp, tp.pos, fp, fn
        # tp.pos is when a variant was called at the same site as a golden variant, but with a different alternate allele

        if query_variant and golden_variant:
        # If the current site is variant in both sets it will be tp or tp.pos
            if query_gt in ["0/1", "1/0", "1|0", "0|1"]:
                num_het_in_golden += 1;
                type_str = "het.golden";
            elif query_gt in ["1/1", "1|1"]:
                num_tp += 1;
                type_str = "tp";

            # if variant in golden_variants:
            #     num_tp += 1;
            #     type_str = "tp";
            # # If the variant is identical to one in the golden variants, it is a true positive
            # # Increment tp count and set type_str to tp
            # else:
            #     num_tp_pos += 1;
            #     type_str = "tp.pos";
            # # If the variant doesn't exist in the golden variant set, but the position does overlap, that
            # # means it is a called variant with a different alternate allele, which I call tp.pos for true positive
            # # at position
        ## End tp block

        if query_variant and not golden_variant:
        # If the current site is variant in the query but not the golden set, it is a false positive
            if query_gt in ["0/1", "1/0", "1|0", "0|1"]:
                num_het_not_in_golden += 1;
                type_str = "het.non.golden";
            elif query_gt in ["1/1", "1|1"]:
                num_fp += 1;
                type_str = "fp";
        ## End fp block

        if not query_variant and golden_variant:
        # If the current site is variant in the golden set but not the query, it is a false negative
            num_fn += 1;
            type_str = "fn";

            variant = ":".join([ line[0], line[1], line[3], line[4].replace(",", ";") ]);
            # Also get the basic site info for a fn to write the details later
        ## End fn block

        if not query_variant and not golden_variant:
            num_tn += 1;
            type_str = "tn";
        ## End tn block

        if query_variant or golden_variant:
        # If the site is variant in either set, we want to get some more details

            fmt = line[9].split(":");
            info = line[7].split(";");
            # Splot the format and info lines

            if info == ['.']:
                no_info += 1;
            # For some reason at some sites no info is found... keep track but skip

            else:         
                dp, mq = "NA", "NA";
                # Set some default values for dp and mq so we can track how many sites have problems later

                for entry in info:
                    if "DP=" in entry:
                        dp = entry.replace("DP=", "");
                    if "MQ=" in entry:
                        mq = entry.replace("MQ=", "");
                # Parse dp and mq from the info field

                if query_variant:
                # We can get more info for called variants

                    pl_ind = 4
                    if "PGT" in line[8]:
                        pl_ind = 6;
                    # The index of pl changes depending on whether the variant is phased

                    cur_details = { "dp" : dp, "mq": mq, "alt.dp" : fmt[1].split(",")[1], "gq" : fmt[3], "pl" : fmt[pl_ind].replace(",", ";") };
                    # Save details to dict to write later

                else:
                    cur_details = { "dp" : dp, "mq": mq, "alt.dp" : "NA", "gq" : "NA", "pl" : "NA" };
                    # Save details to dict to write later
            ## End detail parse block

            outline = [coverage, divergence, heterozygosity, iteration, type_str] + variant.split(":") + [ cur_details[val] for val in snp_headers[9:] ];
            snpfile.write(",".join(outline) + "\n");
            # Write the details to the SNP file
        ## End detailed block
    ## End file loop

    # query_variants, query_non_variants, query_details = readVCF(query_vcf_file, detailed=True);
    # print(len(query_details));
    # # Read the SNPs from the query VCF
    # # query_details is retrieved since details is set to True
    
    # query_variants = set(query_variants);
    # print(bool("19:49023442:T:C" in query_variants));
    # print(bool("19:49023442:T:C" in golden_variants));
    # num_query = str(len(query_variants));
    # PWS("# " + num_query + " variants read", sumfile);
    # # Convert the list of SNPs to a set and count

    ####################

    # PWS("# Counting true positives...", sumfile);
    # tp = set.intersection(golden_variants, query_variants);
    # num_tp = str(len(tp));
    # #print(num_tp);
    # #print(list(tp)[:100])
    # # True positives are the intersect between the golden and query SNP sets

    # for variant in tp:
    #     #print(variant);
    #     #print(query_details[variant]);
    #     outline = [coverage, divergence, "tp"] + variant.split(":") + [ query_details[variant][val] for val in snp_headers[7:] ];
    #     snpfile.write(",".join(outline) + "\n");
    # # Get the details of the query SNPs classified as true positives and write them to the SNP file

    # ####################

    # PWS("# Counting false positives...", sumfile);
    # fp = query_variants - golden_variants
    # num_fp = str(len(fp));
    # # False positives are SNPs are those found in the query but not in the golden set

    # for variant in fp:
    #     outline = [coverage, divergence, "fp"] + variant.split(":") + [ query_details[variant][val] for val in snp_headers[7:] ];
    #     snpfile.write(",".join(outline) + "\n");
    # # Get the details of the query SNPs classified as false positives and write them to the SNP file

    # fp_pos = list(fp);
    # fp_pos = [ ":".join(v.split(":")[:2]) for v in fp_pos ];
    # fp_pos = set(fp_pos);
    # # This creates a list of false positive SNPs with only their location in the genome
    # # to check for SNPs that are called erroneously, but at a correct position

    # ####################

    # PWS("# Counting false negatives...", sumfile);
    # fn = golden_variants - query_variants;
    # num_fn = str(len(fn));
    # #print(num_fn);
    # # False negatives are SNPs are those found in the golden but not in the query set

    # for variant in fn:
    #     outline = [coverage, divergence, "fn"] + variant.split(":") + [ query_details[variant][val] for val in snp_headers[7:] ];
    #     snpfile.write(",".join(outline) + "\n");
    # # Get the details of the query SNPs classified as false positives and write them to the SNP file

    # fn_pos = list(fn);
    # fn_pos = [ ":".join(v.split(":")[:2]) for v in fn_pos ];
    # fn_pos = set(fn_pos);
    # # This creates a list of false negative SNPs with only their location in the genome
    # # to check for SNPs that are called erroneously, but at a correct position

    # ####################

    # tp_pos = fp_pos & fn_pos;
    # num_tp_pos = str(len(tp_pos));
    # #print(num_tp_pos);
    # # SNPs that are in both the false positive and false negative sets BY POSITION are those where a SNP was called
    # # at a correct location, but with the wrong alternate allele, e.g. 19:1:A>T vs. 19:1:A>G
    # # Get the union of the position sets to count these true positives by position

    ####################

    headers = ["region", "coverage", "divergence", "heterozygosity", "iteration", "golden variants", "golden variants same site", "golden invariant", "called variants", "tp", "fp", "fn", "tn", "called het golden", "called het non golden", "no info"];
    PWS(",".join(headers), sumfile);
    # Write the headers for the summary file

    outline = [region, coverage, divergence, heterozygosity, iteration, str(num_golden), str(num_golden_dup_pos), str(num_non_golden), str(num_query), str(num_tp), str(num_fp), str(num_fn), str(num_tn), str(num_het_in_golden), str(num_het_not_in_golden), str(no_info)];
    PWS(",".join(outline), sumfile, newline=False);
    # Write the final counts to the summary file
## End file block and close files

