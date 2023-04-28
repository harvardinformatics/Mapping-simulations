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
    return datetime.datetime.now().strftime(" %m.%d.%Y | %H:%M:%S");

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

def detectCompression(filename):
# Detect compression of a file by examining the first lines in the file

    compression_type = "none";

    magic_dict = {
            b"\x1f\x8b\x08": "gz",
            # b"\x1f\x8b\x08\x08": "gz",
            b"\x42\x5a\x68": "bz2",
            b"\x50\x4b\x03\x04": "zip"
        }
    # An encoded set of possible "magic strings" that start different types of compressed files
    # From: https://www.garykessler.net/library/file_sigs.html
    # \x is the escape code for hex values
    # b converts strings to bytes

    max_len = max(len(x) for x in magic_dict)
    # The number of characters to read from the beginning of the file should be the length of
    # the longest magic string

    file_start = open(filename, "rb").read(max_len);
    # Read the beginning of the file up to the length of the longest magic string

    for magic_string in magic_dict:
        if file_start.startswith(magic_string):
            compression_type = magic_dict[magic_string];
    # Check each magic string against the start of the file

    return compression_type;

####################

def fastaReadSeqs(filename, regions=False, header_sep=False):
# Read a FASTA formatted sequence file
# Great iterator and groupby code from: https://www.biostars.org/p/710/ 
# Returns dictionary with the key:value format as title:sequence.

    from itertools import groupby

    compression = detectCompression(filename);

    if compression == "gz":
        file_stream = gzip.open(filename);
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line.decode()[0] == ">"));
        readstr = lambda s : s.decode().strip();
    elif compression == "none":
        file_stream = open(filename); 
        fa_iter = (x[1] for x in groupby(file_stream, lambda line: line[0] == ">"));
        readstr = lambda s : s.strip();
    # Read the lines of the file depending on the compression level
    # file_stream opens the file as an iterable
    # groupby takes an iterable (file_stream) and a function that indicates the key of the group. It iterates over
    # the iterable and when it encounters a key, it groups all following items to it (perfect for FASTA files).
    # fa_iter is a generator object holding all the grouped iterators.
    # readstr is a function that changes depending on compression level -- for compressed files we also need to decode
    # each string in the iterators below.

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

        if regions and curkey not in regions:
            continue;

        seq = "".join(readstr(s) for s in fa_iter.__next__());
        # The current header should correspond to the current iterator in fa_iter. This gets all those
        # lines and combines them as a string.

        #print(header, len(seq));

        seqdict[curkey] = seq;
        # Save the sequence in seqdict

    return seqdict;

####################

def readVCF(vcffile, regions=False, filter_str="PASS", minimap=False, debug=False):
# Reads the SNPs in a vcf file

    compression = detectCompression(vcffile);

    if compression == "gz":
        opener = gzip.open;
        reader = lambda s : s.decode().strip().split("\t");
    elif compression == "none":
        opener = open; 
        reader = lambda s : s.strip().split("\t");

    variants, variant_list, extra = {}, [], [];
    # The basic info about the SNPs (chr, pos, ref, alt)

    confusing_pos = [];
    # A list of positions that are confusing because NEAT simulated multiple variants at them
    # According to https://github.com/samtools/bcftools/issues/600, only the first variant should
    # be inserted, but I'll just skip them for now
    # Rough count of 4792 out of 1764233 sites

    prev_pos = "";
    x = 1;

    for line in opener(vcffile):
    # Open the file with gzip and read each line
        line = reader(line);
        # Read the current line

        if line[0][0] == "#":
            continue;
        # Skip comment lines

        if line[6] != filter_str:
            continue;
        # Skip lines with variants that are filtered        

        if regions and line[0] not in regions:
            continue;

        pos = ":".join([ line[0], line[1] ]);

        if minimap:
            coord = line[1];
            qstart = line[7].split(";")[1].replace("QSTART=", "");
            if coord != qstart:
                extra += [coord, qstart];
                continue;
        # This seems to happen when minimap predicts an indel, but in some cases just a SNP?
        # Happened ~40 times in my test case so just skip for now

        alt = line[4];
        if "*" in alt:
            continue;
        gt = line[-1].split(":")[0].replace("WP=", "");

        if pos == prev_pos:
            confusing_pos.append(pos);
            x += 1;
            continue;

        variant_list.append(pos);
        variants[pos] = {'region' : line[0], 'pos' : line[1], 'alt' : alt, 'gt' : gt };
        prev_pos = pos;
    ## End file loop

    #print(x);

    return variants, variant_list, set(confusing_pos), extra;
    
#############################################################################

region, coverage, divergence, heterozygosity, iteration, genome_file, iter_file, prev_iter_file, golden_div_vcf, iteration_vcf, minimap_vcf, outdir, outfilename = sys.argv[1:];
# Inputs

# coverage = "30";
# divergence = "0.02";
# heterozygosity = "0.005";
# iteration = "2";
# region = "19";
# genome_file = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/reference-genomes/mm39/Mus_musculus.GRCm39.dna.primary_assembly.fa";
# iter_file = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19-iterative/consensus/gatk/30X/0.02/iter2/mm39-30X-0.02d-0.005h-snps-consensus.fa"
# prev_iter_file = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19-iterative/consensus/gatk/30X/0.02/iter1/mm39-30X-0.02d-0.005h-snps-consensus.fa"
# golden_div_vcf = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19/simulated-reads/30X/0.02/divergent/regions/mm39-19_golden.vcf.gz";
# #golden_het_vcf = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19-iterative/simulated-reads/30X/0.02/heterozygous/0.005/regions/mm39-19_golden.vcf.gz";
# iteration_vcf = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19-iterative/called-variants/gatk/30X/0.02/iter2/mm39-30X-0.02d-0.005h-filtered-snps.vcf.gz";
# minimap_vcf = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19-iterative/consensus/gatk/30X/0.02/iter2/mm39-30X-0.02d-0.005h-snps-consensus-to-ref.vcf";
# outdir = ".";
# outfilename = "test.tsv";
# Test inputs

#region = os.path.basename(golden_div_vcf).split("-")[1];
# The region as a string from the VCF file name

# if not os.path.isdir(outdir):
#     os.makedirs(outdir);

####################

pad = 30;
with open(outfilename, "w") as outfile:
# Open summary file and output file for writing

    runTime("# Compare VCFs", outfile);
    PWS(spacedOut("# Reference:", pad) + genome_file, outfile);
    PWS(spacedOut("# Previous iteration:", pad) + prev_iter_file, outfile);
    PWS(spacedOut("# Region:", pad) + region, outfile);
    PWS(spacedOut("# Golden div VCF:", pad) + golden_div_vcf, outfile);
    #PWS(spacedOut("# Golden het VCF:", pad) + golden_het_vcf, outfile);
    PWS(spacedOut("# Query iteration VCF:", pad) + iteration_vcf, outfile);
    PWS(spacedOut("# Query minimap VCF:", pad) + minimap_vcf, outfile);
    PWS(spacedOut("# Output file:", pad) + outfilename, outfile);
    PWS("# ----------------", outfile);
    # Some basic info for the summary file

    ####################

    PWS("#" + getDateTime() + " Reading original reference sequence...", outfile);
    ref_seq = fastaReadSeqs(genome_file, regions=[region], header_sep=" ");
    PWS("#" + getDateTime() + " " + str(len(ref_seq)) + " region(s) read", outfile);
    # Read ref seq

    ####################

    PWS("#" + getDateTime() + " Reading iteration sequence...", outfile);
    iter_seq = fastaReadSeqs(iter_file, regions=[region], header_sep=" ");
    PWS("#" + getDateTime() + " " + str(len(iter_seq)) + " region(s) read", outfile);
    # Read ref seq

    ####################

    PWS("#" + getDateTime() + " Reading previous iteration sequence...", outfile);
    prev_seq = fastaReadSeqs(prev_iter_file, regions=[region], header_sep=" ");
    PWS("#" + getDateTime() + " " + str(len(prev_seq)) + " region(s) read", outfile);
    # Read ref seq

    ####################

    PWS("#" + getDateTime() + " Reading variants from golden div file...", outfile);
    golden_div_variants, golden_div_variant_list, golden_div_dup_pos, extra = readVCF(golden_div_vcf);
    # Read the SNPs from the golden VCF

    # print(len(golden_div_variants));
    # print(len(golden_div_variant_list));
    # print(len(list(set(golden_div_variant_list))));
    # print(len(golden_div_dup_pos));
    # sys.exit();

    PWS("#" + getDateTime() + " " + str(len(golden_div_variants)) + " variants read", outfile);

    ####################

    PWS("#" + getDateTime() + " Reading variants from iteration file...", outfile);
    iter_div_variants, iter_div_variant_list, iter_div_dup_pos, extra = readVCF(iteration_vcf, regions=[region]);
    # Read the SNPs from the iteration VCF

    # print(len(iter_div_variants));
    # print(len(iter_div_variant_list));
    # print(len(list(set(iter_div_variant_list))));
    # print(len(iter_div_dup_pos));

    PWS("#" + getDateTime() + " " + str(len(iter_div_variants)) + " variants read", outfile);

    ####################

    PWS("#" + getDateTime() + " Reading variants from minimap file...", outfile);
    mmap_div_variants, mmap_div_variant_list, mmap_div_dup_pos, mmap_qstart_mismatches = readVCF(minimap_vcf, regions=[region], filter_str=".", minimap=True);
    # Read the SNPs from the iteration VCF

    # print(len(mmap_div_variants));
    # print(len(mmap_div_variant_list));
    # print(len(list(set(mmap_div_variant_list))));
    # print(len(mmap_div_dup_pos));
    # sys.exit();

    PWS("#" + getDateTime() + " " + str(len(mmap_div_variants)) + " variants read", outfile);
    PWS("#" + getDateTime() + " " + str(len(mmap_qstart_mismatches)/2) + " mis-matched coordinates will be ignored", outfile);

    ####################

    PWS("#" + getDateTime() + " Comparing variants...", outfile);
    meta_headers = [ "coverage", "divergence", "heterozygosity", "iteration" ];
    headers = [ "chr", "pos", "ref", "golden.div", "prev.iter", "iter.gt", "iter.allele.1", "iter.allele.2", "minimap.gt", "minimap.allele.1", "minimap.allele.2" ];
    outfile.write("\t".join(meta_headers + headers) + "\n");
    # Write the headers for the detailed SNP file

    ####################

    meta_outline = [ coverage, divergence, heterozygosity, iteration ]
    
    for i in range(len(ref_seq[region])):
        j = i + 1;
        pos = region + ":" + str(j);

        # prev_match = "Y";
        # if ref_seq[region][i].upper() != prev_seq[region][i].upper():
        #     prev_match = "N"

        if any(pos in varset for varset in [golden_div_variants, iter_div_variants, mmap_div_variants] ):
            outline = { "chr" : region, "pos" : str(i), "ref" : ref_seq[region][i], 'golden.div' : ".",
                        "prev.iter" : prev_seq[region][i], "iter.gt" : "NA", "iter.allele.1" : ".", "iter.allele.2" : ".", 
                        "minimap.gt" : "NA", "minimap.allele.1" : ".", "minimap.allele.2" : "." };

            if pos in golden_div_variants:
                outline['golden.div'] = golden_div_variants[pos]['alt'];
            else:
                outline['golden.div'] = ref_seq[region][i];

            if pos in iter_div_variants:
                outline['iter.allele.1'] = iter_div_variants[pos]['alt'];

                if iter_div_variants[pos]['gt'] in ["1/1", "1|1"]:
                    outline['iter.allele.2'] = iter_div_variants[pos]['alt'];
                    outline['iter.gt'] = "hom.alt";
                else:
                    outline['iter.allele.2'] = prev_seq[region][i];
                    outline['iter.gt'] = "het";
            else:
                outline['iter.allele.1'] = prev_seq[region][i];
                outline['iter.allele.2'] = prev_seq[region][i];
                outline['iter.gt'] = "hom.prev";

                # if pos == "19:3053895":
                #     print(iter_seq[region][i], pos + " " + prev_seq[region][i] + " " + iter_seq[region][i]);

                #if "*" not in iter_seq[region][i] and "*" not in prev_seq[region][i]:
                assert prev_seq[region][i] == iter_seq[region][i], pos + " " + prev_seq[region][i] + " " + iter_seq[region][i];

            if pos in mmap_div_variants:
                outline['minimap.allele.1'] = mmap_div_variants[pos]['alt'];

                if mmap_div_variants[pos]['gt'] in ["1/1", "1|1"]:
                    outline['minimap.allele.2'] = mmap_div_variants[pos]['alt'];
                    outline['minimap.gt'] = "hom.alt";
                else:
                    outline['minimap.allele.2'] = ref_seq[region][i];
                    outline['minimap.gt'] = "het";
            else:
                outline['minimap.allele.1'] = ref_seq[region][i];
                outline['minimap.allele.2'] = ref_seq[region][i];
                outline['minimap.gt'] = "hom.ref";

                # if str(j) not in mmap_qstart_mismatches:
                #     assert ref_seq[region][i] == iter_seq[region][i], str(i) + " " + str(j) + " " + pos + " " + ref_seq[region][i] + " " + iter_seq[region][i];
                # 19:3090394

            final_outline = "\t".join(meta_outline + [ outline[col] for col in headers ]);
            outfile.write(final_outline + "\n");
    ## End site loop
## End file block and close files

