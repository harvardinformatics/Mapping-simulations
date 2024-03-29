#############################################################################
# Compares reads mapped between two BAM files
#############################################################################

import sys
import os
import pysam
import gzip
import datetime
import random
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
# Function to print out some info at the start of the script
    if msg:
        if not msg.startswith("#"):
            msg = "# " + msg;
        PWS(msg, writeout, printout);

    PWS("# PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])), writeout, printout)
    PWS("# Script call:    " + " ".join(sys.argv), writeout, printout)
    PWS("# Runtime:        " + datetime.datetime.now().strftime("%m/%d/%Y %H:%M:%S"), writeout, printout);
    PWS("# ----------------", writeout, printout);

#############################################################################

regions, coverage, divergence, heterozygosity, iteration, ref_index, golden_bam_file, query_bam_file, summary_outfilename, exact_outfilename, unmapped_outfilename = sys.argv[1:];
regions = regions.split(",");
# Inputs

# divergence = "0.02";
# heterozygosity = "0.005";
# iteration = "1";
# coverage = "30";
# regions = ["18"];
# golden_bam_file = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-iterative/simulated-reads/30X/0.02/heterozygous/0.005/mm39_golden.bam";
# query_bam_file = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-iterative/mapped-reads/30X/0.02/0.005/iter2/mm39-30X-0.02d-0.005h-crossmap.sorted.bam";
# ref_index = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/reference-genomes/mm39/Mus_musculus.GRCm39.dna.primary_assembly.chromes.fa.fai";
# summary_outfilename = "test.csv";
# exact_outfilename = "test-exact.bam";
# unmapped_outfilename = "test-unmapped.bam";
#tracefilename = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19/summary-files/20X/0.02/mm39-20X-0.02-compare-bam-trace.txt";
# Test inputs

pad = 20;
with open(summary_outfilename, "w") as outfile:#, open(tracefilename, "w") as tracefile:
    runTime("# Compare BAMs", outfile);
    PWS(spacedOut("# Golden BAM:", pad) + golden_bam_file, outfile);
    PWS(spacedOut("# Query BAM:", pad) + query_bam_file, outfile);
    PWS(spacedOut("# Output file:", pad) + summary_outfilename, outfile);
    PWS("# ----------------", outfile);
    # Runtime info

    ####################

    #window_size = 10000
    PWS("# " + getDateTime() + " Reading reference index...", outfile);
    chr_lens = { line.split("\t")[0] : int(line.split("\t")[1]) for line in open(ref_index, "r") };
    # Reads the reference index

    ####################

    PWS("# " + getDateTime() + " Opening BAM files...", outfile);
    golden_bam = pysam.AlignmentFile(golden_bam_file, "rb");
    query_bam = pysam.AlignmentFile(query_bam_file, "rb");
    # Opens the BAM files

    ####################

    PWS("# " + getDateTime() + " Indexing query BAM by name...", outfile);
    name_indexed_query = pysam.IndexedReads(query_bam);
    name_indexed_query.build();
    # Index the query BAM by name so we can look up reads easily below
    # https://timoast.github.io/blog/2015-10-12-extractreads/
    # Takes about half an hour and uses almost 50GB of RAM on 40GB test BAM

    ####################

    outdict = { "has.match" : 0, "reads.w.primary.only" : 0, "reads.w.secondary.only" : 0, "reads.w.both" : 0, "primary.err" : 0,
                 "exact.map" : 0, "close.map" : 0, "mismapped" : 0, "diff.chr" : 0, "unmapped" : 0, "multi.match" : 0, "pair.err" : 0, "pos.err" : 0, "chr.err" : 0, "missing" : 0 };
    # Just doing counts for now. Could do full distributions of distance between golden and mapped positions, but
    # would have to deal with the memory issues

    PWS("# " + getDateTime() + " Finding reads...", outfile);

    with pysam.AlignmentFile(unmapped_outfilename, "wb", header=golden_bam.header) as unmapped_out, pysam.AlignmentFile(exact_outfilename, "wb", header=golden_bam.header) as exact_out:
        for chrome in chr_lens:
        # Go through the reads by chromosome

            if chrome not in regions:
                continue;
            # For debugging/testing

            PWS("# " + getDateTime() + " Finding reads in " + chrome + "...", outfile);

            num_reads = 0;
            # Counter for progress report

            for read in golden_bam.fetch(chrome):
            # Find every golden read in the current chromosome

                if num_reads % 1000 == 0:
                    print(getDateTime() + " " + str(num_reads));
                num_reads += 1;
                # Progress report

                if read.is_secondary:
                    continue;
                # Only concerned with primary mappings

                if read.is_read1:
                    read_pair = "1";
                else:
                    read_pair = "2";
                # Since NEAT doesn't make it easy to determine which read is which in a pair, get it here

                golden_chr = read.reference_name;
                golden_pos = read.get_reference_positions()[0];
                # Get the chromosome and position where this read should map

                try:
                    matches = name_indexed_query.find(read.query_name);
                except:
                    outdict['missing'] += 1;
                    continue;
                # Find all mappings for the current read in the query bam file

                matches_found = 0;

                has_primary = False;
                has_secondary = False;

                pair_err = False;
                pos_err = False;
                chr_err = False;
                # Keeps track of how many matches we categorize - should only be one primary match

                if matches:
                    outdict['has.match'] += 1;

                for match in matches:
                # Go through every matching read in the query file for the current read

                    if (match.is_read1 and read_pair == "1") or (match.is_read2 and read_pair == "2"):
                    # Make sure they are the same read in the pair

                        classification = "NONE";
                        pair_err = False;

                        if match.is_secondary:
                            has_secondary = True;
                            continue;
                        else:
                            has_primary = True;
                        # Only concerned with primary mappings

                        query_chr = match.reference_name;
                        query_pos_list = match.get_reference_positions();
                        # Get the info for the matching query

                        if golden_chr != query_chr:
                            chr_err = False;
                            classification = "DIFF CHR";
                            outdict['diff.chr'] += 1;
                            matches_found += 1;
                        # If it mapped to a different chromosome completely

                        elif golden_chr == query_chr:
                            chr_err = False;

                            if query_pos_list:
                                pos_err = False;

                                query_pos = query_pos_list[0];
                                # Get the chromosome and position of where the read actually mapped

                                if golden_pos == query_pos:
                                    classification = "EXACT MAP"
                                    outdict['exact.map'] += 1;
                                    matches_found += 1;
                                    exact_out.write(read);
                                # An exact match

                                elif abs(query_pos - golden_pos) <= 150:
                                    classification = "CLOSE MAP"
                                    outdict['close.map'] += 1;
                                    matches_found += 1;
                                # A mapping within one read length of the true mapping

                                else:
                                    classification = "MISMAPPED"
                                    outdict['mismapped'] += 1;
                                    matches_found += 1;
                                # A mapping on the same chromosome, but a different position

                            elif matches_found == 0:
                                pos_err = True;
                            ## End position check block

                        elif matches_found == 0:
                            chr_err = True;
                        ## End chrome check block

                    else:
                        continue;
                    ## End pair check block
                ## End matches loop

                if has_secondary and has_primary:
                    outdict['reads.w.both'] += 1;
                elif has_secondary:
                    print(read.query_name);
                    mate = golden_bam.mate(read);
                    print(mate.query_name);
                    sys.exit();
                    outdict['reads.w.secondary.only'] += 1;
                elif has_primary:
                    outdict['reads.w.primary.only'] += 1;
                else:
                    outdict['primary.err'] += 1;

                if matches_found != 1:
                    if matches_found > 1:
                        outdict['multi.match'] += 1;
                    # Check if multiple matches were found

                    elif matches_found == 0:
                        if pair_err:
                            outdict['pair.err'] += 1;
                        elif pos_err:
                            outdict['pos.err'] += 1;
                        
                        elif chr_err:
                            outdict['chr.err'] += 1;
                        else:
                            outdict['unmapped'] += 1;
                            unmapped_out.write(read);
                    # Check various cases where no match was found
                ## End matches loop
            ## End golden bam read loop
        ## End chrome loop

        headers = ["coverage", "divergence", "heterozygosity", "iteration", 
                    "total.reads", "has.match", "multi.match", 
                    "primary.match.only", "secondary.match.only", "both.matches", "no.matches",
                    "missing.in.query", "primary.err", "pair.err", "pos.err", "chr.err", 
                    "exact.map", "close.map", "mismapped", "diff.chr", "unmapped"];
        PWS(",".join(headers), outfile);

        outline = [coverage, divergence, heterozygosity, iteration, 
                    num_reads, outdict['has.match'], outdict['multi.match'],
                    outdict['reads.w.primary.only'], outdict['reads.w.secondary.only'], outdict['reads.w.both'], outdict['primary.err'],
                    outdict['missing'], outdict['primary.err'], outdict['pair.err'], outdict['pos.err'], outdict['chr.err'], 
                    outdict['exact.map'], outdict['close.map'], outdict['mismapped'], outdict['diff.chr'], outdict['unmapped']];
        outline = [str(col) for col in outline];
        PWS(",".join(outline), outfile, newline=False);
    ## Close unmapped file
## Close summary file

####################



# if random.uniform(0,1) < 0.01:
#     tracefile.write(getDateTime() + "\n");
#     tracefile.write("MATCH CLASS: " + classification + "\n");
#     tracefile.write("READ\n");
#     read_info = [read.query_name, read.is_read1, read.reference_name, read.reference_id];
#     read_info = [ str(r) for r in read_info ];
#     tracefile.write("\t".join(read_info) + "\n");
#     tracefile.write(str(read) + "\n");
#     tracefile.write("MATCH\n");
#     match_info = [match.query_name, match.is_read1, match.reference_name, match.reference_id];
#     match_info = [ str(r) for r in match_info ];
#     tracefile.write("\t".join(match_info) + "\n");
#     tracefile.write(str(match) + "\n");
#     tracefile.write("--------------------------\n");
# Random writing for manual checking.           