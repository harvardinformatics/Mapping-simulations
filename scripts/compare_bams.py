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

regions, coverage, divergence, heterozygosity, iteration, ref_index, golden_bam_file, query_bam_file, summary_outfilename = sys.argv[1:];
regions = regions.split(",");
# Inputs

# divergence = "0.04";
# heterozygosity = "0.005";
# iteration = "1";
# coverage = "30";
# regions = ["18"];
# golden_bam_file = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19/simulated-reads/30X/0.04/heterozygous/0.005/mm39_golden.bam";
# query_bam_file = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-iterative/mapped-reads/30X/0.04/0.005/iter1/mm39-30X-0.04d-0.005h-crossmap.sorted.bam";
# ref_index = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/reference-genomes/mm39/Mus_musculus.GRCm39.dna.primary_assembly.chromes.fa.fai";
# summary_outfilename = "test.csv";
#tracefilename = "/n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19/summary-files/20X/0.02/mm39-20X-0.02-compare-bam-trace.txt";
# I/O options

pad = 20;
with open(summary_outfilename, "w") as outfile:#, open(tracefilename, "w") as tracefile:
    runTime("# Compare BAMs", outfile);
    PWS(spacedOut("# Golden BAM:", pad) + golden_bam_file, outfile);
    PWS(spacedOut("# Query BAM:", pad) + query_bam_file, outfile);
    PWS(spacedOut("# Output file:", pad) + summary_outfilename, outfile);
    PWS("# ----------------", outfile);

    ####################

    #window_size = 10000
    PWS("# " + getDateTime() + " Reading reference index...", outfile);
    chr_lens = { line.split("\t")[0] : int(line.split("\t")[1]) for line in open(ref_index, "r") };
    #print(chr_lens);
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
    # Takes about half an hour and uses almost 50GB of RAm on 40GB test BAM

    ####################

    outdict = {"exact-map" : 0, "close-map" : 0, "mismapped" : 0, "diff-chr" : 0, "unmapped" : 0, "missing" : 0 };
    # Just doing counts for now. Could do full distributions of distance between golden and mapped positions, but
    # would have to deal with the memory issues

    #query_reads = query_bam.fetch();

    PWS("# " + getDateTime() + " Finding reads...", outfile);
    for chrome in chr_lens:
    # Go through the reads by chromosome

        if chrome not in regions:
            continue;
        # For debugging/testing

        PWS("# " + getDateTime() + " Finding reads in " + chrome + "...", outfile);

        num_reads = 0;
        num_reads_missing = 0;

        for read in golden_bam.fetch(chrome):
        # Find every golden read in the current chromosome

            if num_reads % 1000 == 0:
                print(getDateTime() + " " + str(num_reads));
            num_reads += 1;
            # print(read);
            # print(read.query_name);
            # print(read.is_read1);
            # print(read.is_secondary);
            # print("----");

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

            # print(golden_chr);
            # print(golden_pos);

            try:
                matches = name_indexed_query.find(read.query_name);
            except:
                outdict['missing'] += 1;
                continue;
            # Find all mappings for the current read in the query bam file

            matches_found = 0;
            # Keeps track of how many matches we categorize - should only be one primary match

            for match in matches:
            # Go through every matching read in the query file for the current read

                if match.is_secondary:
                    continue;
                # Only concerned with primary mappings

                if (match.is_read1 and read_pair == "1") or (match.is_read2 and read_pair == "2"):
                # Make sure they are the same read in the pair

                    classification = "NONE";
                    # if read.query_name == "mm39-19-19-1369":
                    #     print(read)
                    #     print(read.reference_name, read.reference_id, read.is_read1, read.query_name);
                    #     print(match);
                    #     print(match.reference_name, match.reference_id, match.is_read1, match.query_name);
                    #     print("-------")                      
                    #print(match);
                    # mm39-19-19-1613719
                    query_chr = match.reference_name;
                    #print(query_chr);
                    query_pos_list = match.get_reference_positions();
                    #print(query_pos_list);
                    if query_pos_list:                     
                        query_pos = query_pos_list[0];
                        #print(query_pos);
                        #print("------");
                        # Get the chromosome and position of where the read actually mapped

                        if golden_chr != query_chr:
                            classification = "DIFF CHR";
                            outdict['diff-chr'] += 1;
                            matches_found += 1;

                            # print(read)
                            # print(read.reference_name, read.reference_id, read.is_read1, read.query_name);
                            # print(match);
                            # print(match.reference_name, match.reference_id, match.is_read1, match.query_name);
                            # sys.exit();
                        # If it mapped to a different chromosome completely

                        elif golden_chr == query_chr:
                            if golden_pos == query_pos:
                                classification = "EXACT MAP"
                                outdict['exact-map'] += 1;
                                matches_found += 1;
                            # An exact match

                            elif abs(query_pos - golden_pos) <= 150:
                                classification = "CLOSE MAP"
                                outdict['close-map'] += 1;
                                matches_found += 1;
                            # A mapping within one read length of the true mapping

                            else:
                                classification = "MISMAPPED"
                                outdict['mismapped'] += 1;
                                matches_found += 1;
                        # A mapping on the same chromosome, but a different position

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
            ## End matches loop

            if matches_found == 0:
                outdict['unmapped'] += 1;
            # If there are no primary matches found, increment as unmapped

            # elif matches_found > 1:
            #     print(read)
            #     print(matches);
            # For testing

            # if num_reads > 2000000:
            #     break;
        ## End golden read loop
    ## End chromosome loop
        
    # print(num_reads);
    # print(outdict);

    #outdict = {"exact-map" : 0, "close-map" : 0, "mismapped" : 0, "diff-chr" : 0, "unmapped" : 0 };

    headers = ["coverage", "divergence", "heterozygosity", "iteration", "missing.in.query", "exact.map", "close.map", "mismapped", "diff.chr", "unmapped"];
    PWS(",".join(headers), outfile);

    outline = [coverage, divergence, heterozygosity, iteration, outdict['missing'], outdict['exact-map'], outdict['close-map'], outdict['mismapped'], outdict['diff-chr'], outdict['unmapped']];
    outline = [str(col) for col in outline];
    PWS(",".join(outline), outfile, newline=False);

####################



# for line in gzip.open(golden_bam_file):
#     line = line.decode().strip().split("\t");
#     print(line);

# sys.exit();











# for chrome in chr_lens:
#     if chrome != "19":
#         continue;

#     w_start = 0;
#     w_end = window_size;

#     while True:
#         x = bamfile.count_coverage(contig=chrome, start=w_start, stop=w_end);
#         # print(x);
#         # print(len(x));
#         # print(x[0]);
#         # print(len(x[0]))
#         # print(x[0][0]);
#         # print(sum(x[0]));


#         print(w_start, w_end, sum( [ sum(x[0]), sum(x[1]), sum(x[2]), sum(x[3]) ] ));


#         w_start = w_end;
#         w_end = w_end + window_size;

#         if w_end > chr_lens[chrome]:
#             w_end = chr_lens[chrome];
#             x = bamfile.count_coverage(contig=chrome, start=w_start, stop=w_end);
#             #print(x);
#             print(w_start, w_end, sum( [ sum(x[0]), sum(x[1]), sum(x[2]), sum(x[3]) ] ));
#             break;

# sys.exit();



# # read_iter = bamfile.pileup("19", 3088636, 3088637);
# # #print(len(list(read_iter)));
# # for x in read_iter:
# #     print(x);
# #     print('---')


# # for read in bamfile.fetch("19", 3088636, 3088836):
# #     print(read);


# i = 1;
# for pileupcolumn in bamfile.pileup(contig="19", start=3088636, stop=3088637):
#     print (i, "coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n));
#     i += 1;
#     # for pileupread in pileupcolumn.pileups:
#     #     if not pileupread.is_del and not pileupread.is_refskip:
#     #         # query position is None if is_del or is_refskip is set.
#     #         print ('\tbase in read %s = %s' %
#     #               (pileupread.alignment.query_name,
#     #                pileupread.alignment.query_sequence[pileupread.query_position]))


# # x = bamfile.get_index_statistics()
# # for y in x:
# #     print(y);

# x = bamfile.count_coverage(contig="19", start=3088636, stop=3088637);
# print(x);
# # Total coverage for all positions in given interval