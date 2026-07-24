#############################################################################

import os
import sys
import yaml
from collections import defaultdict

#############################################################################

def trackerCount(trackerfile, data_type):
    tracker_counts = {};
    tracker = {};
    # Init vars to count tracks

    for line in open(trackerfile):
        line = line.strip().split("\t");
        # Parse the current line in the file

        region = line[0];
        overlap = line[-1];
        # Read the region and the number of overlaps
        ## CHECK

        if region not in tracker:
            tracker[region] = { };
        # Initialize the dict if the region hasn't been seen before     

        if overlap != "0":
            if data_type == "var":
                name = "-".join([ line[0], line[1], line[2] ]);
                iteration = int(line[3][-1]);
            elif data_type == "read":
                name = line[3] ## CHECK
                iteration = int(line[12][-1]);
            # The bed files for reads and snps are slightly different, so get name of feature and iteration here depending

            if name not in tracker[region]:
                tracker[region][name] = 0;
            # Initialize classification for this feature as 0

            tracker[region][name] += iteration;
            # Add the iteration number to the tracker for this feature
            # Tracker codes:
            # 0: feature not found in iter2 or 3
            # 2: feature found in iter2 only
            # 3: feature found in iter3 only
            # 5: feature found in both iters 2 and 3
    ## End tracker file line loop

    for region in tracker:
        tracker_counts[region] = { 'num.none' : 0, 'num.iter2.uniq' : 0, 'num.iter3.uniq' : 0, 'num.shared' : 0 };
        # Init the tracker counts for this region

        for name in tracker[region]:
            if tracker[region][name] == 0:
                tracker_counts[region]['num.none'] += 1;
            elif tracker[region][name] == 2:
                tracker_counts[region]['num.iter2.uniq'] += 1;
            elif tracker[region][name] == 3:
                tracker_counts[region]['num.iter3.uniq'] += 1;
            elif tracker[region][name] == 5:
                tracker_counts[region]['num.shared'] += 1;
            # Count the number of iterations this feature is found in according to the code above

        ## End name loop for tracker counts
    ## End region loop for tracker counts

    return tracker_counts;

#####################

def overlapCount(bedfilename):
    overlap_counts = {};
    overlaps = {};
    # Init vars to count and list overlaps

    for line in open(bedfilename):
        line = line.strip().split("\t");
        # Parse the current line in the file

        region = line[0];
        overlap = line[-1];
        # Read the region and the number of overlaps

        if region not in overlap_counts:
            overlap_counts[region] = { 'overlap' : 0, 'no-overlap' : 0 };
            overlaps[region] = [];
        # Initialize the dicts if the region hasn't been seen before

        if overlap == "0":
            overlap_counts[region]['no-overlap'] += 1;
        else:
            overlap_counts[region]['overlap'] += 1;
            overlaps[region].append("-".join([ line[0], line[1], line[2] ]));
        # Count whether an overlap has occurred

    return overlap_counts, overlaps;

#####################

def parseOverlaps(repeat_ols, gene_ols, total_features):
    repeat_ols = set(repeat_ols);
    gene_ols = set(gene_ols);
    # Convert overlap lists to sets for easy comparisons

    shared = len(repeat_ols.intersection(gene_ols));
    # Count the total number of FNs that overlap with a gene AND a repeat                        

    repeat_uniq = len(repeat_ols - gene_ols);
    # Count the total number of FNs that overlap ONLY with a repeat

    gene_uniq = len(gene_ols - repeat_ols);
    # Count the total number of FNs that overlap ONLY with a gene

    none = total_features - len(repeat_ols.union(gene_ols));
    # Count the number of FNs that overlap with NEITHER a gene or a repeat

    assert shared+repeat_uniq+gene_uniq+none == total_features;
    # Make sure that the numbers match up

    return repeat_uniq, gene_uniq, shared, none;

#############################################################################

config_file = "../simulation-configs/mm39-iterative.yaml";
outfilename = "../data/mm39-30X-0.005h-annotations.tsv";
trackeroutfile = "../data/mm39-30X-0.005h-tracker.tsv";
# Input and output files

var_types = ["tps", "fns"];
map_types = ["exact", "unmapped"];

with open(config_file, "r") as stream:
    try:
        config = yaml.safe_load(stream);
    except yaml.YAMLError as exc:
        print(exc);
# Read the sim config file

#sim_name = config['ref_str'] + "-" + "-".join(config['regions']) + "-iterative";

headers = ["type", "region", "coverage", "divergence", "heterozygosity", "iteration", "total", "num.repeat.uniq", "num.gene.uniq", "num.shared", "num.none" ];
tracker_headers = ["type", "region", "coverage", "divergence", "heterozygosity", "total", "num.lost", "num.iter2.uniq", "num.iter3.uniq", "num.shared" ];
# Headers for the output file

with open(outfilename, "w") as outfile, open(trackeroutfile, "w") as tracker:
    outfile.write("\t".join(headers) + "\n");
    tracker.write("\t".join(tracker_headers) + "\n");
    # Write the headers to the output file

    for cov in config['covs']:
        for div in config['divs']:
            if div == "0.00":
                continue;
            # Skip diviergence level 0
            for het in config['hets']:
                for n in range(1,config['iterations']+1):
                    for region in config['regions']:
                        for var_type in var_types:

                            if n == 1:
                                tracker_name = [ config['ref_str'], region, cov + "X", div + "d", het + "h", "compare-vcf-" + var_type + "-tracker.bed" ];
                                tracker_name = "-".join(tracker_name);
                                tracker_path = os.path.join(config['sim_outdir'], "summary-files", str(cov) + "X", str(div), "regions", tracker_name);
                                # Read the tracker bed file

                                tracker_counts = trackerCount(tracker_path, "var");
                                # Count the overlaps between iterations

                                total = tracker_counts[region]['num.none'] + tracker_counts[region]['num.iter2.uniq'] + tracker_counts[region]['num.iter3.uniq'] + tracker_counts[region]['num.shared'];
                                # Count the total number of SNPs

                                tracker_outline = [ var_type, region, cov, div, het, total, tracker_counts[region]['num.none'], tracker_counts[region]['num.iter2.uniq'], tracker_counts[region]['num.iter3.uniq'], tracker_counts[region]['num.shared'] ];
                                tracker_outline = [ str(o) for o in tracker_outline ];
                                tracker.write("\t".join(tracker_outline) + "\n");
                                # Output tracker info to file
                            ## Do tracker stuff on first iteration

                            ##########

                            repeat_name = [ config['ref_str'], "iter" + str(n), region, cov + "X", div + "d", het + "h", "compare-vcf-" + var_type + "-repeats.bed" ];
                            repeat_name = "-".join(repeat_name);
                            repeat_path = os.path.join(config['sim_outdir'], "summary-files", str(cov) + "X", str(div), "regions", repeat_name);
                            repeat_ol_counts, repeat_ols = overlapCount(repeat_path);
                            # Read the FN repeat overlaps

                            gene_name = [ config['ref_str'], "iter" + str(n), region, cov + "X", div + "d", het + "h", "compare-vcf-" + var_type + "-genes.bed" ];
                            gene_name = "-".join(gene_name);
                            gene_path = os.path.join(config['sim_outdir'], "summary-files", str(cov) + "X", str(div), "regions", gene_name);
                            gene_ol_counts, gene_ols = overlapCount(gene_path);
                            # Read the FN gene overlaps

                            total = repeat_ol_counts[region]['overlap'] + repeat_ol_counts[region]['no-overlap'];
                            assert total == gene_ol_counts[region]['overlap'] + gene_ol_counts[region]['no-overlap'];
                            # Count the total number of FNs and make sure they are consistent in both bed files

                            repeat_uniq, gene_uniq, shared, none = parseOverlaps(repeat_ols[region], gene_ols[region], total);
                            # Parse the overlaps

                            outline = [ var_type, region, cov, div, het, n, total, repeat_uniq, gene_uniq, shared, none ];
                            outline = [ str(o) for o in outline ];
                            outfile.write("\t".join(outline) + "\n");
                            # Compile and write the output line for this set of params
                        ## End var type loop
                    ## End region loop for variants

                    #####################

                    for map_type in map_types:

                        if n == 1:
                            tracker_name = [ config['ref_str'], cov + "X", div + "d", het + "h", "compare-bam-" + map_type + "-tracker.bed" ];
                            tracker_name = "-".join(tracker_name);
                            tracker_path = os.path.join(config['sim_outdir'], "summary-files", str(cov) + "X", str(div), tracker_name);
                            # Get tracker file name

                            tracker_counts = trackerCount(tracker_path, "read");
                            # Count the tracks

                            for region in config['regions']:
                                total = tracker_counts[region]['num.none'] + tracker_counts[region]['num.iter2.uniq'] + tracker_counts[region]['num.iter3.uniq'] + tracker_counts[region]['num.shared'];

                                tracker_outline = [ map_type, region, cov, div, het, total, tracker_counts[region]['num.none'], tracker_counts[region]['num.iter2.uniq'], tracker_counts[region]['num.iter3.uniq'], tracker_counts[region]['num.shared'] ];
                                tracker_outline = [ str(o) for o in tracker_outline ];
                                tracker.write("\t".join(tracker_outline) + "\n");
                            ## Output by region
                        ## Do tracker stuff on first iteration

                        ##########

                        repeat_name = [ config['ref_str'], cov + "X", div + "d", het + "h", "iter" + str(n), "compare-bam-" + map_type + "-repeats.bed" ]
                        repeat_name = "-".join(repeat_name);
                        repeat_path = os.path.join(config['sim_outdir'], "summary-files", str(cov) + "X", str(div), repeat_name);
                        repeat_ol_counts, repeat_ols = overlapCount(repeat_path);
                        # Read the bam repeat overlaps                             

                        gene_name = [ config['ref_str'], cov + "X", div + "d", het + "h", "iter" + str(n), "compare-bam-" + map_type + "-genes.bed" ]
                        gene_name = "-".join(gene_name);
                        gene_path = os.path.join(config['sim_outdir'], "summary-files", str(cov) + "X", str(div), gene_name);                        
                        gene_ol_counts, gene_ols = overlapCount(gene_path);
                        # Read the bam gene overlaps

                        for region in config['regions']:
                            total = repeat_ol_counts[region]['overlap'] + repeat_ol_counts[region]['no-overlap'];
                            assert total == gene_ol_counts[region]['overlap'] + gene_ol_counts[region]['no-overlap'];
                            # Count the total number of FNs and make sure they are consistent in both bed files

                            repeat_uniq, gene_uniq, shared, none = parseOverlaps(repeat_ols[region], gene_ols[region], total);
                            # Parse the overlaps

                            outline = [ map_type, region, cov, div, het, n, total, repeat_uniq, gene_uniq, shared, none ];
                            outline = [ str(o) for o in outline ];
                            outfile.write("\t".join(outline) + "\n");
                            # Compile and write the output line for this set of params

                            #sys.exit();
                        ## End region loop for reads
                    ## End map type loop
                ## End iteration loop
            ## End het loop
        ## End div loop
    ## End cov loop
## Close output file                

#############################################################################                    


