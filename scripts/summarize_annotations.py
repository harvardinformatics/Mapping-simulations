#############################################################################

import os
import sys
import yaml
from collections import defaultdict

#############################################################################

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
# Headers for the output file

with open(outfilename, "w") as outfile:
    outfile.write("\t".join(headers) + "\n");
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