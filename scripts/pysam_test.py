# samtools view -h /n/holylfs05/LABS/informatics/Users/gthomas/Mapping-simulations-data/mm39-18-19/simulated-reads/30X/0.08/heterozygous/0.005/mm39_golden.bam | head -n107 | samtools view -b > test.bam
# samtools index test.bam

import pysam

bam_file = "test.bam";
outfile = "pysam.bam";

bam = pysam.AlignmentFile(bam_file, "rb");

with pysam.AlignmentFile(outfile, "wb", header=bam.header) as bam_out:
    for read in bam.fetch():
        print(read);
        print(read.reference_name);
        bam_out.write(read);