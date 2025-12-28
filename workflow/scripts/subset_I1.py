import gzip, pysam

def filter_common_reads(path_i1, path_r1, path_I1_out):
    with pysam.FastqFile(path_i1) as vcf_i1, \
         pysam.FastqFile(path_r1) as vcf_r1, \
         gzip.open(path_I1_out, 'wt') as out:

        ## Initialization
        iter_i1 = iter(vcf_i1)
        iter_r1 = iter(vcf_r1)

        # Get Reads
        current_i1 = next(iter_i1, None)
        current_r1 = next(iter_r1, None)
        while current_i1 is not None and current_r1 is not None:
            if current_i1.name == current_r1.name:
                # Match found - write to output
                out.write(f"@{current_i1.name}\n")
                out.write(f"{current_i1.sequence}\n")
                out.write("+\n")
                out.write(f"{current_i1.quality}\n")

                # Advance both iterators
                current_i1 = next(iter_i1, None)
                current_r1 = next(iter_r1, None)
            else:
                # Read missing in r1, skip in i1
                current_i1 = next(iter_i1, None)


if __name__ == "__main__":
    filter_common_reads(path_i1=snakemake.input.fastq_i1,
                        path_r1=snakemake.input.fastq_r1,
                        path_I1_out=snakemake.output.fastq_i1)
