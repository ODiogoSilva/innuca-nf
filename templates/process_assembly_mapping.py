#!/usr/bin/env python3

"""
process_assembly_mapping template for nextflow

Purpose
-------

This module is intended to process the coverage report from the assembly
mapping process

Expected input
--------------
fastq_id: Pair of FastQ file paths
    .: 'SampleA'
assembly: Fasta assembly file
    .: 'SH10761A.assembly.fasta'
coverage: TSV file with the average coverage for each assembled contig
    .: 'coverage.tsv'
min_assembly_coverage: Minimum coverage for assembled contigs. Can be 'auto'.
    .: 'auto' or '10'

Generated output
----------------

"""


FASTQ_ID = '$fastq_id'
ASSEMBLY_FILE = '$assembly'
COVERAGE_FILE = '$coverage'
MIN_ASSEMBLY_COVERAGE = '$min_assembly_coverage'


def main():

    print(FASTQ_ID)
    print(ASSEMBLY_FILE)
    print(COVERAGE_FILE)
    print(MIN_ASSEMBLY_COVERAGE)


main()
