#!/usr/bin/env python3

"""
fastqc template for nextflow

Purpose
-------

This module is intended to run FastQC on paired-end FastQ files.

Expected input
--------------
fastq_id: Sample Identification string
    .: 'SampleA'
fastq_pair: Pair of FastQ file paths
    .: 'SampleA_1.fastq.gz SampleA_2.fastq.gz'

Generated output
----------------
pair_[1,2]_data: File containing FastQC report at the nucleotide level
    .: 'pair_1_data' and 'pair_2_data'
pair_[1,2]_summary: File containing FastQC report for each category
    .: 'pair_1_summary' and 'pair_2_summary'
"""

import os
import subprocess

from subprocess import PIPE
from os.path import exists, join

FASTQ_PAIR = '$fastq_pair'.split()
FASTQ_ID = '$fastq_id'
ADAPTER_FILE = eval('$ad')


def convert_adatpers(adapter_fasta):
    """Generates an adapter file for FastQC from a fasta file

    Parameters
    ----------
    adapter_fasta : str
        Path to Fasta file with adapter sequences

    Returns
    -------

    """

    adapter_out = "fastqc_adapters.tab"

    try:

        with open(adapter_fasta) as fh, \
                open(adapter_out) as adap_fh:

            for line in fh:
                if line.startswith(">"):

                    head = line[1:].strip()
                    sequence = next(fh).strip()

                    adap_fh.write("{}\\t{}\\n".format(head, sequence))

        return adapter_out

    except FileNotFoundError:
        return


def main():

    # If an adapter file was provided, convert it to FastQC format
    if ADAPTER_FILE:
        adapter_file = convert_adatpers(ADAPTER_FILE)
    else:
        adapter_file = None

    # Setting command line for FastQC
    cli = [
        "fastqc",
        "--extract",
        "--nogroup",
        "--format",
        "fastq",
        "--threads",
        "$task.cpus"
    ]

    # Add adapters file to command line, if it exists
    if adapter_file:
        cli += ["--adapters {}".format(adapter_file)]

    # Add FastQ files at the end of command line
    cli += FASTQ_PAIR

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # Check if the FastQC output was correctly generated.
    with open("fastq_status", "w") as fh:
        for fastq in FASTQ_PAIR:
            fpath = join(fastq.split(".")[0] + "_fastqc", "fastqc_data.txt")
            # If the FastQC output does not exist, pass the STDERR to
            # the output status channel and exit
            if not exists(fpath):
                fh.write(str(stderr))
                return

        # If the output directories exist, write 'pass' to the output status
        # channel
        fh.write("pass")

    # Both FastQC have been correctly executed. Get the relevant FastQC
    # output files for the output channel
    for i, fastq in enumerate(FASTQ_PAIR):
        # Get results for each pair
        fastqc_dir = fastq.split(".")[0] + "_fastqc"

        summary_file = join(fastqc_dir, "summary.txt")
        fastqc_data_file = join(fastqc_dir, "fastqc_data.txt")

        os.rename(fastqc_data_file, "pair_{}_data".format(i + 1))
        os.rename(summary_file, "pair_{}_summary".format(i + 1))


main()
