#!/usr/bin/env python3

"""
trimmomatic template for nextflow

Purpose
-------

This module is intended execute trimmomatic on sets of paired end FastQ files

Expected input
--------------
fastq_id: Pair of FastQ file paths
    .: 'SampleA'
fastq_pair: Pair of FastQ file paths
    .: 'SampleA_1.fastq.gz SampleA_2.fastq.gz'
trim_range: Crop range detected using FastQC
    .: ''
opts: List of options for trimmomatic
    .: '["5:20", "3", "3", "55"]'
    .: '[trim_sliding_window, trim_leading, trim_trailing, trim_min_length]'
phred: List of guessed phred values for each sample
    .: '[SampleA: 33, SampleB: 33]'

Generated output
----------------

"""

#TODO: More control over read trimming
#TODO: Add option to remove adapters
#TODO: What to do when there is encoding failure

import os
import subprocess

from subprocess import PIPE

FASTQ_ID = '$fastq_id'
FASTQ_PAIR = '$fastq_pair'.split()
TRIM_RANGE = '$trim_range'.split()
TRIM_OPTS = [x.strip() for x in '$opts'.strip("[]").split(",")]
PHRED = '$phred'

TRIM_PATH = "/NGStools/Trimmomatic-0.36/trimmomatic.jar"


def main():

    # Create base CLI
    cli = [
        "java",
        "-jar",
        TRIM_PATH.strip(),
        "PE",
    ]

    # If the phred encoding was detected, provide it
    try:
        # Check if the provided PHRED can be converted to int
        phred = int(PHRED)
        phred_flag = "-phred{}".format(str(phred))
        cli += [phred_flag]
    # Could not detect phred encoding. Do not add explicit encoding to
    # trimmomatic and let it guess
    except ValueError:
        pass

    # Add input samples to CLI
    cli += FASTQ_PAIR

    # Add output file names
    output_names = []
    for fastq in FASTQ_PAIR:
        output_names.append("{}_P.fastq.gz".format(fastq.split(".")[0]))
        output_names.append("{}_U.fastq.gz".format(fastq.split(".")[0]))
    cli += output_names

    # Add trimmomatic options
    cli += [
        "CROP:{}".format(TRIM_RANGE[1]),
        "HEADCROP:{}".format(TRIM_RANGE[0]),
        "SLIDINGWINDOW:{}".format(TRIM_OPTS[0]),
        "LEADING:{}".format(TRIM_OPTS[1]),
        "TRAILING:{}".format(TRIM_OPTS[2]),
        "MINLEN:{}".format(TRIM_OPTS[3]),
        "TOPHRED33"
    ]

    p = subprocess.Popen(cli, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    # Check if trimmomatic ran successfully. If not, write the error message
    # to the status channel and exit.
    with open("trimmomatic_status", "w") as fh:
        if p.returncode != 0:
            fh.write(str(stderr))
            return
        else:
            fh.write("pass")

main()
