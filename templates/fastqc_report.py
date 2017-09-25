#!/usr/bin/env python3

"""
fastqc_report template for nextflow

Purpose
-------

This module is intended parse the results of FastQC for paired end FastQ
samples

Expected input
--------------
fastq_id: Sample Identification string
    .: 'SampleA'
result_p1: Path to FastQC result files for pair 1
    .: 'SampleA_1_fastqc/summary.txt', 'SampleA_1_fastqc/fastqc_data.txt'
result_p2 Path to FastQC result files for pair 2
    .: 'SampleA_2_fastqc/summary.txt', 'SampleA_2_fastqc/fastqc_data.txt'

Generated output
----------------


"""

RESULT_P1 = '$result_p1'.split()
RESULT_P2 = '$result_p2'.split()
FASTQ_ID = '$fastq_id'


def get_trim_index(biased_list):
    """

    Parameters
    ----------
    biased_list

    Returns
    -------

    """

    # Return index 0 if there are no biased positions
    if set(biased_list) == {False}:
        return 0

    # Iterate over the biased_list array. Keep the iteration going until
    # we find a biased position with the two following positions unbiased
    # (e.g.: True, False, False).
    # When this condition is verified, return the last biased position
    # index for subsequent trimming.
    for i, val in enumerate(biased_list):
        if val and set(biased_list[i+1:i+3]) == {False}:
            return i + 1

    # If the previous iteration could not find and index to trim, it means
    # that the whole list is basically biased. Return the length of the
    # biased_list
    return len(biased_list)


def trim_range(data_file):
    """

    Parameters
    ----------
    data_file

    Returns
    -------

    """

    # Target string for nucleotide bias assessment
    target_nuc_bias = ">>Per base sequence content	fail"
    # This flag will become True when gathering base proportion data
    # from file.
    gather = False

    # This variable will store a boolean array on the biased/unbiased
    # positions. Biased position will be True, while unbiased positions
    # will be False
    biased = []

    with open(data_file) as fh:

        for line in fh:
            # Start assessment of nucleotide bias
            if line.startswith(target_nuc_bias):
                # Skip comment line
                next(fh)
                gather = True
            # Stop assessment when reaching end of target module
            elif line.startswith(">>END_MODULE") and gather:
                break
            elif gather:
                # Get proportions of each nucleotide
                g, a, t, c = [float(x) for x in line.strip().split()[1:]]
                # Get 'GC' and 'AT content
                gc = (g + 0.1) / (c + 0.1)
                at = (a + 0.1) / (t + 0.1)
                # Assess bias
                if 0.8 <= gc <= 1.2 and 0.8 <= at <= 1.2:
                    biased.append(False)
                else:
                    biased.append(True)

    # Split biased list in half to get the 5' and 3' ends
    biased_5end, biased_3end = biased[:int(len(biased))],\
                               biased[int(len(biased)):][::-1]

    trim_nt = [0, 0]
    # Assess number of nucleotides to clip at 5' end
    trim_nt[0] = get_trim_index(biased_5end)
    # Assess number of nucleotides to clip at 3' end
    trim_nt[1] = len(biased) - get_trim_index(biased_3end)

    return trim_nt


def get_sample_trim(p1_data, p2_data):
    """

    Parameters
    ----------
    p1_data
    p2_data

    Returns
    -------

    """

    sample_ranges = [trim_range(x) for x in [p1_data, p2_data]]

    # Get the optimal trim position for 5' end
    optimal_5trim = max([x[0] for x in sample_ranges])
    # Get optimal trim position for 3' end
    optimal_3trim = min(x[1] for x in sample_ranges)

    return optimal_5trim, optimal_3trim


def main():

    # Get optimal trimming range for sample, based on the per base sequence
    # content
    optimal_trim = get_sample_trim(RESULT_P1[0], RESULT_P2[0])

    with open("optimal_trim", "w") as fh:
        fh.write("{}".format(" ".join([str(x) for x in optimal_trim])))


main()
