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
    .: 'SampleA_1_fastqc/data', 'SampleA_1_fastqc/summary'
result_p2 Path to FastQC result files for pair 2
    .: 'SampleA_2_fastqc/data', 'SampleA_2_fastqc/summary'

Generated output
----------------
fastqc_health: Stores the health check for the current sample. If it passes
    all checks, it contains only the string 'pass'. Otherwise, contains
    the summary categories and their respective results
    .: 'pass'
optimal_trim: Stores a tuple with the optimal trimming positions for 5' and
    3' ends of the reads.
    .: '15 151'

"""

from collections import OrderedDict


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


def get_summary(summary_file):
    """

    Parameters
    ----------
    summary_file

    Returns
    -------

    """

    summary_info = OrderedDict()

    with open(summary_file) as fh:
        for line in fh:
            # Skip empty lines
            if not line.strip():
                continue
            # Populate summary info
            fields = [x.strip() for x in line.split("\t")]
            summary_info[fields[1]] = fields[0]

    return summary_info


def check_summary_health(summary_file):
    """Checks the health of a sample from the FastQC summary file.

    Parses the FastQC summary file and tests whether the sample is good
    or not. There are four categories that cannot fail, and two that
    must pass in order to the sample pass this check

    Parameters
    ----------
    summary_file: str
        Path to FastQC summary file

    Returns
    -------
    _ : Boolean
        Returns True if the sample passes all tests. False if not.
    summary_info : dict
        A dictionary with the FastQC results for each category.

    """

    # Store the summary categories that cannot fail. If they fail, do not
    # proceed with this sample
    fail_sensitive = [
        "Per base sequence quality",
        "Overrepresented sequences",
        "Sequence Length Distribution",
        "Per sequence GC content"
    ]

    # Store summary categories that must pass. If they do not, do not proceed
    # with that sample
    must_pass = [
        "Per base N content",
        "Adapter Content"
    ]

    # Get summary dictionary
    summary_info = get_summary(summary_file)

    for cat, test in summary_info.items():

        # Check for fail sensitive
        if cat in fail_sensitive and test == "FAIL":
            return False, summary_info

        # Check for must pass
        if cat in must_pass and test != "PASS":
            return False, summary_info

    # Passed all tests
    return True, summary_info


def main():

    with open("fastqc_health", "w") as health_fh, \
            open("optimal_trim", "w") as trim_fh:

        # Perform health check according to the FastQC summary report for
        # each pair. If both pairs pass the check, send the 'pass' information
        # to the 'fastqc_health' channel. If at least one fails, send the
        # summary report.
        for fastqc_summary in [RESULT_P1[1], RESULT_P2[1]]:

            health, summary_info = check_summary_health(fastqc_summary)

            # If one of the health flags returns False, send the summary report
            # through the status channel
            if not health:
                for k, v in summary_info.items():
                    health_fh.write("{}: {}\\n".format(k, v))
                    trim_fh.write("fail")
                return
            else:
                health_fh.write("pass")

        # Get optimal trimming range for sample, based on the per base sequence
        # content
        optimal_trim = get_sample_trim(RESULT_P1[0], RESULT_P2[0])

        trim_fh.write("{}".format(" ".join([str(x) for x in optimal_trim])))


main()
