#!/usr/bin/env python3

"""
trimmomatic_report template for nextflow

Purpose
-------

This module is intended parse the results of the Trimmomatic log for a set
of samples

Expected input
--------------
log_files: Trimmomatic log files
    .: 'Sample1_trimlog.txt Sample2_trimlog.txt'

Generated output
----------------

"""


from collections import OrderedDict


if __file__.endswith(".command.sh"):
    LOG_FILES = '$log_files'.split()


def main(log_files):

    log_storage = OrderedDict()

    template = {
        # Total trimmed base pairs
        "total_trim": 0,
        # Total length after trimming
        "total_len": 0,
        # Total trimmed at 5' end
        "5trim": 0,
        # Total trimmed at 3' end
        "3trim": 0
    }

    for log in log_files:

        log_id = log.split("_")[0]

        # Copy the template dictionary to the current sample
        log_storage[log_id] = dict(template.items())

        with open(log) as fh:

            for line in fh:

                # This will split the log fields into:
                # 0. read length after trimming
                # 1. amount trimmed from the start
                # 2. last surviving base
                # 3. amount trimmed from the end
                fields = [int(x) for x in line.strip().split()[2:]]

                log_storage[log_id]["5trim"] += fields[1]
                log_storage[log_id]["3trim"] += fields[3]
                log_storage[log_id]["total_trim"] += fields[1] + fields[3]
                log_storage[log_id]["total_len"] += fields[0]


if __name__ == '__main__':
    main(LOG_FILES)