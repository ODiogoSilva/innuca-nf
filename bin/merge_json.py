#!/usr/bin/env python3

import sys
import json

core_file, f1, f2 = sys.argv[1:]


def get_core_genes(core_file):

    with open(core_file) as fh:
        core_genes = [x.strip() for x in fh.readlines()[1:]
                      if x.strip() != ""]

    return core_genes


def filter_core_genes(locus_array, info_array, core_genes):

    core_array = []

    for gene, info in zip(*[info_array, locus_array]):
        if gene in core_genes:
            core_array.append(info)

    return core_array


def assess_quality(core_array, core_genes):

    locus_not_found = core_array.count("LNF")
    perc = float(locus_not_found) / float(len(core_genes))

    # Fail sample with higher than 2% missing loci
    with open(".status", "w") as fh:
        if perc > 0.02:
            status = "fail"
        elif perc > 0.003:
            status = "warning"
        else:
            status = "pass"

        fh.write(status)

    return status, perc


def main():
    core_genes = get_core_genes(core_file)

    with open(f1) as f1h, open(f2) as f2h:

        j1 = json.load(f1h)
        j2 = json.load(f2h)

        current_result = j1["sample_polished.assembly.fasta"]
        current_array = j1["header"]
        core_results = filter_core_genes(current_result, current_array,
                                         core_genes)
        status, perc = assess_quality(core_results, core_genes)

        res = {"cagao": [j1, j2], "status": status, 'lnfPercentage': perc}

        with open(".report.json", "w") as fh:
            fh.write(json.dumps(res))


main()
