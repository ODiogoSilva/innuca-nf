#!/usr/bin/env python

import argparse


from process_templates import HeaderSkeleton as hs
from process_templates.Process import IntegrityCoverage, FastQC, Trimmomatic, \
    Spades, ProcessSpades, AssemblyMapping, Pilon, CheckCoverage, Mlst, \
    Abricate, Prokka


class ProcessError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class ChannelError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class NextflowGenerator:

    process_map = {
        "integrity_coverage": IntegrityCoverage,
        "check_coverage": CheckCoverage,
        "fastqc": FastQC,
        "trimmomatic": Trimmomatic,
        "fastqc_trimmomatic": Trimmomatic,
        "spades": Spades,
        "process_spades": ProcessSpades,
        "assembly_mapping": AssemblyMapping,
        "pilon": Pilon,
        "mlst": Mlst,
        "abricate": Abricate,
        "prokka": Prokka
    }
    """
    dict: Maps the process ids to the corresponding template interface class
    """

    def __init__(self, process_list, nextflow_file, process_ids=None):

        # Check if all specified processes are available
        for p in process_list:
            if p not in self.process_map:
                raise ValueError(
                    "The process '{}' is not available".format(p))

        if process_ids:
            if len(process_ids) != len(process_list):
                raise ProcessError(
                    "The provided list of process ids must match the length"
                    " of the process list."
                )

        else:
            process_ids = [None] * len(process_list)

        self.processes = [
            self.process_map[p](template=p, pid=pid) for p, pid in
            zip(process_list, process_ids)
        ]
        """
        list: Stores the process interfaces in the specified order
        """

        self.nf_file = nextflow_file
        """
        str: Path to file where the pipeline will be generated
        """

        self.template = ""
        """
        str: String that will harbour the pipeline code
        """

        self.secondary_channels = {}
        """
        dict: Stores secondary channel links
        """

        self._check_pipeline_requirements()

    def _check_pipeline_requirements(self):

        pipeline_types = [x.ptype for x in self.processes]
        pipeline_names = [x.template for x in self.processes]

        # Check if the pipeline contains at least one integrity_coverage
        # process
        if pipeline_names.index("integrity_coverage") != 0:
            raise ProcessError("The pipeline must contain at least one"
                               "integrity coverage at the start of the"
                               "pipeline")

        # Check if the pipeline contains an assembly process
        if pipeline_types.count("assembly") != 1:
            raise ProcessError("The pipeline must contain one and only one"
                               " assembly process")

    def _build_header(self):

        self.template += hs.header + hs.start_channel

    def _set_channels(self):

        previous_channel = None

        for idx, p in enumerate(self.processes):

            if not previous_channel:
                # Set the first output type
                previous_channel = p
            else:
                # Check if the connecting processes can be linked by their
                # input/output types
                if p.ptype == "annotation":
                    pass
                elif previous_channel.output_type != p.input_type:
                    raise ChannelError(
                        "The output of the '{}' process ({}) cannot link with"
                        " the input of the '{}' process ({}). Please check the"
                        " order of the processes.".format(
                            previous_channel.template,
                            previous_channel.output_type,
                            p.template,
                            p.input_type
                        ))

                previous_channel = p

            # Make sure that the process id starts at 1
            pidx = idx + 1

            # Check if the current process has a start of a secondary
            # side channel
            if p.link_start:
                for l in p.link_start:
                    self.secondary_channels[l] = {"p": p, "end": []}

            # check if the current process receives a secondary side channel.
            # If so, add to the links list of that side channel
            if p.link_end:
                for l in p.link_end:
                    if l["link"] in self.secondary_channels:
                        self.secondary_channels[l["link"]]["end"].append(
                            "{}_{}".format(l["alias"], pidx))

            p.set_channels(**{"pid": pidx})

    def _set_secondary_channels(self):

        for source, vals in self.secondary_channels.items():

            # Skip if there are no receiving ends for this secondary channel
            if not vals["end"]:
                continue

            vals["p"].set_secondary_channels(source, vals["end"])

    def build(self):

        # Generate regular nextflow header that sets up the shebang, imports
        # and all possible initial channels
        self._build_header()

        self._set_channels()

        self._set_secondary_channels()

        for p in self.processes:
            self.template += p.template_str

        with open(self.nf_file, "w") as fh:
            fh.write(self.template)


def get_args():

    parser = argparse.ArgumentParser(
        description="Nextflow pipeline generator")

    parser.add_argument("-t", "--tasks", nargs="+", dest="tasks",
                        type=get_tuples,
                        help="Space separated tasks of the pipeline")
    parser.add_argument("-o", dest="output_nf",
                        help="Name of the pipeline file")

    args = parser.parse_args()

    return args


def get_tuples(task):

    task_id = task.split(":")

    if len(task_id) != 2:
        raise argparse.ArgumentTypeError(
            "Tasks arguments must be in a format of <task>:<task id> "
            "(e.g.: 'spades:5')")

    # print(tasks_ids)
    return task_id


def main(args):

    # pipeline = [
    #     "integrity_coverage",
    #     # "check_coverage",
    #     # "fastqc_trimmomatic",
    #     "fastqc",
    #     "trimmomatic",
    #     # "trimmomatic",
    #     # "fastqc",
    #     # "check_coverage",
    #     # "trimmomatic",
    #     # "fastqc_trimmomatic",
    #     # "fastqc",
    #     "spades",
    #     # "process_spades",
    #     "assembly_mapping",
    #     "pilon",
    #     "mlst",
    #     "abricate",
    #     "prokka"
    # ]

    # Get process names
    process_names = [x[0] for x in args.tasks]

    # Get process ids
    process_ids = [x[1] for x in args.tasks]

    # nfg = NextflowGenerator(args.tasks, args.output_nf)
    nfg = NextflowGenerator(process_list=process_names,
                            process_ids=process_ids,
                            nextflow_file=args.output_nf)

    nfg.build()


if __name__ == '__main__':

    args = get_args()

    # print(args)

    main(args)
