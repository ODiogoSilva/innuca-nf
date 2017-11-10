#!/usr/bin/env python

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

    def __init__(self, process_list, nextflow_file):

        # Check if all specified processes are available
        for p in process_list:
            if p not in self.process_map:
                raise ValueError(
                    "The process '{}' is not available".format(p))

        self.processes = [self.process_map[p](p) for p in process_list]
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

        print(self.secondary_channels)

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


def main():

    pipeline = [
        "integrity_coverage",
        # "check_coverage",
        # "fastqc_trimmomatic",
        "fastqc",
        "trimmomatic",
        # "trimmomatic",
        # "fastqc",
        # "check_coverage",
        # "trimmomatic",
        # "fastqc_trimmomatic",
        # "fastqc",
        "spades",
        # "process_spades",
        "assembly_mapping",
        "pilon",
        "mlst",
        "abricate",
        "prokka"
    ]

    nfg = NextflowGenerator(pipeline, "custom_pipe.nf")

    nfg.build()


if __name__ == '__main__':
    main()
