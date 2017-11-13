from process_templates.ProcessSkeleton import Skeletons


class Process:

    def __init__(self, ptype=None, template=None, pid=None):

        accepted_types = [
            "pre_assembly", "assembly", "post_assembly", "annotation"]

        if ptype not in accepted_types:
            raise ValueError(
                "{} is not an accepted process type".format(ptype))

        self.skels = Skeletons()

        self.pid = pid
        self.ptype = ptype
        self.template = template

        if ptype == "pre_assembly":
            self.main_in_str = self.main_out_str = "MAIN_fq"
        elif ptype == "assembly":
            self.main_in_str = "MAIN_fq"
            self.main_out_str = "MAIN_assembly"
        else:
            self.main_in_str = self.main_out_str = "MAIN_assembly"

        self.link_start = [self.main_out_str]
        self.link_end = []

        self._input_channel = None
        self._output_channel = None
        self.output_type = None

        self.kwargs = None
        self.forks = []

        print(self.template)
        print(self.pid)

    @property
    def template_str(self):

        if self.kwargs:
            x = self.skels[self.template].format(**self.kwargs)

            if self.forks:
                x += "\n".join(self.forks)

            return x

        else:
            raise Exception("Channels must be setup first using the "
                            "set_channels method")

    def set_channels(self, **kwargs):
        """ General purpose method that sets the template placeholders

        In most cases, the kwargs dict will contain at least the ``pid``
        (process_id) key.

        Parameters
        ----------
        kwargs : dict
            Dictionary with the placeholders for the template string.
        """

        self.pid = kwargs.get("pid")

        self._input_channel = "{}_{}".format(self.main_in_str, self.pid)

        if self.output_type:
            self._output_channel = "{}_{}".format(self.main_out_str,
                                                  self.pid + 1)

        self.kwargs = {**kwargs, **{"input_channel": self._input_channel,
                                    "output_channel": self._output_channel}}

    def set_secondary_channels(self, source, channel_list):
        """ Forks the start of a secondary channel to multiple channels

        The fork will be based on the pid_list, which should be a list of
        numbers.

        Parameters
        ----------
        source : str
            String with the name of the source channel
        channel_list : list
            List of channels that will receive a fork of the secondary
            channel
        """

        if source == "MAIN_assembly":
            self.kwargs["output_channel"] = ",".join(channel_list)
            return

        # Handle the special case where the main channel is forked
        elif source == "MAIN_fq":
            # Update previous output_channel
            self.kwargs["output_channel"] = "_{}".format(self._output_channel)
            # Set source to modified output channel
            source = self.kwargs["output_channel"]
            # Modify original channel_list for non duplicate main channels
            channel_list = ["_{}".format(x) for x in channel_list]
            # Add original output channel to channel_list
            channel_list.append(self._output_channel)
        else:
            source = "{}_{}".format(source, self.pid)

        if len(channel_list) == 1:
            self.forks.append("\n{}.into{{ {} }}\n".format(source,
                                                           channel_list[0]))

        else:
            self.forks.append("\n{}.into{{ {} }}\n".format(
                source, ";".join(channel_list)))


class IntegrityCoverage(Process):
    """Process template interface for first integrity_coverage process

    This process is mandatory as the first step in a pipeline

    The required template fields are:

        - `pid`

    """

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.link_start.extend(["SIDE_phred", "SIDE_max_len"])


class CheckCoverage(Process):
    """Process template interface for additional integrity_coverage process

    This process can be added before the assembly process multiple times

    The required template fields are:

        - `pid`
        - `input_channel`
        - `output_channel`

    """

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.link_start.extend(["SIDE_max_len"])


class FastQC(Process):
    """FastQC process template interface

    The required template fields are:

        -

    """

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"


class Trimmomatic(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.link_end.append({"link": "SIDE_phred", "alias": "SIDE_phred"})


class Spades(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "assembly"

        self.link_end.append({"link": "SIDE_max_len", "alias": "SIDE_max_len"})


class ProcessSpades(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="post_assembly",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = "assembly"


class AssemblyMapping(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="post_assembly",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = "assembly"

        self.link_end.append({"link": "MAIN_fq", "alias": "MAIN_assembly"})


class Pilon(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="post_assembly",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = "assembly"


class Mlst(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="annotation",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = None

        self.link_start = None
        self.link_end.append({"link": "MAIN_assembly",
                              "alias": "MAIN_assembly"})


class Abricate(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="annotation",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = None

        self.link_start = None
        self.link_end.append({"link": "MAIN_assembly",
                              "alias": "MAIN_assembly"})


class Prokka(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="annotation",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = None

        self.link_start = None
        self.link_end.append({"link": "MAIN_assembly",
                              "alias": "MAIN_assembly"})

