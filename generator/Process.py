import os
import jinja2
import logging

from os.path import dirname, join, abspath

logger = logging.getLogger("main.{}".format(__name__))


class Process:

    def __init__(self, ptype, template, process_id=None):

        accepted_types = [
            "pre_assembly",
            "assembly",
            "post_assembly",
            "annotation",
            "status"
        ]

        if ptype not in accepted_types:
            raise ValueError(
                "{} is not an accepted process type".format(ptype))

        self.pid = None
        """
        int: Process ID
        """

        self.process_id = process_id

        self.ptype = ptype
        """
        str: Process type. See :py:attr:`accepted_types`.
        """

        self.template = template
        """
        str: Template name for the current process. This string will be used
        to fetch the correct template in :py:attr:`_skels`.
        """
        self._template_path = None
        self._set_template(template)

        self.input_type = None
        """
        str: Type of expected input data. Used to verify the connection between
        two processes is viable.
        """

        self.output_type = None
        """
        str: Type of output data. Used to verify the connection between
        two processes is viable.
        """

        self.dependencies = []
        """
        list: Contains the dependencies of the current process in the form
        of the :py:attr:`Process.template` attribute.
        """

        self._main_in_str = None
        """
        str: String used to specify the prefix of main input channel.
        """

        self._main_out_str = None
        """
        str: String used to specify the prefix of main output channel.
        """

        self._input_channel = None
        """
        str: Place holder of the main input channel for the current process.
        This attribute can change dynamically depending on the forks and
        secondary channels in the final pipeline.
        """

        self._output_channel = None
        """
        str: Place holder of the main output channel for the current process.
        This attribute can change dynamically depending on the forks and
        secondary channels in the final pipeline.
        """

        self._set_main_channel_name(ptype)

        self.link_start = [self._main_out_str]
        """
        list: List of strings with the starting points for secondary channels.
        When building the pipeline, these strings will be matched with equal
        strings in the :py:attr:`link_end` attribute of other Processes.
        """

        self.link_end = []
        """
        list: List of dictionaries containing the a string of the ending point
        for a secondary channel. Each dictionary should contain at least
        two key/vals: {"link": <link string>, "alias":<string for template>}
        """

        self.status_channels = ["STATUS"]
        self.status_strs = []
        """
        str: Name of the status channel for the current process. These strings
        will be provided to the StatusCompiler process to collect and
        compile status reports
        """

        self.forks = []
        """
        list: List of strings with the literal definition of the forks for
        the current process, ready to be added to the template string.
        """

        self._context = None
        """
        dict: Dictionary with the keyword placeholders for the string template
        of the current process.
        """

    def _set_template(self, template):
        """Sets the path to the appropriate jinja template file

        When a Process instance is initialized, this method will fetch
        the location of the appropriate template file, based on the
        ``template`` argument. It will raise an exception is the template
        file is not found. Otherwise, it will set the
        :py:attr:`Process.template_path` attribute.
        """

        # Set template directory
        tpl_dir = join(dirname(abspath(__file__)), "templates")

        # Set template file path
        tpl_path = join(tpl_dir, template + ".nf")

        if not os.path.exists(tpl_path):
            raise Exception("Template {} does not exist".format(tpl_path))

        self._template_path = join(tpl_dir, template + ".nf")

    def _set_main_channel_name(self, ptype):
        """

        Returns
        -------

        """

        if ptype == "pre_assembly":
            self._main_in_str = self._main_out_str = "MAIN_fq"
        elif ptype == "assembly":
            self._main_in_str = "MAIN_fq"
            self._main_out_str = "MAIN_assembly"
        else:
            self._main_in_str = self._main_out_str = "MAIN_assembly"

    @staticmethod
    def render(template, context):
        """Wrapper to the jinja2 render method from a template file

        Parameters
        ----------
        template : str
            Path to template file.
        context : dict
            Dictionary with kwargs context to populate the template
        """

        path, filename = os.path.split(template)

        return jinja2.Environment(
            loader=jinja2.FileSystemLoader(path or './')
        ).get_template(filename).render(context)

    @property
    def template_str(self):

        if not self._context:
            raise Exception("Channels must be setup first using the "
                            "set_channels method")

        logger.debug("Setting context for template {}: {}".format(
            self.template, self._context
        ))

        x = self.render(self._template_path, self._context)

        return x

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

        self._input_channel = "{}_{}".format(self._main_in_str, self.pid)

        for i in self.status_channels:
            self.status_strs.append("{}_{}".format(i, self.pid))

        if self.output_type:
            self._output_channel = "{}_{}".format(self._main_out_str,
                                                  self.pid + 1)

        self._context = {**kwargs, **{"input_channel": self._input_channel,
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

        logger.debug("Setting secondary channel for source '{}': {}".format(
            source, channel_list))

        # Handle the case where the main channel is forked
        if source.startswith("MAIN"):
            # Update previous output_channel to prevent overlap with
            # subsequent main channels. This is done by adding a "_" at the
            # beginning of the channel name
            self._context["output_channel"] = "_{}".format(
                self._output_channel)
            # Set source to modified output channel
            source = self._context["output_channel"]
            # Add the next first main channel to the channel_list
            channel_list.append(self._output_channel)
        # Handle forks from non main channels
        else:
            source = "{}_{}".format(source, self.pid)

        # When there is only one channel to fork into, use the 'set' operator
        # instead of 'into'
        if len(channel_list) == 1:
            self.forks.append("\n{}.set{{ {} }}\n".format(source,
                                                           channel_list[0]))
        else:
            self.forks.append("\n{}.into{{ {} }}\n".format(
                source, ";".join(channel_list)))

        logger.debug("Setting forks attribute to: {}".format(self.forks))

        self._context = {**self._context, **{"forks": "\n".join(self.forks)}}


class Status(Process):

    def __init__(self, **kwargs):

        super().__init__(**kwargs)

    def set_status_channels(self, channel_list):

        if len(channel_list) == 1:
            logger.debug("Setting only one status channel: {}".format(
                channel_list[0]))
            self._context = {"status_channels": channel_list[0]}

        else:
            first_status = channel_list[0]
            lst = ",".join(channel_list[1:])

            s = "{}.mix({})".format(first_status, lst)

            logger.debug("Status channel string: {}".format(s))

            self._context = {"status_channels": s}


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


class FastqcTrimmomatic(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="pre_assembly",
                         **kwargs)

        self.input_type = "fastq"
        self.output_type = "fastq"

        self.link_end.append({"link": "SIDE_phred", "alias": "SIDE_phred"})

        self.status_channels = ["STATUS_fastqc", "STATUS_trim"]


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

        self.status_channels = ["STATUS_am", "STATUS_amp"]

        self.link_end.append({"link": "MAIN_fq", "alias": "_MAIN_assembly"})


class Pilon(Process):

    def __init__(self, **kwargs):

        super().__init__(ptype="post_assembly",
                         **kwargs)

        self.input_type = "assembly"
        self.output_type = "assembly"

        self.dependencies = ["assembly_mapping"]


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


class StatusCompiler(Status):

    def __init__(self, **kwargs):

        super().__init__(ptype="status",
                         **kwargs)

        self.link_start = None

