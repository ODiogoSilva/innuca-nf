#!/usr/bin/nextflow

import Helper
import CheckParams

// Pipeline version
if (workflow.commitId){
    version = "0.1 $workflow.revision"
} else {
    version = "0.1 (local version)"
}

params.help = false
if (params.help){
    Help.print_help(version, params)
    exit 0
}

CheckParams.check(params)

nsamples = file(params.fastq).size()
Help.start_info(version, nsamples, "$workflow.start", "$workflow.profile")

// CHANNEL NOMENCLATURE //
// IN_* : Input channels, created from params options
// MAIN_* : The main pipeline channel. Must always contain a fastq_id and the sequence data
// LOG_* : Logging or reporting channels that are meant for terminal processes
// SIDE_* : Channels that are meant to merge with the MAIN channel downstream in the pipeline
// STATUS_* : Contain only the report of the status for any given process

// SETTING CHANNELS //
// GENERAL PARAMS //

// Channel for FastQ files
IN_fastq_raw = Channel.fromFilePairs(params.fastq)
// Channel for expected genome size
IN_genome_size = Channel
                .value(params.genomeSize)
// Channel for minimum coverage threshold
IN_min_coverage = Channel
                .value(params.minCoverage)

// FASTQC CHANNELS //
// Channel for adapters file
IN_adapters = Channel
                .value(params.adapters)

// TRIMMOMATIC CHANNELS //
IN_trimmomatic_opts = Channel
                .value([params.trimSlidingWindow,
                        params.trimLeading,
                        params.trimTrailing,
                        params.trimMinLength])

// SPADES CHANNELS //
IN_spades_opts = Channel
                .value([params.spadesMinCoverage,
                        params.spadesMinKmerCoverage])
IN_spades_kmers = Channel
                .value(params.spadesKmers)

IN_process_spades_opts = Channel
                .value([params.spadesMinContigLen,
                        params.spadesMinKmerCoverage])

// ASSEMBLY MAPPING CHANNELS //
IN_assembly_mapping_opts = Channel
                .value(params.minAssemblyCoverage)

/** INTEGRITY_COVERAGE - MAIN
This process will check the integrity, encoding and get the estimated
coverage for each FastQ pair. Corrupted FastQ files will also be detected
and filtered here.
*/
process integrity_coverage {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 1"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1"
    }

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1

	input:
	set fastq_id, file(fastq_pair) from IN_fastq_raw
	val gsize from IN_genome_size
	val cov from IN_min_coverage
	// This channel is for the custom options of the integrity_coverage.py
	// script. See the script's documentation for more information.
	val opts from Channel.value('')

	output:
	set fastq_id,
	    file(fastq_pair),
	    file('*_encoding'),
	    file('*_phred'),
	    file('*_coverage') into MAIN_integrity
	file('*_report') into LOG_report_coverage1

	script:
	template "integrity_coverage.py"

}

// TRIAGE OF CORRUPTED SAMPLES
LOG_corrupted = Channel.create()
MAIN_PreCoverageCheck = Channel.create()
// Corrupted samples have the 2nd value with 'corrupt'
MAIN_integrity.choice(LOG_corrupted, MAIN_PreCoverageCheck) {
    a -> a[2].text == "corrupt" ? 0 : 1
}

// TRIAGE OF LOW COVERAGE SAMPLES
MAIN_fastqc_in = Channel.create()
SIDE_phred = Channel.create()

MAIN_PreCoverageCheck
// Low coverage samples have the 4th value of the Channel with 'fail'
    .filter{ it[4].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'SIDE_phred'
    .separate(MAIN_fastqc_in, SIDE_phred){
        a -> [ [a[0], a[1]], [a[0], a[3].text] ]
    }

/** REPORT_COVERAGE - PLUG-IN
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_initial.csv'
*/
process report_coverage {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage/'

    input:
    file(report) from LOG_report_coverage1.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_initial.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_initial.csv
    cat $report >> estimated_coverage_initial.csv
    """
}

/** REPORT_CORRUPT - PLUG-IN
This process will report the corrupted samples and write the results to
'reports/corrupted/corrupted_samples.txt'
*/
process report_corrupt {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/corrupted/'

    input:
    val fastq_id from LOG_corrupted.collect{it[0]}

    output:
    file 'corrupted_samples.txt'

    """
    echo ${fastq_id.join(",")} | tr "," "\n" >> corrupted_samples.txt
    """

}

/** FASTQC - MAIN
This process will perform the fastQC analysis for each sample. In this run,
the output files (summary and data) of FastQC are sent to the output channel
as pair_1* and pair_2* files.
*/
process fastqc {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 2"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 2"
    }

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from MAIN_fastqc_in
    val ad from Channel.value('None')

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into MAIN_fastqc_out
    set fastq_id, val("fastqc"), file(".status") into STATUS_fastqc

    when:
    params.stopAt != "fastqc"

    script:
    template "fastqc.py"
}

/** FASTQC_REPORT - MAIN
This process will parse the result files from a FastQC analyses and output
the optimal_trim information for Trimmomatic
*/
process fastqc_report {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 3"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 3"
    }

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc/run_1/', pattern: '*summary.txt', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out
    val opts from Channel.value("--ignore-tests")

    output:
    set fastq_id, file(fastq_pair), '.status', 'optimal_trim' into MAIN_fastqc_trim
    file '*_trim_report' into LOG_trim
    file "*_status_report" into LOG_fastqc_report
    file "${fastq_id}_*_summary.txt" optional true

    script:
    template "fastqc_report.py"

}


/** TRIM_REPORT - PLUG-IN
This will collect the optimal trim points assessed by the fastqc_report
process and write the results of all samples in a single csv file
*/
process trim_report {

    publishDir 'reports/fastqc/', mode: 'copy'

    input:
    file trim from LOG_trim.collect()

    output:
    file "FastQC_trim_report.csv"

    """
    echo Sample,Trim begin, Trim end >> FastQC_trim_report.csv
    cat $trim >> FastQC_trim_report.csv
    """
}


process compile_fastqc_status {

    publishDir 'reports/fastqc/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report.collect()

    output:
    file 'FastQC_1run_report.csv'

    """
    echo Sample, Failed? >> FastQC_1run_report.csv
    cat $rep >> FastQC_1run_report.csv
    """

}


/** TRIMMOMATIC - MAIN
This process will execute trimmomatic. Currently, the main channel requires
information on the trim_range and phred score.
*/
process trimmomatic {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 4"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 4"
    }

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair), trim_range, phred from MAIN_fastqc_trim.phase(SIDE_phred).map{ [it[0][0], it[0][1], file(it[0][3]).text, it[1][1]] }
    val opts from IN_trimmomatic_opts

    output:
    set fastq_id, "${fastq_id}_*P*" optional true into MAIN_trimmomatic_out, SIDE_bowtie_in
    set fastq_id, val("trimmomatic"), file(".status") into STATUS_trimmomatic
    file '*_trimlog.txt' optional true into LOG_trimmomatic

    when:
    params.stopAt != "trimmomatic"

    script:
    template "trimmomatic.py"

}


process trimmomatic_report {

    publishDir 'reports/trimmomatic/'

    input:
    file log_files from LOG_trimmomatic.collect()

    output:
    file 'trimmomatic_report.csv'

    script:
    template "trimmomatic_report.py"

}

/** INTEGRITY_COVERAGE_2 - MAIN
This process will estimate the coverage of the processed FastQ files after
the trimmomatic trimming. Note that the encoding guessing is turned-of
by using the '-e' option in the opts variable.
*/
process integrity_coverage_2 {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 5"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 5"
    }

    tag { fastq_id }
    cpus 1

	input:
	set fastq_id, file(fastq_pair) from MAIN_trimmomatic_out
	val gsize from IN_genome_size
	val cov from IN_min_coverage
	// Use -e option for skipping encoding guess
	val opts from Channel.value('-e')

	output:
	set fastq_id,
	    file(fastq_pair),
	    file('*_coverage'),
	    file('*_max_len') into MAIN_integrity2
	file('*_report') into LOG_report_coverage2

	script:
	template "integrity_coverage.py"
}

// Checking for coverage again after trimmomatic trimming.
// Low coverage samples have the 2nd value of the Channel with 'fail'
MAIN_fastqc_in2 = Channel.create()
SIDE_max_len = Channel.create()

MAIN_integrity2
// Low coverage samples have the 2nd value of the Channel with 'fail'
    .filter{ it[2].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'SIDE_phred'
    .separate(MAIN_fastqc_in2, SIDE_max_len){
        a -> [ [a[0], a[1]], [a[0], a[3]]]
    }


/** REPORT_COVERAGE_2 - PLUG-IN
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_second.csv'
*/
process report_coverage_2 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage/'

    input:
    file(report) from LOG_report_coverage2.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_second.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_second.csv
    cat $report >> estimated_coverage_second.csv
    """
}


/** FASTQC_2 - MAIN
This process will perform the second fastQC analysis for each sample.
In this run, the output files of FastQC are sent to the output channel
*/
process fastqc2 {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 6"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 6"
    }

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from MAIN_fastqc_in2
    val ad from IN_adapters

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into MAIN_fastqc_out2
    set fastq_id, val("fastqc2"), file(".status") into STATUS_fastqc2

    when:
    params.stopAt != "fastqc2"

    script:
    template "fastqc.py"
}


process fastqc2_report {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 7"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 7"
    }

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc/run_2/', pattern: '*summary.txt', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out2
    val opts from Channel.value("")

    output:
    set fastq_id, file(fastq_pair), '.status' into MAIN_fastqc_report
    file "*_status_report" into LOG_fastqc_report2
    file "${fastq_id}_*_summary.txt" optional true

    script:
    template "fastqc_report.py"

}


process compile_fastqc_status2 {

    publishDir 'reports/fastqc/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report2.collect()

    output:
    file 'FastQC_2run_report.csv'

    """
    echo Sample, Failed? >> FastQC_2run_report.csv
    cat $rep >> FastQC_2run_report.csv
    """

}

MAIN_spades_in = Channel.create()

MAIN_fastqc_report
        .filter{ it[2].text == "pass" }
        .map{ [it[0], it[1]] }
        .into(MAIN_spades_in)

/** SPADES - MAIN
This process performs the FastQ assembly using SPAdes. Besides the FastQ
files, this process requires an estimate of the maximum contig len
(inferred in the integrity_coverage process), and user specified
options for SPAdes (see spades.py template).
*/
process spades {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 8"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 8"
    }

    tag { fastq_id }
    publishDir 'results/assembly/spades/', pattern: '*_spades.assembly.fasta', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), max_len from MAIN_spades_in.phase(SIDE_max_len).map{ [it[0][0], it[0][1], file(it[1][1]).text] }
    val opts from IN_spades_opts
    val kmers from IN_spades_kmers

    output:
    set fastq_id, file('*_spades.assembly.fasta') optional true into MAIN_spades_out, LOG_spades
    set fastq_id, val("spades"), file(".status") into STATUS_spades

    when:
    params.stopAt != "spades"

    script:
    template "spades.py"

}

/** SPADES_REPORT - PLUG IN
Plug-in process that provides an assembly report for each sample
*/
process spades_report {

    tag { fastq_id }

    input:
    set fastq_id, file(assembly) from LOG_spades

    output:
    file "*_assembly_report.csv" into LOG_spades_report

    script:
    template "assembly_report.py"
}


/** COMPILE_ASSEMBLY_REPORT - PLUG IN
Plug-in process that compiles the results of the spades_report process for
all samples
*/
process compile_spades_report {

    publishDir "reports/assembly/spades/", mode: 'copy'

    input:
    file(report) from LOG_spades_report.collect()

    output:
    file "spades_assembly_report.csv"

    """
    echo Sample,Number of contigs,Average contig size,N50,Total assembly length,GC content,Missing data > spades_assembly_report.csv
    cat $report >> spades_assembly_report.csv
    """
}


/** PROCESS_SPADES - MAIN
Processes and filters the SPAdes assembly according to user-specified options
*/
process process_spades {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 9"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 9"
    }

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1
    publishDir "reports/assembly/spades_filter", pattern: '*.report.csv', mode: 'copy'

    input:
    set fastq_id, file(assembly) from MAIN_spades_out
    val opts from IN_process_spades_opts
    val gsize from IN_genome_size

    output:
    set fastq_id, file('*.assembly.fasta') into MAIN_spades_filtered
    set fastq_id, val("process_spades"), file(".status") into STATUS_process_spades
    file '*.report.csv'

    when:
    params.stopAt != "process_spades"

    script:
    template "process_spades.py"

}

// Merge the bowtie_input channel from Trimmomatic with the processed FastQ
// files and the spades assembly file.
// The resulting channel will consist of:
// fastq_id, fastq_1, fastq_2, assembly_file
MAIN_am_in = Channel.create()
SIDE_bowtie_in
        .phase(MAIN_spades_filtered)
        .map{ [it[0][0], it[0][1][0], it[0][1][1], it[1][1]] }
        .into(MAIN_am_in)


/** ASSEMBLY_MAPPING - MAIN
Performs the mapping of pairs of FastQ reads into an assembled genome.
It outputs a table with the coverage estimates for each contig in the
assembly as well as the sorted BAM file.
*/
process assembly_mapping {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 10"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 10"
    }

    tag { fastq_id }
    echo false

    input:
    set fastq_id, file(fastq_1), file(fastq_2), file(assembly) from MAIN_am_in

    output:
    set fastq_id, file(assembly), 'coverages.tsv', 'sorted.bam', 'sorted.bam.bai' optional true into MAIN_am_out
    set fastq_id, val("assembly_mapping"), file(".status") into STATUS_am

    when:
    params.stopAt != "assembly_mapping"

    script:
    """
    {
        bowtie2-build --threads ${task.cpus} $assembly genome_index
        bowtie2 -q --very-sensitive-local --threads ${task.cpus} -x genome_index -1 $fastq_1 -2 $fastq_2 -S mapping.sam
        samtools sort -o sorted.bam -O bam -@ ${task.cpus} mapping.sam && rm *.sam
        samtools index sorted.bam
        parallel -j ${task.cpus} samtools depth -ar {} sorted.bam \\> {}.tab  ::: \$(grep ">" $assembly | cut -c 2-)
        # Insert 0 coverage count in empty files. See Issue #2
        find . -size 0 -print0 | xargs -0 -I{} sh -c 'echo -e 0"\t"0"\t"0 > "{}"'
        parallel -j ${task.cpus} echo -n {.} '"\t"' '&&' cut -f3 {} '|' paste -sd+ '|' bc >> coverages.tsv  ::: *.tab
        rm *.tab
        if [ -f "coverages.tsv" ]
        then
            echo pass > .status
        else
            echo fail > .status
        fi
    } || {
        echo fail > .status
    }
    """
}


/** PROCESS_ASSEMBLY_MAPPING -  MAIN
Processes the results from the assembly_mapping process and filters the
assembly contigs based on coverage and length thresholds.
*/
process process_assembly_mapping {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 11"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 11"
    }

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly), file(coverage), file(bam_file), file(bam_index) from MAIN_am_out
    val min_assembly_coverage from IN_assembly_mapping_opts
    val gsize from IN_genome_size

    output:
    set fastq_id, '*_filtered.assembly.fasta', 'filtered.bam', 'filtered.bam.bai' into MAIN_pilon_in
    set fastq_id, val("process_am"), file(".status") into STATUS_process_am

    script:
    template "process_assembly_mapping.py"

}


/** PILON - MAIN
Executes the pilon software on a given assembly, with the sorted BAM file
resulting from the mapping of the raw reads into the assembly.
*/
process pilon {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 12"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 12"
    }

    tag { fastq_id }
    echo false
    publishDir 'results/assembly/pilon/', mode: 'copy'

    input:
    set fastq_id, file(assembly), file(bam_file), file(bam_index) from MAIN_pilon_in

    output:
    set fastq_id, '*_polished.assembly.fasta' into MAIN_pilon_out, MAIN_pilon
    set fastq_id, val("pilon"), file(".status") into STATUS_pilon

    script:
    """
    {
        pilon_mem=${String.valueOf(task.memory).substring(0, String.valueOf(task.memory).length() - 1).replaceAll("\\s", "")}
        java -jar -Xms256m -Xmx\${pilon_mem} /NGStools/pilon-1.22.jar --genome $assembly --frags $bam_file --output ${fastq_id}_polished.assembly --changes --vcf --threads $task.cpus
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


// Post assembly processes that required an assembly file
MAIN_mlst_in = Channel.create()
MAIN_prokka_in = Channel.create()
MAIN_abricate_in = Channel.create()
MAIN_chewbbaca = Channel.create()
// For last assembly channel
MAIN_pilon_out.into{ MAIN_mlst_in;MAIN_prokka_in;MAIN_abricate_in;MAIN_chewbbaca }


process pilon_report {

    tag { fastq_id }

    input:
    set fastq_id, file(assembly) from MAIN_pilon

    output:
    file "*_assembly_report.csv" into MAIN_pilon_report

    script:
    template "assembly_report.py"

}


process compile_pilon_report {

    publishDir "reports/assembly/pilon/", mode: 'copy'

    input:
    file(report) from MAIN_pilon_report.collect()

    output:
    file "pilon_assembly_report.csv"

    """
    echo Sample,Number of contigs,Average contig size,N50,Total assembly length,GC content,Missing data > pilon_assembly_report.csv
    cat $report >> pilon_assembly_report.csv
    """
}


/** MLST - PLUG-IN
Executs MLST on a given assembly
*/
process mlst {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 13"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 13"
    }

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly) from MAIN_mlst_in

    output:
    file '*.mlst.txt' into MAIN_mlst_out
    set fastq_id, val("mlst"), file(".status") into STATUS_mlst

    when:
    params.mlstRun  == true && params.annotationRun

    script:
    """
    {
        mlst $assembly >> ${fastq_id}.mlst.txt
        echo pass > .status
    } || {
        echo fail > .status
    }
    """
}


process compile_mlst {

    publishDir "results/annotation/mlst/"

    input:
    file res from MAIN_mlst_out.collect()

    output:
    file "mlst_report.tsv"

    when:
    params.mlstRun == true && params.annotationRun

    script:
    """
    cat $res >> mlst_report.tsv
    """
}


process abricate {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 14"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 14"
    }

    tag { "${fastq_id} ${db}" }
    publishDir "results/annotation/abricate/${fastq_id}"

    input:
    set fastq_id, file(assembly) from MAIN_abricate_in
    each db from params.abricateDatabases

    output:
    file '*.tsv'
    set fastq_id, val("abricate_${db}"), file(".status") into STATUS_abricate

    when:
    params.abricateRun == true && params.annotationRun

    script:
    """
    {
        abricate --db $db $assembly > ${fastq_id}_abr_${db}.tsv
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


process prokka {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 15"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 15"
    }

    tag { fastq_id }
    publishDir "results/annotation/prokka/${fastq_id}"

    input:
    set fastq_id, file(assembly) from MAIN_prokka_in

    output:
    file "${fastq_id}/*"
    set fastq_id, val("prokka"), file(".status") into STATUS_prokka

    when:
    params.prokkaRun == true && params.annotationRun

    script:
    """
    {
        prokka --outdir $fastq_id --cpus $task.cpus --centre UMMI --compliant \
               --increment 10 $assembly
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


process chewbbaca {

    // Send POST request to platform
    if ( params.platformHTTP != null ) {
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 16"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 16"
    }

    maxForks 1
    tag { fastq_id }
    echo false
    scratch true
    publishDir "results/chewbbaca/${fastq_id}"

    input:
    set fastq_id, file(assembly) from MAIN_chewbbaca
    each file(schema) from Channel.fromPath(params.schema_path)

    output:
    file 'chew_results'
    set fastq_id, val("chewbbaca"), file(".status") into STATUS_chewbbaca

    when:
    params.chewbbacaRun == true

    script:
    """
    {
        echo $assembly >> input_file.txt
        chewBBACA.py AlleleCall -i input_file.txt -g $schema -o chew_results --json --cpu $task.cpus -t "Streptococcus agalactiae"
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


// LISTENER PROCESSES
// The next set of processes are intended to be of general use to several
// processes for reporting/status purposes. They basically listen to
// channels from arbitrary process during the pipeline execution
// for those purposes. Therefore, they must be defined at the end.


/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { fastq_id }

    input:
    set fastq_id, task_name, status from STATUS_fastqc.mix(STATUS_trimmomatic,
                                                           STATUS_fastqc2,
                                                           STATUS_spades,
                                                           STATUS_process_spades,
                                                           STATUS_am,
                                                           STATUS_process_am,
                                                           STATUS_pilon,
                                                           STATUS_mlst,
                                                           STATUS_abricate,
                                                           STATUS_prokka,
                                                           STATUS_chewbbaca)

    output:
    file 'status_*' into master_status

    """
    echo $fastq_id, $task_name, \$(cat $status) > status_${fastq_id}_${task_name}
    """
}

process compile_status {

    publishDir 'reports/status'

    input:
    file status from master_status.collect()

    output:
    file 'master_status.csv'

    """
    cat $status >> master_status.csv
    """
}

workflow.onComplete{
    Help.complete_info(workflow)
}