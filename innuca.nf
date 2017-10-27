#!/usr/bin/nextflow

import Helper
import CheckParams

// Pipeline version
version = "0.1"

params.help = false
if (params.help){
    Help.print_help(version, params)
    exit 0
}

CheckParams.check(params)

nsamples = file(params.fastq).size()
Help.start_info(version, nsamples, "$workflow.start", "$workflow.profile")

// SETTING CHANNELS //
// GENERAL PARAMS //

// Channel for FastQ files
fastq_raw = Channel.fromFilePairs(params.fastq)
// Channel for expected genome size
genome_size = Channel
                .value(params.genomeSize)
// Channel for minimum coverage threshold
min_coverage = Channel
                .value(params.minCoverage)

// FASTQC CHANNELS //
// Channel for adapters file
adapters = Channel
                .value(params.adapters)

// TRIMMOMATIC CHANNELS //
trimmomatic_opts = Channel
                .value([params.trimSlidingWindow,
                        params.trimLeading,
                        params.trimTrailing,
                        params.trimMinLength])

// SPADES CHANNELS //
spades_opts = Channel
                .value([params.spadesMinCoverage,
                        params.spadesMinKmerCoverage])
spades_kmers = Channel
                .value(params.spadesKmers)

process_spades_opts = Channel
                .value([params.spadesMinContigLen,
                        params.spadesMinKmerCoverage])

// ASSEMBLY MAPPING CHANNELS //
assembly_mapping_opts = Channel
                .value(params.minAssemblyCoverage)

/** INTEGRITY_COVERAGE - MAIN
This process will check the integrity, encoding and get the estimated
coverage for each FastQ pair. Corrupted FastQ files will also be detected
and filtered here.
*/
process integrity_coverage {

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1

	input:
	set fastq_id, file(fastq_pair) from fastq_raw
	val gsize from genome_size
	val cov from min_coverage
	// This channel is for the custom options of the integrity_coverage.py
	// script. See the script's documentation for more information.
	val opts from Channel.value('')

	output:
	set fastq_id,
	    file(fastq_pair),
	    file('*_encoding'),
	    file('*_phred'),
	    file('*_coverage') into integrity_processed
	file('*_report') into cov_report

	script:
	template "integrity_coverage.py"

}

// TRIAGE OF CORRUPTED SAMPLES
corrupted = Channel.create()
sample_ok = Channel.create()
// Corrupted samples have the 2nd value with 'corrupt'
integrity_processed.choice(corrupted, sample_ok) {
    a -> a[2].text == "corrupt" ? 0 : 1
}

// TRIAGE OF LOW COVERAGE SAMPLES
sample_good = Channel.create()
sample_listen = Channel.create()
sample_phred = Channel.create()

sample_ok
// Low coverage samples have the 4th value of the Channel with 'fail'
    .filter{ it[4].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'sample_phred'
    .separate(sample_good, sample_listen, sample_phred){
        a -> [ [a[0], a[1]], [a[0], a[1]], [a[0], a[3].text] ]
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
    file(report) from cov_report.filter{ it.text != "corrupt" }.collect()

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
    val fastq_id from corrupted.collect{it[0]}

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

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from sample_good
    val ad from Channel.value('None')

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into fastqc_processed
    set fastq_id, val("fastqc"), file("fastq_status") into fastqc_status

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

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc/', pattern: '*summary.txt', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from fastqc_processed

    output:
    set fastq_id, file(fastq_pair), 'fastqc_health', 'optimal_trim' into fastqc_trim
    file 'report' into trim_rep
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
    file trim from trim_rep.collect()

    output:
    file "FastQC_trim_report.csv"

    """
    echo Sample,Trim begin, Trim end >> FastQC_trim_report.csv
    cat $trim >> FastQC_trim_report.csv
    """
}


// Triage of samples with bad health according to FastQC report
fail_fastqc_report = Channel.create()
pass_fastqc_report = Channel.create()

fastqc_trim.choice(fail_fastqc_report, pass_fastqc_report) {
    a -> a[2].text == "pass" ? 0 : 1
}


/** TRIMMOMATIC - MAIN
This process will execute trimmomatic. Currently, the main channel requires
information on the trim_range and phred score.
*/
process trimmomatic {

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair), trim_range, phred from pass_fastqc_report.phase(sample_phred).map{ [it[0][0], it[0][1], file(it[0][3]).text, it[1][1]] }
    val opts from trimmomatic_opts

    output:
    set fastq_id, "${fastq_id}_*P*" optional true into trimmomatic_processed, bowtie_input
    set fastq_id, val("trimmomatic"), file("trimmomatic_status") into trimmomatic_status
    file '*_trimlog.txt' optional true into trimmomatic_log

    when:
    params.stopAt != "trimmomatic"

    script:
    template "trimmomatic.py"

}


process trimmomatic_report {

    publishDir 'reports/trimmomatic/'

    input:
    file log_files from trimmomatic_log.collect()

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

    tag { fastq_id }
    cpus 1

	input:
	set fastq_id, file(fastq_pair) from trimmomatic_processed
	val gsize from genome_size
	val cov from min_coverage
	// Use -e option for skipping encoding guess
	val opts from Channel.value('-e')

	output:
	set fastq_id,
	    file(fastq_pair),
	    file('*_coverage'),
	    file('*_max_len') into integrity_processed_2
	file('*_report') into cov_report_2

	script:
	template "integrity_coverage.py"
}

// Checking for coverage again after trimmomatic trimming.
// Low coverage samples have the 2nd value of the Channel with 'fail'
sample_good_2 = Channel.create()
sample_max_len = Channel.create()
sample_listen_2 = Channel.create()

integrity_processed_2
// Low coverage samples have the 2nd value of the Channel with 'fail'
    .filter{ it[2].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'sample_phred'
    .separate(sample_good_2, sample_listen_2, sample_max_len){
        a -> [ [a[0], a[1]], [a[0], a[1]], [a[0], a[3]]]
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
    file(report) from cov_report_2.filter{ it.text != "corrupt" }.collect()

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
process fastqc {

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from sample_good_2
    val ad from adapters

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into fastqc_processed_2
    set fastq_id, val("fastqc2"), file("fastq_status") into fastqc_status_2

    when:
    params.stopAt != "fastqc2"

    script:
    template "fastqc.py"
}


/** SPADES - MAIN
This process performs the FastQ assembly using SPAdes. Besides the FastQ
files, this process requires an estimate of the maximum contig len
(inferred in the integrity_coverage process), and user specified
options for SPAdes (see spades.py template).
*/
process spades {

    tag { fastq_id }
    publishDir 'results/assemblies/spades/', pattern: '*_spades.assembly.fasta', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), max_len from fastqc_processed_2.phase(sample_max_len).map{ [it[0][0], it[0][1], file(it[1][1]).text] }
    val opts from spades_opts
    val kmers from spades_kmers

    output:
    set fastq_id, file('*_spades.assembly.fasta') optional true into spades_processed, s_report
    set fastq_id, val("spades"), file("spades_status") into spades_status

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
    set fastq_id, file(assembly) from s_report

    output:
    file "*_assembly_report.csv" into sm_report

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
    file(report) from sm_report.collect()

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

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1
    publishDir "reports/assembly/spades_filter", pattern: '*.report.csv', mode: 'copy'

    input:
    set fastq_id, file(assembly) from spades_processed
    val opts from process_spades_opts
    val gsize from genome_size

    output:
    set fastq_id, file('*.assembly.fasta') into spades_assembly
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
assembly_mapping_input = Channel.create()
bowtie_input
        .phase(spades_assembly)
        .map{ [it[0][0], it[0][1][0], it[0][1][1], it[1][1]] }
        .into(assembly_mapping_input)


/** ASSEMBLY_MAPPING - MAIN
Performs the mapping of pairs of FastQ reads into an assembled genome.
It outputs a table with the coverage estimates for each contig in the
assembly as well as the sorted BAM file.
*/
process assembly_mapping {

    tag { fastq_id }
    echo false

    input:
    set fastq_id, file(fastq_1), file(fastq_2), file(assembly) from assembly_mapping_input

    output:
    set fastq_id, file(assembly), 'coverages.tsv', 'sorted.bam', 'sorted.bam.bai' optional true into mapping_coverage
    set fastq_id, val("assembly_mapping"), file("assembly_mapping_status") into assembly_mapping_status

    when:
    params.stopAt != "assembly_mapping"

    script:
    """
    bowtie2-build --threads ${task.cpus} $assembly genome_index
    bowtie2 -q --very-sensitive-local --threads ${task.cpus} -x genome_index -1 $fastq_1 -2 $fastq_2 -S mapping.sam
    samtools sort -o sorted.bam -O bam -@ ${task.cpus} mapping.sam && rm *.sam
    samtools index sorted.bam
    parallel -j ${task.cpus} samtools depth -ar {} sorted.bam \\> {}.tab  ::: \$(grep ">" $assembly | cut -c 2-)
    parallel -j ${task.cpus} echo -n {.} '"\t"' '&&' cut -f3 {} '|' paste -sd+ '|' bc >> coverages.tsv  ::: *.tab
    rm *.tab
    if [ -f "coverages.tsv" ]
    then
        echo pass > assembly_mapping_status
    else
        echo fail > assembly_mapping_status
    fi
    """
}


/** PROCESS_ASSEMBLY_MAPPING -  MAIN
Processes the results from the assembly_mapping process and filters the
assembly contigs based on coverage and length thresholds.
*/
process process_assembly_mapping {

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly), file(coverage), file(bam_file), file(bam_index) from mapping_coverage
    val min_assembly_coverage from assembly_mapping_opts
    val gsize from genome_size

    output:
    set fastq_id, '*_filtered.assembly.fasta', 'filtered.bam', 'filtered.bam.bai' into processed_assembly_mapping

    script:
    template "process_assembly_mapping.py"

}


/** PILON - MAIN
Executes the pilon software on a given assembly, with the sorted BAM file
resulting from the mapping of the raw reads into the assembly.
*/
process pilon {

    tag { fastq_id }
    echo false
    publishDir 'results/assemblies/pilon/', mode: 'copy'

    input:
    set fastq_id, file(assembly), file(bam_file), file(bam_index) from processed_assembly_mapping

    output:
    set fastq_id, '*_polished.assembly.fasta' into pilon_processed, p_report


    """
    java -jar /NGStools/pilon-1.22.jar --genome $assembly --frags $bam_file --output ${fastq_id}_polished.assembly --changes --vcf --threads $task.cpus
    """

}


// Post assembly processes that required an assembly file
mlst_input = Channel.create()
prokka_input = Channel.create()
abricate_input = Channel.create()
// For last assembly channel
pilon_processed.into{ mlst_input;prokka_input;abricate_input }


process pilon_report {

    tag { fastq_id }

    input:
    set fastq_id, file(assembly) from p_report

    output:
    file "*_assembly_report.csv" into pm_report

    script:
    template "assembly_report.py"

}


process compile_spades_report {

    publishDir "reports/assembly/pilon/", mode: 'copy'

    input:
    file(report) from pm_report.collect()

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

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly) from mlst_input

    output:
    file '*.mlst.txt' into mlst_result

    when:
    params.mlstRun  == true

    script:
    """
    mlst $assembly >> ${fastq_id}.mlst.txt
    """
}


process compile_mlst {

    publishDir "results/mlst/"

    input:
    file res from mlst_result.collect()

    output:
    file "mlst_report.tsv"

    when:
    params.mlstRun == true

    script:
    """
    cat $res >> mlst_report.tsv
    """
}


process abricate {

    tag { fastq_id }
    publishDir "results/abricate/${fastq_id}"

    input:
    set fastq_id, file(assembly) from abricate_input
    each db from params.abricateDatabases

    output:
    file '*.tsv'

    when:
    params.abricateRun == true

    script:
    """
    abricate --db $db $assembly > ${fastq_id}_abr_${db}.tsv
    """

}


process prokka {

    tag { fastq_id }

    input:
    set fastq_id, file(assembly) from prokka_input

    """
    prokka --outdir $fastq_id --cpus $task.cpus --centre UMMI --compliant \
           --locustag ${fastq_id}p --increment 10 $assembly
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
    set fastq_id, task_name, status from fastqc_status.mix(trimmomatic_status,
                                                           fastqc_status_2,
                                                           spades_status,
                                                           assembly_mapping_status)

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