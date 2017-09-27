#!/usr/bin/nextflow

// SETTING CHANNELS //
nsamples = file(params.fastq_files).size()
// Channel for FastQ files
fastq_raw = Channel.fromFilePairs(params.fastq_files)

// Channel for expected genome size
genome_size = Channel
                .value(params.genome_size)
// Channel for minimum coverage threshold
min_coverage = Channel
                .value(params.min_coverage)
// Channel for adapters file
adapters = Channel
                .value("None")
// Channels for Trimmomatic options
trimmomatic_opts = Channel
                .value([params.trim_sliding_window,
                        params.trim_leading,
                        params.trim_trailing,
                        params.tim_min_length])

spades_opts = Channel
                .value([params.spades_min_coverage,
                        params.spades_min_kmer_coverage])
spades_kmers = Channel
                .value(params.spades_kmers)

/** integrity_coverage
This process will check the integrity, encoding and get the estimated
coverage for each FastQ pair
*/
process integrity_coverage {

    tag { fastq_id }

	input:
	set fastq_id, file(fastq_pair) from fastq_raw
	val gsize from genome_size
	val cov from min_coverage

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

// Corrupted samples have the 2nd value with 'Corrupt'
integrity_processed.choice(corrupted, sample_ok) {
    a -> a[2].text == "Corrupt" ? 0 : 1
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

sample_listen.ifEmpty{ exit 1, "No samples left after checking FastQ integrity and estimating coverage. Exiting." }

/** report_coverage
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_initial.csv'
*/
process report_coverage {

    tag { report }
    publishDir 'reports/coverage/'

    input:
    file(report) from cov_report.filter{ it.text != "Corrupt" }.collect()

    output:
    file 'estimated_coverage_initial.csv'

    """
    echo 'Sample,Estimated coverage,Test' > estimated_coverage_initial.csv
    cat $report >> estimated_coverage_initial.csv
    """
}

/** report_corrupt
This process will report the corrupted samples and write the results to
'reports/corrupted/corrupted_samples.txt'
*/
process report_corrupt {

    tag { fastq_id }
    publishDir 'reports/corrupted/'

    input:
    val fastq_id from corrupted.map{it[0]}

    output:
    file 'corrupted_samples.txt'

    """
    echo $fastq_id >> corrupted_samples.txt
    """

}

/** fastqc
This process will perform the fastQC analysis for each sample. In this run,
the output files of FastQC are sent to the output channel
*/
process fastqc {

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from sample_good
    val ad from adapters

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into fastqc_listen, fastqc_processed
    file "fastq_status" into fastqc_status

    script:
    template "fastqc.py"
}

fastqc_listen.ifEmpty{ exit 1, "No samples left after running FastQC. Exiting." }

/** fastqc_report
This process will parse the result files from a FastQC analyses and output
the optimal_trim information for Trimmomatic
*/
process fastqc_report {

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from fastqc_processed

    output:
    set fastq_id, file(fastq_pair), 'fastqc_health', 'optimal_trim' into fastqc_trim

    script:
    template "fastqc_report.py"

}

// Triage of samples with bad health
fail_fastqc_report = Channel.create()
pass_fastqc_report = Channel.create()

fastqc_trim.choice(fail_fastqc_report, pass_fastqc_report) {
    a -> a[2].text == "pass" ? 0 : 1
}


/** trimmomatic
This process will execute trimmomatic
*/
process trimmomatic {

    tag { fastq_id }
    cpus 1

    input:
    set fastq_id, file(fastq_pair), trim_range, phred from pass_fastqc_report.phase(sample_phred).map{ [it[0][0], it[0][1], file(it[0][3]).text, it[1][1]] }
    val opts from trimmomatic_opts

    output:
    set fastq_id, "${fastq_id}_*P*" optional true into trimmomatic_listen, trimmomatic_processed
    file "trimmomatic_status" into trimmomatic_status

   script:
   template "trimmomatic.py"

}

trimmomatic_listen.ifEmpty{ exit 1, "No samples left after running Trimmomatic. Exiting." }

process integrity_coverage_2 {

    tag { fastq_id }

	input:
	set fastq_id, file(fastq_pair) from trimmomatic_processed
	val gsize from genome_size
	val cov from min_coverage

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

sample_listen_2.ifEmpty{ exit 1, "No samples left after running second estimated coverage. Exiting." }


/** report_coverage_2
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_second.csv'
*/
process report_coverage_2 {

    tag { report }
    publishDir 'reports/coverage/'

    input:
    file(report) from cov_report_2.filter{ it.text != "Corrupt" }.collect()

    output:
    file 'estimated_coverage_second.csv'

    """
    echo 'Sample,Estimated coverage,Test' > estimated_coverage_second.csv
    cat $report >> estimated_coverage_second.csv
    """
}


/** fastqc_2
This process will perform the second fastQC analysis for each sample.
In this run, the output files of FastQC are sent to the output channel
*/
process fastqc {

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from sample_good_2
    val ad from adapters

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into fastqc_listen_2, fastqc_processed_2
    file "fastq_status" into fastqc_status_2

    script:
    template "fastqc.py"
}

fastqc_listen_2.ifEmpty{ exit 1, "No samples left after running FastQC. Exiting." }


process spades {

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair), max_len from fastqc_processed_2.phase(sample_max_len).map{ [it[0][0], it[0][1], file(it[1][1]).text] }
    val opts from spades_opts
    val kmers from spades_kmers

    output:
    set fastq_pair, file('contigs.fasta') optional true into spades_listen, spades_processed
    file "spades_status" into spades_status

    script:
    template "spades.py"

}


spades_listen.ifEmpty{ exit 1, "No samples left after running Spades. Exiting." }

