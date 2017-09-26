#!/usr/bin/nextflow

// PARAMETERS //
// These can be defined when executing the script
params.fastq_files = "data/*_{1,2}.*"
params.genome_size = 2.1
params.min_coverage = 15
// Trimmomatic
params.trim_sliding_window = '5:20'
params.trim_leading = '3'
params.trim_trailing = '3'
params.tim_min_length = '55'

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
	    file('*_coverage') into integrity
	file('*_report') into cov_report

	script:
	template "integrity_coverage.py"

}

// TRIAGE OF CORRUPTED SAMPLES
corrupted = Channel.create()
sample_ok = Channel.create()

// Corrupted samples have the 2nd value with 'Corrupt'
integrity.choice(corrupted, sample_ok) {
    a -> a[2].text == "Corrupt" ? 0 : 1
}

// TRIAGE OF LOW COVERAGE SAMPLES
sample_good = Channel.create()
sample_phred = Channel.create()

sample_ok
// Low coverage samples have the 4th value of the Channel with 'None'
    .filter{ it[4].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'sample_phred'
    .separate(sample_good, sample_phred){
        a -> [ [a[0], a[1]], [a[0], a[3].text] ]
    }


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
    container 'odiogosilva/fastqc:0.11.5'

    input:
    set fastq_id, file(fastq_pair) from sample_good
    val ad from adapters

    output:
    set fastq_id, file(fastq_pair), "fastq_status", file('pair_1*'), file('pair_2*') into fastqc_processed

    script:
    template "fastqc.py"

}

// TRIAGE OF UNSUCCESSFUL FASTQC SAMPLE ANALYSES
bad_fastqc = Channel.create()
ok_fastqc = Channel.create()

// Corrupted samples have the 2nd value with 'Corrupt'
fastqc_processed.choice(bad_fastqc, ok_fastqc) {
    a -> a[2].text == "pass" ? 1 : 0
}

// Prepare channel for processing of FastQC results
fastqc_results = Channel.create()

ok_fastqc
        .map{ [it[0], it[1], it[3], it[4]] }
        .into(fastqc_results)

/** fastqc_report
This process will parse the result files from a FastQC analyses and output
the optimal_trim information for Trimmomatic
*/
process fastqc_report {

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from fastqc_results

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
    container 'odiogosilva/trimmomatic'

    input:
    set fastq_id, file(fastq_pair), trim_range, phred from pass_fastqc_report.phase(sample_phred).map{ [it[0][0], it[0][1], file(it[0][3]).text, it[1][1]] }
    val opts from trimmomatic_opts

    output:
    set fastq_id, file(fastq_pair), "trimmomatic_status" into trimmomatic_processed

   script:
   template "trimmomatic.py"

}
