#!/usr/bin/nextflow

// SETTING CHANNELS //
// GENERAL PARAMS //
nsamples = file(params.fastq_files).size()
// Channel for FastQ files
fastq_raw = Channel.fromFilePairs(params.fastq_files)
// Channel for expected genome size
genome_size = Channel
                .value(params.genome_size)
// Channel for minimum coverage threshold
min_coverage = Channel
                .value(params.min_coverage)

// FASTQC CHANNELS //
// Channel for adapters file
adapters = Channel
                .value(params.adapters)

// TRIMMOMATIC CHANNELS //
trimmomatic_opts = Channel
                .value([params.trim_sliding_window,
                        params.trim_leading,
                        params.trim_trailing,
                        params.tim_min_length])

// SPADES CHANNELS //
spades_opts = Channel
                .value([params.spades_min_coverage,
                        params.spades_min_kmer_coverage])
spades_kmers = Channel
                .value(params.spades_kmers)

process_spades_opts = Channel
                .value([params.spades_min_contig_len,
                        params.spades_min_kmer_coverage])

// ASSEMBLY MAPPING CHANNELS //
assembly_mapping_opts = Channel
                .value(params.min_assembly_coverage)

/** integrity_coverage
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

sample_listen.ifEmpty{ exit 1, "No samples left after checking FastQ integrity and estimating coverage. Exiting." }

/** report_coverage
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

/** report_corrupt
This process will report the corrupted samples and write the results to
'reports/corrupted/corrupted_samples.txt'
*/
process report_corrupt {

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1
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
the output files (summary and data) of FastQC are sent to the output channel
as pair_1* and pair_2* files.
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

/** trim_report
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


/** trimmomatic
This process will execute trimmomatic. Currently, the main channel requires
information on the trim_range and phred score.
*/
process trimmomatic {

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair), trim_range, phred from pass_fastqc_report.phase(sample_phred).map{ [it[0][0], it[0][1], file(it[0][3]).text, it[1][1]] }
    val opts from trimmomatic_opts

    output:
    set fastq_id, "${fastq_id}_*P*" optional true into trimmomatic_listen, trimmomatic_processed, bowtie_input
    file "trimmomatic_status" into trimmomatic_status

    script:
    template "trimmomatic.py"

}

trimmomatic_listen.ifEmpty{ exit 1, "No samples left after running Trimmomatic. Exiting." }

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

sample_listen_2.ifEmpty{ exit 1, "No samples left after running second estimated coverage. Exiting." }


/** report_coverage_2
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
    set fastq_id, file('contigs.fasta') optional true into spades_listen, spades_processed, s_report
    file "spades_status" into spades_status

    script:
    template "spades.py"

}


process spades_report {

    tag { fastq_id }

    input:
    set fastq_id, file(assembly) from s_report

    output:
    set val('spades'), "*_assembly_report.csv" into ma_report

    script:
    template "assembly_report.py"

}

process compile_assembly_report {

    input:
    set assembler, file(report) from ma_report.collect()
    publishDir "reports/assembly/$assembler/"

    output:
    file "${assembler}_assembly_report.csv"

    """
    echo Sample,Number of contigs,Average contig size,N50,Total assembly length,GC content,Missing data > ${assembler}_assembly_report.csv
    cat $report >> ${assembler}_assembly_report.csv
    """

}

spades_listen.ifEmpty{ exit 1, "No samples left after running Spades. Exiting." }


process process_spades {

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly) from spades_processed
    val opts from process_spades_opts
    val gsize from genome_size

    output:
    set fastq_id, file('*.assembly.fasta') into spades_assembly
    file '*.report.fasta' into spades_report

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


process assembly_mapping {

    tag { fastq_id }
    echo false

    input:
    set fastq_id, file(fastq_1), file(fastq_2), file(assembly) from assembly_mapping_input

    output:
    set fastq_id, file(assembly), 'coverages.tsv', 'sorted.bam', 'sorted.bam.bai' into mapping_coverage

    """
    bowtie2-build --threads ${task.cpus} $assembly genome_index
    bowtie2 -q --very-sensitive-local --threads ${task.cpus} -x genome_index -1 $fastq_1 -2 $fastq_2 -S mapping.sam
    samtools sort -o sorted.bam -O bam -@ ${task.cpus} mapping.sam && rm *.sam
    samtools index sorted.bam
    parallel -j ${task.cpus} samtools depth -ar {} sorted.bam \\> {}.tab  ::: \$(grep ">" $assembly | cut -c 2-)
    parallel -j ${task.cpus} echo -n {.} '"\t"' '&&' cut -f3 {} '|' paste -sd+ '|' bc >> coverages.tsv  ::: *.tab
    rm *.tab
    """
}


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


process pilon {

    tag { fastq_id }
    echo false
    publishDir 'assemblies/', mode: 'copy'

    input:
    set fastq_id, file(assembly), file(bam_file), file(bam_index) from processed_assembly_mapping

    output:
    set fastq_id, '*_polished.assembly.fasta' into pilon_processed


    """
    java -jar /NGStools/pilon-1.22.jar --genome $assembly --frags $bam_file --output ${fastq_id}_polished.assembly --changes --vcf --threads $task.cpus
    """

}


process mlst {

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly) from pilon_processed

    """
    mlst $assembly
    """

}
