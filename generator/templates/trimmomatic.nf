
process trimmomatic {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair), phred from {{ input_channel }}.join(SIDE_phred_{{ pid }})
    val trim_range from Channel.value("None")
    val opts from IN_trimmomatic_opts

    output:
    set fastq_id, "${fastq_id}_*P*" optional true into {{ output_channel }}
    set fastq_id, val("trimmomatic_{{ pid }}"), file(".status") into STATUS_{{ pid }}
    file '*_trimlog.txt' optional true into LOG_trimmomatic_{{ pid }}

    when:
    params.stopAt != "trimmomatic"

    script:
    template "trimmomatic.py"

}


process trimmomatic_report {

    publishDir 'reports/trimmomatic/'

    input:
    file log_files from LOG_trimmomatic_{{ pid }}.collect()

    output:
    file 'trimmomatic_report.csv'

    script:
    template "trimmomatic_report.py"

}

{{ forks }}
