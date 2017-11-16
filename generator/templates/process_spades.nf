
process process_spades {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1
    publishDir "reports/assembly/spades_filter", pattern: '*.report.csv', mode: 'copy'

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    val opts from IN_process_spades_opts
    val gsize from IN_genome_size

    output:
    set fastq_id, file('*.assembly.fasta') into {{ output_channel }}
    set fastq_id, val("process_spades"), file(".status") into STATUS_{{ pid }}
    file '*.report.csv'

    when:
    params.stopAt != "process_spades"

    script:
    template "process_spades.py"

}

{{ forks }}

