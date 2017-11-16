
process prokka {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }
    publishDir "results/annotation/prokka/${fastq_id}"

    input:
    set fastq_id, file(assembly) from {{ input_channel }}

    output:
    file "${fastq_id}/*"
    set fastq_id, val("prokka"), file(".status") into STATUS_{{ pid }}

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


