
process abricate {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { "${fastq_id} ${db}" }
    publishDir "results/annotation/abricate/${fastq_id}"

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    each db from params.abricateDatabases

    output:
    file '*.tsv'
    set fastq_id, val("abricate_${db}"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}

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