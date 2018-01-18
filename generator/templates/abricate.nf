
process abricate {

    // Send POST request to platform
    {% with id_modifier=".1" %}
    {% include "post.txt" ignore missing %}
    {% endwith %}

    tag { "${fastq_id} ${db}" + " getStats"}
    publishDir "results/annotation/abricate/${fastq_id}"

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    each db from params.abricateDatabases

    output:
    file '*.tsv' into abricate_out_{{ pid }}
    set fastq_id, val("abricate_${db}"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}

    when:
    params.abricateRun == true && params.annotationRun

    script:
    """
    {
        # Run abricate
        abricate --db $db $assembly > ${fastq_id}_abr_${db}.tsv
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}


process process_abricate {

    // Send POST request to platform
    {% with id_modifier=".2" %}
    {% include "post.txt" ignore missing %}
    {% endwith %}

    input:
    file abricate_file from abricate_out_{{ pid }}.collect()

    script:
    template "process_abricate.py"


}



