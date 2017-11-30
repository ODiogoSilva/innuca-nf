process chewbbaca {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    maxForks 1
    tag { fastq_id }
    scratch true
    publishDir "results/chewbbaca/${fastq_id}"
    if (params.chewbbacaQueue != null) {
        queue '${params.schemaPath}'
    }

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    each file(schema) from Channel.fromPath(params.schema_path)

    output:
    file 'chew_results'
    set fastq_id, val("chewbbaca"), file(".status") into STATUS_{{ pid }}

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