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
    each file(schema) from Channel.fromPath(params.schemaPath)

    output:
    file 'chew_results'
    set fastq_id, val("chewbbaca"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}
    file '.report.json'

    when:
    params.chewbbacaRun == true

    script:
    """
    {
        echo $assembly >> input_file.txt
        chewBBACA.py AlleleCall -i input_file.txt -g $schema -o chew_results --json --cpu $task.cpus -t "Streptococcus agalactiae"
        jq -cs '.[0] * .[1]' chew_results/*/results* >> .report.json
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}