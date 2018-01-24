
process mlst {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly) from {{ input_channel }}
    val expectedSpecies from Channel.value(params.mlstSpecies)

    output:
    file '*.mlst.txt' into MAIN_mlst_out_{{ pid }}
    set fastq_id, val("mlst"), file(".status"), file(".warning"), file(".fail") into STATUS_{{ pid }}

    when:
    params.mlstRun  == true && params.annotationRun

    script:
    """
    {
        mlst $assembly >> ${fastq_id}.mlst.txt
        json_str="{'expectedSpecies':\'$expectedSpecies\', 'species': '\$(cat *.mlst.txt | cut -f2)'}"
        echo \$json_str > .report.json
        echo pass > .status
    } || {
        echo fail > .status
    }
    """
}

process compile_mlst {

    publishDir "results/annotation/mlst/"

    input:
    file res from MAIN_mlst_out_{{ pid }}.collect()

    output:
    file "mlst_report.tsv"

    when:
    params.mlstRun == true && params.annotationRun

    script:
    """
    cat $res >> mlst_report.tsv
    """
}


