

process seq_typing {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }
    errorStrategy "ignore"

    input:
    set fastq_id, file(fastq_pair) from SIDE_SeqType_raw_{{ pid }}
    each refO from IN_referenceO
    each refH from IN_referenceH

    output:
    file "seq_typing.report.txt"

    script:
    """
    # Prevents read-only issues
    mkdir rematch_temp
    cp -r /NGStools/ReMatCh rematch_temp
    export PATH="\$(pwd)/rematch_temp/ReMatCh:\$PATH"

    seq_typing.py -f ${fastq_pair[0]} ${fastq_pair[1]} -r $refO $refH -o ./ -j $task.cpus --extraSeq 0 --mapRefTogether
    json_str="{'typing':{'seqtyping':'\$(cat seq_typing.report.txt)'}}"
    echo \$json_str > .report.json

    rm -r rematch_temp
    """

}

