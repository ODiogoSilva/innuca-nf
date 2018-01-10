

process seq_typing {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from SIDE_SeqType_raw_{{ pid }}

    script:
    """
    echo boi
    """

}

