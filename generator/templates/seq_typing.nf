

process seq_typing {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }

    input:
    set fastq_id, file(fastq_pair) from SIDE_SeqType_raw_{{ pid }}
    each refO from IN_referenceO
    each refH from IN_referenceH

    script:
    """
    seq_typing.py -f ${fastq_pair[0]} ${fastq_pair[1]} -r $refO $refH -o ./ -j 3 --extraSeq 0 --mapRefTogether
    """

}

