
process patho_typing {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id + " getStats" }

    input:
    set fastq_id, file(fastq_pair) from SIDE_PathoType_raw_{{ pid }}
    val species from IN_pathoSpecies

    script:
    """
    # Prevents read-only issues
    mkdir rematch_temp
    cp -r /NGStools/ReMatCh rematch_temp
    export PATH="\$(pwd)/rematch_temp/ReMatCh:\$PATH"

    patho_typing.py -f ${fastq_pair[0]} ${fastq_pair[1]} -o ./ -j $task.cpus --trueCoverage --species $species
    """

}

