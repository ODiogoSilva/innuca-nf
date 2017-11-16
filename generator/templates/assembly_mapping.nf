
process assembly_mapping {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }
    echo false

    input:
    set fastq_id, file(assembly), file(fastq) from {{ input_channel }}.join(_{{ input_channel }})

    output:
    set fastq_id, file(assembly), 'coverages.tsv', 'sorted.bam', 'sorted.bam.bai' optional true into MAIN_am_out_{{ pid }}
    set fastq_id, val("assembly_mapping"), file(".status") into STATUS_am_{{ pid }}

    when:
    params.stopAt != "assembly_mapping"

    script:
    """
    {
        bowtie2-build --threads ${task.cpus} $assembly genome_index
        bowtie2 -q --very-sensitive-local --threads ${task.cpus} -x genome_index -1 ${fastq[0]} -2 ${fastq[1]} -S mapping.sam
        samtools sort -o sorted.bam -O bam -@ ${task.cpus} mapping.sam && rm *.sam
        samtools index sorted.bam
        parallel -j ${task.cpus} samtools depth -ar {} sorted.bam \\> {}.tab  ::: \$(grep ">" $assembly | cut -c 2-)
        # Insert 0 coverage count in empty files. See Issue #2
        find . -size 0 -print0 | xargs -0 -I{} sh -c 'echo -e 0"\t"0"\t"0 > "{}"'
        parallel -j ${task.cpus} echo -n {.} '"\t"' '&&' cut -f3 {} '|' paste -sd+ '|' bc >> coverages.tsv  ::: *.tab
        rm *.tab
        if [ -f "coverages.tsv" ]
        then
            echo pass > .status
        else
            echo fail > .status
        fi
    } || {
        echo fail > .status
    }
    """
}


/** PROCESS_ASSEMBLY_MAPPING -  MAIN
Processes the results from the assembly_mapping process and filters the
assembly contigs based on coverage and length thresholds.
*/
process process_assembly_mapping {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly), file(coverage), file(bam_file), file(bam_index) from MAIN_am_out_{{ pid }}
    val min_assembly_coverage from IN_assembly_mapping_opts
    val gsize from IN_genome_size

    output:
    set fastq_id, '*_filtered.assembly.fasta', 'filtered.bam', 'filtered.bam.bai' into {{ output_channel }}
    set fastq_id, val("process_am"), file(".status") into STATUS_amp_{{ pid }}

    script:
    template "process_assembly_mapping.py"

}

{{ forks }}

