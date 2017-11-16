
process pilon {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { fastq_id }
    echo false
    publishDir 'results/assembly/pilon/', mode: 'copy'

    input:
    set fastq_id, file(assembly), file(bam_file), file(bam_index) from {{ input_channel }}

    output:
    set fastq_id, '*_polished.assembly.fasta' into {{ output_channel }}, pilon_report_{{ pid }}
    set fastq_id, val("pilon"), file(".status") into STATUS_{{ pid }}

    script:
    """
    {
        pilon_mem=${String.valueOf(task.memory).substring(0, String.valueOf(task.memory).length() - 1).replaceAll("\\s", "")}
        java -jar -Xms256m -Xmx\${pilon_mem} /NGStools/pilon-1.22.jar --genome $assembly --frags $bam_file --output ${fastq_id}_polished.assembly --changes --vcf --threads $task.cpus
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}

process pilon_report {

    tag { fastq_id }

    input:
    set fastq_id, file(assembly) from pilon_report_{{ pid }}

    output:
    file "*_assembly_report.csv" into pilon_report_out_{{ pid }}

    script:
    template "assembly_report.py"

}


process compile_pilon_report {

    publishDir "reports/assembly/pilon/", mode: 'copy'

    input:
    file(report) from pilon_report_out_{{ pid }}.collect()

    output:
    file "pilon_assembly_report.csv"

    """
    echo Sample,Number of contigs,Average contig size,N50,Total assembly length,GC content,Missing data > pilon_assembly_report.csv
    cat $report >> pilon_assembly_report.csv
    """
}

{{ forks }}
