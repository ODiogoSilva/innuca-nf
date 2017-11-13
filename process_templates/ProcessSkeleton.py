class Skeletons:

    def __getitem__(self, item):
        return getattr(self, item)

    integrity_coverage = '''
process integrity_coverage {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(fastq_pair) from IN_fastq_raw
    val gsize from IN_genome_size
    val cov from IN_min_coverage
    // This channel is for the custom options of the integrity_coverage.py
    // script. See the script's documentation for more information.
    val opts from Channel.value('')

    output:
    set fastq_id,
        file(fastq_pair),
        file('*_encoding'),
        file('*_phred'),
        file('*_coverage'),
        file('*_max_len') into MAIN_integrity
    file('*_report') into LOG_report_coverage1

    script:
    template "integrity_coverage.py"

}}

// TRIAGE OF CORRUPTED SAMPLES
LOG_corrupted = Channel.create()
MAIN_PreCoverageCheck = Channel.create()
// Corrupted samples have the 2nd value with 'corrupt'
MAIN_integrity.choice(LOG_corrupted, MAIN_PreCoverageCheck) {{
    a -> a[2].text == "corrupt" ? 0 : 1
}}

// TRIAGE OF LOW COVERAGE SAMPLES
{output_channel} = Channel.create()
SIDE_phred_{pid} = Channel.create()
SIDE_max_len_{pid} = Channel.create()

MAIN_PreCoverageCheck
// Low coverage samples have the 4th value of the Channel with 'fail'
    .filter{{ it[4].text != "fail" }}
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'SIDE_phred'
    .separate({output_channel}, SIDE_phred_{pid}, SIDE_max_len_{pid}){{
        a -> [ [a[0], a[1]], [a[0], a[3].text], [a[0], a[5].text] ]
    }}

/** REPORT_COVERAGE - PLUG-IN
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_initial.csv'
*/
process report_coverage {{

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage/'

    input:
    file(report) from LOG_report_coverage1.filter{{ it.text != "corrupt" }}.collect()

    output:
    file 'estimated_coverage_initial.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_initial.csv
    cat $report >> estimated_coverage_initial.csv
    """
}}

/** REPORT_CORRUPT - PLUG-IN
This process will report the corrupted samples and write the results to
'reports/corrupted/corrupted_samples.txt'
*/
process report_corrupt {{

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/corrupted/'

    input:
    val fastq_id from LOG_corrupted.collect{{it[0]}}

    output:
    file 'corrupted_samples.txt'

    """
    echo ${{fastq_id.join(",")}} | tr "," "\\n" >> corrupted_samples.txt
    """

}}

'''

    check_coverage = '''
process integrity_coverage_2 {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}
    cpus 1

    input:
    set fastq_id, file(fastq_pair) from {input_channel}
    val gsize from IN_genome_size
    val cov from IN_min_coverage
    // Use -e option for skipping encoding guess
    val opts from Channel.value('-e')

    output:
    set fastq_id,
        file(fastq_pair),
        file('*_coverage'),
        file('*_max_len') into MAIN_integrity_{pid}
    file('*_report') into LOG_report_coverage_{pid}

    script:
    template "integrity_coverage.py"
}}

// Low coverage samples have the 2nd value of the Channel with 'fail'
{output_channel} = Channel.create()
SIDE_max_len_{pid} = Channel.create()

MAIN_integrity_{pid}
// Low coverage samples have the 2nd value of the Channel with 'fail'
    .filter{{ it[2].text != "fail" }}
    .separate({output_channel}, SIDE_max_len_{pid}){{
        a -> [ [a[0], a[1]], [a[0], a[3].text]]
    }}


/** REPORT_COVERAGE_2 - PLUG-IN
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_second.csv'
*/
process report_coverage_2 {{

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage/'

    input:
    file(report) from LOG_report_coverage_{pid}.filter{{ it.text != "corrupt" }}.collect()

    output:
    file 'estimated_coverage_second.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_second.csv
    cat $report >> estimated_coverage_second.csv
    """
}}

'''

    fastqc_trimmomatic = '''
process fastqc {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}

    input:
    set fastq_id, file(fastq_pair) from {input_channel}
    val ad from Channel.value('None')

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into MAIN_fastqc_out_{pid}
    set fastq_id, val("fastqc"), file(".status") into STATUS_{pid}

    when:
    params.stopAt != "fastqc"

    script:
    template "fastqc.py"
}}


process fastqc_report {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid}"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid}"
    }}

    tag {{ fastq_id }}
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc/run_1/', pattern: '*summary.txt', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_{pid}
    val opts from Channel.value("--ignore-tests")

    output:
    set fastq_id, file(fastq_pair), '.status', 'optimal_trim' into MAIN_fastqc_trim_{pid}
    file '*_trim_report' into LOG_trim_{pid}
    file "*_status_report" into LOG_fastqc_report_{pid}
    file "${{fastq_id}}_*_summary.txt" optional true

    script:
    template "fastqc_report.py"

}}

process compile_fastqc_status {{

    publishDir 'reports/fastqc/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report_{pid}.collect()

    output:
    file 'FastQC_1run_report.csv'

    """
    echo Sample, Failed? >> FastQC_1run_report.csv
    cat $rep >> FastQC_1run_report.csv
    """

}}


process trimmomatic {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId 4"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 4"
    }}

    tag {{ fastq_id }}

    input:
    set fastq_id, file(fastq_pair), trim_range, phred from MAIN_fastqc_trim_{pid}.phase(SIDE_phred_{pid}).map{{ [it[0][0], it[0][1], file(it[0][3]).text, it[1][1]] }}
    val opts from IN_trimmomatic_opts

    output:
    set fastq_id, "${{fastq_id}}_*P*" optional true into {output_channel}
    set fastq_id, val("trimmomatic"), file(".status") into STATUS_trimmomatic_{pid}
    file '*_trimlog.txt' optional true into LOG_trimmomatic_{pid}

    when:
    params.stopAt != "trimmomatic"

    script:
    template "trimmomatic.py"

}}


process trimmomatic_report {{

    publishDir 'reports/trimmomatic/'

    input:
    file log_files from LOG_trimmomatic_{pid}.collect()

    output:
    file 'trimmomatic_report.csv'

    script:
    template "trimmomatic_report.py"

}}

'''

    trimmomatic = """
process trimmomatic {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}

    input:
    set fastq_id, file(fastq_pair), phred from {input_channel}.phase(SIDE_phred_{pid}).map{{ [it[0][0], it[0][1], it[1][1]] }}
    val trim_range from Channel.value("None")
    val opts from IN_trimmomatic_opts

    output:
    set fastq_id, "${{fastq_id}}_*P*" optional true into {output_channel}
    set fastq_id, val("trimmomatic"), file(".status") into STATUS_trimmomatic_{pid}
    file '*_trimlog.txt' optional true into LOG_trimmomatic_{pid}

    when:
    params.stopAt != "trimmomatic"

    script:
    template "trimmomatic.py"

}}   

process trimmomatic_report {{

    publishDir 'reports/trimmomatic/'

    input:
    file log_files from LOG_trimmomatic_{pid}.collect()

    output:
    file 'trimmomatic_report.csv'

    script:
    template "trimmomatic_report.py"

}}

"""

    fastqc = '''
process fastqc2 {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}

    input:
    set fastq_id, file(fastq_pair) from {input_channel}
    val ad from IN_adapters

    output:
    set fastq_id, file(fastq_pair), file('pair_1*'), file('pair_2*') optional true into MAIN_fastqc_out_{pid}
    set fastq_id, val("fastqc2"), file(".status") into STATUS_fastqc_{pid}

    when:
    params.stopAt != "fastqc2"

    script:
    template "fastqc.py"
}}


process fastqc2_report {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid}"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid}"
    }}

    tag {{ fastq_id }}
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc/run_2/', pattern: '*summary.txt', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_{pid}
    val opts from Channel.value("")

    output:
    set fastq_id, file(fastq_pair), '.status' into MAIN_fastqc_report_{pid}
    file "*_status_report" into LOG_fastqc_report_{pid}
    file "${{fastq_id}}_*_summary.txt" optional true

    script:
    template "fastqc_report.py"

}}


process compile_fastqc_status2 {{

    publishDir 'reports/fastqc/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report_{pid}.collect()

    output:
    file 'FastQC_2run_report_{pid}.csv'

    """
    echo Sample, Failed? >> FastQC_2run_report_{pid}.csv
    cat $rep >> FastQC_2run_report_{pid}.csv
    """

}}

{output_channel} = Channel.create()

MAIN_fastqc_report_{pid}
        .filter{{ it[2].text == "pass" }}
        .map{{ [it[0], it[1]] }}
        .into({output_channel})

'''

    spades = """
process spades {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}
    publishDir 'results/assembly/spades/', pattern: '*_spades.assembly.fasta', mode: 'copy'

    input:
    set fastq_id, file(fastq_pair), max_len from {input_channel}.phase(SIDE_max_len_{pid}).map{{ [it[0][0], it[0][1], it[1][1]] }}
    val opts from IN_spades_opts
    val kmers from IN_spades_kmers

    output:
    set fastq_id, file('*_spades.assembly.fasta') optional true into {output_channel}
    set fastq_id, val("spades"), file(".status") into STATUS_spades_{pid}

    when:
    params.stopAt != "spades"

    script:
    template "spades.py"

}}

"""

    process_spades = """
process process_spades {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}
    // This process can only use a single CPU
    cpus 1
    publishDir "reports/assembly/spades_filter", pattern: '*.report.csv', mode: 'copy'

    input:
    set fastq_id, file(assembly) from {input_channel}
    val opts from IN_process_spades_opts
    val gsize from IN_genome_size

    output:
    set fastq_id, file('*.assembly.fasta') into {output_channel}
    set fastq_id, val("process_spades"), file(".status") into STATUS_process_spades_{pid}
    file '*.report.csv'

    when:
    params.stopAt != "process_spades"

    script:
    template "process_spades.py"

}}

"""

    assembly_mapping = '''
process assembly_mapping {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}
    echo false

    input:
    set fastq_id, file(fastq_1), file(fastq_2), file(assembly) from {input_channel}.phase(_{input_channel}).map{{ [ it[0][0], it[1][1][0], it[1][1][1], it[0][1] ] }}

    output:
    set fastq_id, file(assembly), 'coverages.tsv', 'sorted.bam', 'sorted.bam.bai' optional true into MAIN_am_out_{pid}
    set fastq_id, val("assembly_mapping"), file(".status") into STATUS_am_{pid}

    when:
    params.stopAt != "assembly_mapping"

    script:
    """
    {{
        bowtie2-build --threads ${{task.cpus}} $assembly genome_index
        bowtie2 -q --very-sensitive-local --threads ${{task.cpus}} -x genome_index -1 $fastq_1 -2 $fastq_2 -S mapping.sam
        samtools sort -o sorted.bam -O bam -@ ${{task.cpus}} mapping.sam && rm *.sam
        samtools index sorted.bam
        parallel -j ${{task.cpus}} samtools depth -ar {{}} sorted.bam \\\\> {{}}.tab  ::: \$(grep ">" $assembly | cut -c 2-)
        # Insert 0 coverage count in empty files. See Issue #2
        find . -size 0 -print0 | xargs -0 -I{{}} sh -c 'echo -e 0"\\t"0"\\t"0 > "{{}}"'
        parallel -j ${{task.cpus}} echo -n {{.}} '"\\t"' '&&' cut -f3 {{}} '|' paste -sd+ '|' bc >> coverages.tsv  ::: *.tab
        rm *.tab
        if [ -f "coverages.tsv" ]
        then
            echo pass > .status
        else
            echo fail > .status
        fi
    }} || {{
        echo fail > .status
    }}
    """
}}

process process_assembly_mapping {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid}"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid}"
    }}

    tag {{ fastq_id }}
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly), file(coverage), file(bam_file), file(bam_index) from MAIN_am_out_{pid}
    val min_assembly_coverage from IN_assembly_mapping_opts
    val gsize from IN_genome_size

    output:
    set fastq_id, '*_filtered.assembly.fasta', 'filtered.bam', 'filtered.bam.bai' into {output_channel}
    set fastq_id, val("process_am"), file(".status") into STATUS_process_am_{pid}

    script:
    template "process_assembly_mapping.py"

}}


'''

    pilon = '''
process pilon {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}
    
    tag {{ fastq_id }}
    echo false
    publishDir 'results/assembly/pilon/', mode: 'copy'

    input:
    set fastq_id, file(assembly), file(bam_file), file(bam_index) from {input_channel}

    output:
    set fastq_id, '*_polished.assembly.fasta' into {output_channel}, pilon_report_{pid}
    set fastq_id, val("pilon"), file(".status") into STATUS_pilon_{pid}

    script:
    """
    {{
        pilon_mem=${{String.valueOf(task.memory).substring(0, String.valueOf(task.memory).length() - 1).replaceAll("\\\\s", "")}}
        java -jar -Xms256m -Xmx\\${{pilon_mem}} /NGStools/pilon-1.22.jar --genome $assembly --frags $bam_file --output ${{fastq_id}}_polished.assembly --changes --vcf --threads $task.cpus
        echo pass > .status
    }} || {{
        echo fail > .status
    }}
    """

}}

process pilon_report {{

    tag {{ fastq_id }}

    input:
    set fastq_id, file(assembly) from pilon_report_{pid}

    output:
    file "*_assembly_report.csv" into pilon_report_out_{pid}

    script:
    template "assembly_report.py"

}}


process compile_pilon_report {{

    publishDir "reports/assembly/pilon/", mode: 'copy'

    input:
    file(report) from pilon_report_out_{pid}.collect()

    output:
    file "pilon_assembly_report.csv"

    """
    echo Sample,Number of contigs,Average contig size,N50,Total assembly length,GC content,Missing data > pilon_assembly_report.csv
    cat $report >> pilon_assembly_report.csv
    """
}}

'''

    mlst = '''
process mlst {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}
    // This process can only use a single CPU
    cpus 1

    input:
    set fastq_id, file(assembly) from {input_channel}

    output:
    file '*.mlst.txt' into MAIN_mlst_out_{pid}
    set fastq_id, val("mlst"), file(".status") into STATUS_mlst_{pid}

    when:
    params.mlstRun  == true && params.annotationRun

    script:
    """
    {{
        mlst $assembly >> ${{fastq_id}}.mlst.txt
        echo pass > .status
    }} || {{
        echo fail > .status
    }}
    """
}}

process compile_mlst {{

    publishDir "results/annotation/mlst/"

    input:
    file res from MAIN_mlst_out_{pid}.collect()

    output:
    file "mlst_report.tsv"

    when:
    params.mlstRun == true && params.annotationRun

    script:
    """
    cat $res >> mlst_report.tsv
    """
}}

'''

    abricate = '''
process abricate {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ "${{fastq_id}} ${{db}}" }}
    publishDir "results/annotation/abricate/${{fastq_id}}"

    input:
    set fastq_id, file(assembly) from {input_channel}
    each db from params.abricateDatabases

    output:
    file '*.tsv'
    set fastq_id, val("abricate_${{db}}"), file(".status") into STATUS_abricate_{pid}

    when:
    params.abricateRun == true && params.annotationRun

    script:
    """
    {{
        abricate --db $db $assembly > ${{fastq_id}}_abr_${{db}}.tsv
        echo pass > .status
    }} || {{
        echo fail > .status
    }}
    """

}}

'''

    prokka = '''
process prokka {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    tag {{ fastq_id }}
    publishDir "results/annotation/prokka/${{fastq_id}}"

    input:
    set fastq_id, file(assembly) from {input_channel}

    output:
    file "${{fastq_id}}/*"
    set fastq_id, val("prokka"), file(".status") into STATUS_prokka_{pid}

    when:
    params.prokkaRun == true && params.annotationRun

    script:
    """
    {{
        prokka --outdir $fastq_id --cpus $task.cpus --centre UMMI --compliant \
               --increment 10 $assembly
        echo pass > .status
    }} || {{
        echo fail > .status
    }}
    """

}}

'''

    chewbbaca = '''
process chewbbaca {{

    // Send POST request to platform
    if ( params.platformHTTP != null ) {{
        beforeScript "startup_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId {pid} $params.platformHTTP"
    }}

    maxForks 1
    tag {{ fastq_id }}
    echo false
    scratch true
    publishDir "results/chewbbaca/${{fastq_id}}"

    input:
    set fastq_id, file(assembly) from {input_channel}
    each file(schema) from Channel.fromPath(params.schema_path)

    output:
    file 'chew_results'
    set fastq_id, val("chewbbaca"), file(".status") into {status_channel}

    when:
    params.chewbbacaRun == true

    script:
    """
    {{
        echo $assembly >> input_file.txt
        chewBBACA.py AlleleCall -i input_file.txt -g $schema -o chew_results --json --cpu $task.cpus -t "Streptococcus agalactiae"
        echo pass > .status
    }} || {{
        echo fail > .status
    }}
    """

}}

'''