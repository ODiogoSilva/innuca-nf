
/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { fastq_id }

    input:
    set fastq_id, task_name, status from {{ status_channels }}

    output:
    file 'status_*' into master_status

    """
    echo $fastq_id, $task_name, \$(cat $status) > status_${fastq_id}_${task_name}
    """
}

process compile_status {

    publishDir 'reports/status'

    input:
    file status from master_status.collect()

    output:
    file 'master_status.csv'

    """
    cat $status >> master_status.csv
    """
}

