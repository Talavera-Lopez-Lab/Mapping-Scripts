process DOWNLOAD_DATA {
    publishDir "${params.outdir}/fastq", mode: 'copy'
    
    cpus 4
    memory { 8.GB * task.attempt }
    time { 24.hour * task.attempt }
    
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3
    
    input:
    path metadata
    
    output:
    path "fastq/*.fastq.gz", emit: fastq_files
    path "fastq/download_report.csv", emit: report
    path "fastq/download.log", emit: log
    
    script:
    """
    mkdir -p fastq
    
    python ${projectDir}/bin/download_data.py \
        --metadata ${metadata} \
        --output-dir fastq \
        --connections ${task.cpus} \
        --log-file fastq/download.log
    """
}