process DOWNLOAD_DATA {
    publishDir "/mnt/LaCIE/annaM/gut_project/raw_fastq_files/Yu_2021/results/fastq", mode: 'copy'
    
    input:
    path metadata
    
    output:
    path "*.fastq.gz", emit: fastq_files
    path "download_report.csv", emit: report
    path "download.log", emit: log
    
    script:
    """
    python ${projectDir}/bin/download_data.py \
        --metadata ${metadata} \
        --output-dir \$PWD \
        --connections ${task.cpus} \
        --log-file download.log
    """
}