process MAP_READS {
    publishDir "${params.outdir}/mapped", mode: 'copy'
    
    cpus 32
    memory { 64.GB * task.attempt }
    time { 24.hour * task.attempt }
    
    errorStrategy { task.attempt <= 3 ? 'retry' : 'finish' }
    maxRetries 3
    
    input:
    tuple val(sample_id), path(fastq_files), path(metadata)
    
    output:
    tuple val(sample_id), path("${sample_id}/output"), emit: mapped_data
    path "${sample_id}/mapping.log", emit: logs
    
    script:
    """
    python ${projectDir}/bin/map_reads.py \
        --fastq-dir . \
        --output-dir ${sample_id} \
        --index-dir ${params.index_dir} \
        --whitelist-dir ${params.whitelist_dir} \
        --metadata ${metadata} \
        --threads ${task.cpus} \
        --log-file ${sample_id}/mapping.log
    """
}