#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { DOWNLOAD_DATA } from './modules/download'
include { MAP_READS } from './modules/mapping'

// Default parameters
params.metadata = null
params.outdir = 'results'
params.index_dir = "${projectDir}/reference/star_index"
params.whitelist_dir = "${projectDir}/reference/whitelists"

// Validate inputs
def validateInputs() {
    if (!params.metadata) {
        error "Must provide metadata file with --metadata"
    }
    if (!file(params.index_dir).exists()) {
        error "STAR index directory not found: ${params.index_dir}"
    }
    if (!file(params.whitelist_dir).exists()) {
        error "Whitelist directory not found: ${params.whitelist_dir}"
    }
}

// Main workflow
workflow {
    // Validate inputs
    validateInputs()
    
    // Create channel from metadata
    metadata_ch = channel.fromPath(params.metadata)
    
    // Download data
    DOWNLOAD_DATA(metadata_ch)
    
    // Create a channel combining the metadata and fastq files
    mapping_input = DOWNLOAD_DATA.out.fastq_files
        .flatten()
        .map { file -> 
            def sample_name = file.name.split('_')[0]
            return tuple(sample_name, file)
        }
        .groupTuple()
        .combine(metadata_ch)  // Combine with metadata channel
    
    // Run mapping
    MAP_READS(mapping_input)
}

// Completion handler
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Duration: $workflow.duration"
    log.info "Number of samples processed: ${workflow.stats.succeedCount}"
    if (workflow.success) {
        log.info "Pipeline completed successfully"
    } else {
        log.info "Pipeline completed with errors"
    }
1}