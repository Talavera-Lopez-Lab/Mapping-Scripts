#!/usr/bin/env python3
import click
import pandas as pd
import subprocess
import os
import glob
from pathlib import Path
import logging
from typing import Optional, List, Dict, Tuple

def setup_logging(log_file: Optional[str] = None):
    """Set up logging configuration."""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        logging.basicConfig(
            level=logging.INFO,
            format=log_format,
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

def get_fastq_files(fastq_dir: str, assay_name: str) -> Tuple[str, str]:
    """
    Get R1 and R2 fastq files based on Assay Name.
    
    Args:
        fastq_dir: Directory containing FASTQ files
        assay_name: Assay Name from metadata
    
    Returns:
        Tuple of (R1 path, R2 path)
    """
    fastq_dir = Path(fastq_dir)
    
    # Look for files matching the pattern
    r1_pattern = f"{assay_name}_*R1*_001.fastq.gz"
    r2_pattern = f"{assay_name}_*R2*_001.fastq.gz"
    
    r1_files = list(fastq_dir.glob(r1_pattern))
    r2_files = list(fastq_dir.glob(r2_pattern))
    
    if not r1_files:
        raise FileNotFoundError(f"R1 file not found for pattern: {r1_pattern}")
    if not r2_files:
        raise FileNotFoundError(f"R2 file not found for pattern: {r2_pattern}")
    
    # Use the first matching file if multiple matches found
    r1_path = r1_files[0]
    r2_path = r2_files[0]
    
    return str(r1_path), str(r2_path)

def detect_parameters(r1_file: str, r2_file: str, whitelist_dir: str) -> dict:
    """Detect correct parameters for STAR mapping based on sample data."""
    logging.info("Starting parameter detection")
    output_dir = Path("test_params")
    output_dir.mkdir(exist_ok=True)
    
    try:
        # Create test files with subset of reads
        test_r1 = output_dir / "test.R1.fastq"
        test_r2 = output_dir / "test.R2.fastq"
        
        # Extract 200k reads for testing
        subprocess.run(f"zcat {r1_file} | head -800000 | seqtk sample -s100 - 200000 > {test_r1}", 
                      shell=True, check=True)
        subprocess.run(f"zcat {r2_file} | head -800000 | seqtk sample -s100 - 200000 > {test_r2}", 
                      shell=True, check=True)
        
        # Check whitelist matches
        whitelists = {
            "3M-february-2018.txt": 0,
            "737K-august-2016.txt": 0,
            "737K-arc-v1.txt": 0,
            "737K-april-2014_rc.txt": 0
        }
        
        # Read length statistics
        with open(test_r1) as f:
            r1_lengths = [len(next(f).strip()) for _ in range(0, 800000, 4)]
        r1_len = sum(r1_lengths) // len(r1_lengths)
        r1_uniform = len(set(r1_lengths)) == 1
        
        logging.info(f"R1 length: {r1_len}, uniform: {r1_uniform}")
        
        # Check whitelist matches
        for wl_file in whitelists:
            wl_path = Path(whitelist_dir) / wl_file
            if wl_path.exists():
                with open(wl_path) as f:
                    whitelist = set(line.strip() for line in f)
                with open(test_r1) as f:
                    barcodes = [next(f).strip()[:16] for _ in range(0, 800000, 4)]
                whitelists[wl_file] = sum(1 for bc in barcodes if bc in whitelist)
        
        # Determine best whitelist
        best_whitelist = max(whitelists.items(), key=lambda x: x[1])
        logging.info(f"Whitelist matches: {whitelists}")
        
        # Set parameters based on whitelist
        params = {
            "whitelist": str(Path(whitelist_dir) / best_whitelist[0]),
            "cb_len": 16 if "3M" in best_whitelist[0] or "737K-august-2016" in best_whitelist[0] else 14,
            "umi_len": 12 if "3M" in best_whitelist[0] else 10,
            "strand": "Forward",
            "barcodeReadLength": 0 if not r1_uniform else None
        }
        
        logging.info(f"Detected parameters: {params}")
        return params
        
    finally:
        # Cleanup
        if test_r1.exists():
            test_r1.unlink()
        if test_r2.exists():
            test_r2.unlink()
        if output_dir.exists():
            output_dir.rmdir()

def run_star_mapping(fastq_r1: str, fastq_r2: str, output_dir: str, params: dict,
                    index_dir: str, threads: int = 32) -> None:
    """Run STAR mapping with detected parameters."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Ensure output_path ends with a separator
    output_prefix = str(output_path) + os.sep
    
    command = [
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", index_dir,
        "--readFilesIn", fastq_r2, fastq_r1,
        "--runDirPerm", "All_RWX",
        "--soloCBwhitelist", params['whitelist'],
        "--soloType", "CB_UMI_Simple",
        "--soloCBstart", "1",
        "--soloCBlen", str(params['cb_len']),
        "--soloUMIstart", str(params['cb_len'] + 1),
        "--soloUMIlen", str(params['umi_len']),
        "--soloStrand", params['strand'],
        "--soloCBmatchWLtype", "1MM_multi_Nbase_pseudocounts",
        "--soloUMIfiltering", "MultiGeneUMI_CR",
        "--soloUMIdedup", "1MM_CR",
        "--clipAdapterType", "CellRanger4",
        "--outFilterScoreMin", "30",
        "--soloFeatures", "Gene", "GeneFull", "Velocyto",
        "--readFilesCommand", "zcat",
        "--soloOutFileNames", "output/", "features.tsv", "barcodes.tsv", "matrix.mtx",
        "--outFileNamePrefix", output_prefix
    ]
    
    if params.get('barcodeReadLength') is not None:
        command.extend(["--soloBarcodeReadLength", "0"])
    
    logging.info(f"Running STAR command in directory: {output_path}")
    logging.info(" ".join(command))
    
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        
        # Verify output directory exists
        output_dir_path = output_path / "output"
        if not output_dir_path.exists():
            raise RuntimeError(f"STAR completed but output directory not found: {output_dir_path}")
            
    except subprocess.CalledProcessError as e:
        logging.error(f"STAR mapping failed with error:")
        logging.error(f"stdout: {e.stdout}")
        logging.error(f"stderr: {e.stderr}")
        raise

@click.command()
@click.option('--fastq-dir', required=True, help='Directory containing FASTQ files')
@click.option('--output-dir', required=True, help='Output directory for mapped files')
@click.option('--index-dir', required=True, help='STAR index directory')
@click.option('--whitelist-dir', required=True, help='Directory containing whitelist files')
@click.option('--metadata', required=True, help='Metadata file path')
@click.option('--threads', default=32, help='Number of threads for STAR')
@click.option('--log-file', default=None, help='Log file path')
def main(fastq_dir: str, output_dir: str, index_dir: str, whitelist_dir: str,
         metadata: str, threads: int, log_file: str):
    """Process multiple samples with STAR."""
    setup_logging(log_file)
    
    try:
        # Read metadata
        metadata_df = pd.read_csv(metadata, sep="\t")
        logging.info(f"Metadata columns: {metadata_df.columns.tolist()}")
        
        if 'Assay Name' not in metadata_df.columns:
            raise ValueError("Metadata must contain 'Assay Name' column")
        
        # Process each unique Assay Name
        for assay_name in metadata_df['Assay Name'].unique():
            try:
                logging.info(f"Processing sample: {assay_name}")
                
                # Create output directory structure
                sample_output_dir = Path(output_dir) / "StarMapped" / assay_name
                
                # Get fastq files using Assay Name pattern
                r1_file, r2_file = get_fastq_files(fastq_dir, assay_name)
                logging.info(f"Found FASTQ files: R1={r1_file}, R2={r2_file}")
                
                # Detect parameters
                params = detect_parameters(r1_file, r2_file, whitelist_dir)
                logging.info(f"Detected parameters for {assay_name}: {params}")
                
                # Run mapping
                run_star_mapping(r1_file, r2_file, sample_output_dir, params, index_dir, threads)
                
                logging.info(f"Successfully processed {assay_name}")
                
            except Exception as e:
                logging.error(f"Error processing {assay_name}: {str(e)}")
                raise
            
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        raise

if __name__ == '__main__':
    main()