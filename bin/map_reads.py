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

def get_fastq_files(fastq_dir: str, metadata_row: pd.Series) -> Tuple[str, str]:
    """
    Get R1 and R2 fastq files from metadata row.
    
    Args:
        fastq_dir: Directory containing FASTQ files
        metadata_row: Row from metadata DataFrame containing file information
    
    Returns:
        Tuple of (R1 path, R2 path)
    """
    fastq_dir = Path(fastq_dir)
    
    # Get filenames from metadata
    r1_file = metadata_row.get('Comment[read1 file]')
    r2_file = metadata_row.get('Comment[read2 file]')
    
    if pd.isna(r1_file) or pd.isna(r2_file):
        raise ValueError(f"Missing read1 or read2 file information in metadata: {metadata_row.name}")
    
    r1_path = fastq_dir / r1_file
    r2_path = fastq_dir / r2_file
    
    if not r1_path.exists():
        raise FileNotFoundError(f"R1 file not found: {r1_path}")
    if not r2_path.exists():
        raise FileNotFoundError(f"R2 file not found: {r2_path}")
    
    return str(r1_path), str(r2_path)

#!/usr/bin/env python3
# (Previous imports remain the same)

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
            wl_path = os.path.join(whitelist_dir, wl_file)  # Using os.path.join instead of /
            if os.path.exists(wl_path):
                with open(wl_path) as f:
                    whitelist = set(line.strip() for line in f)
                with open(test_r1) as f:
                    barcodes = [next(f).strip()[:16] for _ in range(0, 800000, 4)]
                whitelists[wl_file] = sum(1 for bc in barcodes if bc in whitelist)
        
        # Determine best whitelist
        best_whitelist = max(whitelists.items(), key=lambda x: x[1])
        logging.info(f"Whitelist matches: {whitelists}")
        
        # Set parameters based on whitelist using os.path.join for path construction
        params = {
            "whitelist": os.path.join(whitelist_dir, best_whitelist[0]),  # Using os.path.join
            "cb_len": 16 if "3M" in best_whitelist[0] or "737K-august-2016" in best_whitelist[0] else 14,
            "umi_len": 12 if "3M" in best_whitelist[0] else 10,
            "strand": "Forward",  # Default
            "barcodeReadLength": 0 if not r1_uniform else None
        }
        
        # Verify whitelist file exists
        if not os.path.exists(params["whitelist"]):
            raise FileNotFoundError(f"Whitelist file not found: {params['whitelist']}")
        
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
    
    # Verify all input files exist
    for filepath in [fastq_r1, fastq_r2, params['whitelist'], index_dir]:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Required file not found: {filepath}")
    
    command = [
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", index_dir,
        "--readFilesIn", fastq_r2, fastq_r1,  # Note the order: R2 then R1
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
        "--soloOutFileNames", "output/", "features.tsv", "barcodes.tsv", "matrix.mtx"
    ]
    
    if params.get('barcodeReadLength') is not None:
        command.extend(["--soloBarcodeReadLength", "0"])
    
    command.extend(["--outFileNamePrefix", str(output_path / "mapping_")])
    
    # Log the exact command for debugging
    logging.info("STAR command:")
    logging.info(" ".join(command))
    
    # Run STAR
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
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
        
        # Group by sample name if needed
        if 'Extract Name' in metadata_df.columns:
            sample_groups = metadata_df.groupby('Extract Name')
        else:
            # If no Extract Name, treat each row as a separate sample
            metadata_df['_sample_id'] = [f"sample_{i}" for i in range(len(metadata_df))]
            sample_groups = metadata_df.groupby('_sample_id')
        
        for sample_name, sample_data in sample_groups:
            try:
                logging.info(f"Processing sample: {sample_name}")
                
                # Use first row for this sample
                sample_row = sample_data.iloc[0]
                
                # Get fastq files from metadata
                r1_file, r2_file = get_fastq_files(fastq_dir, sample_row)
                logging.info(f"Found FASTQ files: R1={r1_file}, R2={r2_file}")
                
                # Detect parameters
                params = detect_parameters(r1_file, r2_file, whitelist_dir)
                logging.info(f"Detected parameters for {sample_name}: {params}")
                
                # Run mapping
                sample_output_dir = os.path.join(output_dir, sample_name)
                run_star_mapping(r1_file, r2_file, sample_output_dir, params, index_dir, threads)
                
                logging.info(f"Successfully processed {sample_name}")
                
            except Exception as e:
                logging.error(f"Error processing {sample_name}: {str(e)}")
                raise
            
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        raise

if __name__ == '__main__':
    main()