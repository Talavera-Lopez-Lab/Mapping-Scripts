#!/usr/bin/env python3
import click
import pandas as pd
import subprocess
import os
from pathlib import Path
import logging
from typing import Optional, List, Dict
from urllib.parse import quote

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

def sanitize_url(url: str) -> str:
    """Sanitize URL by encoding special characters."""
    if '://' in url:
        protocol, path = url.split('://', 1)
    else:
        protocol, path = 'ftp', url
    
    parts = path.split('/')
    encoded_parts = [quote(part) for part in parts]
    encoded_path = '/'.join(encoded_parts)
    
    return f"{protocol}://{encoded_path}"

def check_file_exists(file_path: Path) -> bool:
    """Check if file exists and has non-zero size."""
    return file_path.exists() and file_path.stat().st_size > 0

def download_file(url: str, output_path: Path, connections: int = 4) -> bool:
    """Download a file using axel with multiple connections."""
    try:
        # Check if file already exists
        if check_file_exists(output_path):
            logging.info(f"File already exists and is complete: {output_path}")
            return True

        # Sanitize URL
        sanitized_url = sanitize_url(url)
        
        # Create temporary shell script for download
        script_path = output_path.parent / f"download_{output_path.name}.sh"
        with open(script_path, 'w') as f:
            f.write(f"""#!/bin/bash
axel -n {connections} -o "{output_path}" "{sanitized_url}"
""")
        
        # Make script executable
        os.chmod(script_path, 0o755)
        
        # Run the download script
        result = subprocess.run(
            [str(script_path)],
            check=True,
            capture_output=True,
            text=True,
            encoding='latin-1'
        )
        
        # Clean up script
        script_path.unlink()
        
        # Verify file exists and has size > 0
        if check_file_exists(output_path):
            return True
        else:
            logging.error(f"Download completed but file is empty or missing: {output_path}")
            return False
            
    except subprocess.CalledProcessError as e:
        logging.error(f"Download failed for {url}: {e.stderr}")
        return False
    except Exception as e:
        logging.error(f"Unexpected error downloading {url}: {str(e)}")
        return False

def process_metadata(metadata_df: pd.DataFrame) -> Dict[str, List[str]]:
    """Process metadata to get all FASTQ URIs for each sample."""
    sample_urls = {}
    
    fastq_columns = [col for col in metadata_df.columns if 'FASTQ_URI' in col]
    if not fastq_columns:
        fastq_columns = ['Comment[FASTQ_URI]']
    
    logging.info(f"Found FASTQ columns: {fastq_columns}")
    
    for _, row in metadata_df.iterrows():
        sample_name = row.get('Extract Name', row.get('Source Name', None))
        if not sample_name:
            continue
            
        urls = []
        for col in fastq_columns:
            if pd.notna(row.get(col)):
                url = row[col].strip()
                if url:
                    urls.append(url)
        
        if urls:
            if sample_name in sample_urls:
                sample_urls[sample_name].extend(urls)
            else:
                sample_urls[sample_name] = urls
    
    return sample_urls

@click.command()
@click.option('--metadata', required=True, help='Metadata file path (SDRF format)')
@click.option('--output-dir', required=True, help='Directory to save downloaded files')
@click.option('--connections', default=4, help='Number of concurrent connections for axel')
@click.option('--log-file', default=None, help='Path to log file')
@click.option('--retry-count', default=3, help='Number of download retries')
def main(metadata: str, output_dir: str, connections: int, log_file: str, retry_count: int):
    """Download FastQ files from metadata file with multiple FASTQ URIs per sample."""
    # Convert output_dir to absolute path
    output_path = Path(output_dir).resolve()
    # Remove extra 'fastq' directory from path if it exists
    if output_path.name == 'fastq':
        output_path = output_path.parent
    
    # Create output directory
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Setup logging
    setup_logging(log_file)
    
    try:
        logging.info(f"Reading metadata from {metadata}")
        data = pd.read_csv(metadata, sep='\t')
        logging.info(f"Metadata columns: {data.columns.tolist()}")
        
        # Create download report file
        report_path = output_path / 'download_report.csv'
        download_results = []
        
        # Get sample URLs
        sample_urls = process_metadata(data)
        
        if not sample_urls:
            raise ValueError("No valid FASTQ URIs found in metadata file")
        
        # Download files
        total_files = sum(len(urls) for urls in sample_urls.values())
        current_file = 0
        
        for sample_name, urls in sample_urls.items():
            logging.info(f"Processing sample: {sample_name}")
            
            for url in urls:
                current_file += 1
                filename = url.split('/')[-1]
                file_path = output_path / filename
                
                logging.info(f"Processing {filename} ({current_file}/{total_files})")
                
                # Check if file already exists
                if check_file_exists(file_path):
                    logging.info(f"File already exists, skipping: {file_path}")
                    success = True
                    attempts = 0
                else:
                    # Try multiple times
                    success = False
                    attempts = 0
                    for attempt in range(retry_count):
                        attempts = attempt + 1
                        if attempt > 0:
                            logging.info(f"Retry attempt {attempt + 1} for {filename}")
                        success = download_file(url, file_path, connections)
                        if success:
                            break
                
                # Record result
                download_results.append({
                    'sample_name': sample_name,
                    'filename': filename,
                    'url': url,
                    'success': success,
                    'file_path': str(file_path),
                    'attempts': attempts
                })
        
        # Save download report
        pd.DataFrame(download_results).to_csv(report_path, index=False)
        
        # Summary
        successful = sum(1 for result in download_results if result['success'])
        logging.info(f"Download process completed. "
                    f"Successfully downloaded/verified {successful}/{total_files} files.")
        
        if successful == 0:
            raise RuntimeError("No files were downloaded successfully")
            
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        raise

if __name__ == '__main__':
    main()