# ArrayExpress to STARsolo Mapping Pipeline

This repository contains a pipeline for downloading and processing data from ArrayExpress, mapping reads using STARsolo. The workflow uses Nextflow, and steps are automated to identify mapping parameters such as UMI length, barcode length, and whitelist.

## Installation

### Step 1: Install Nextflow

Download and install Nextflow to manage the pipeline:
```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p $HOME/.local/bin/
mv nextflow $HOME/.local/bin/
```

### Step 2: Set up STAR Index

To run the mapping step with STARsolo, you’ll need to provide a STAR index. We do not supply a STAR index file due to:
1. Variation in species-specific requirements.
2. Regular updates necessary with new reference genomes.
3. Custom annotation needs, such as long non-coding RNA or PolyA features.

To create a STAR index, refer to the instructions in the `preprocessing` folder. Once created, create a soft link to your STAR index file as follows:
```bash
ln -s /path/to/your/star_index reference/star_index
```

### Step 3: Install Required Tools

Ensure you have the following tools installed before proceeding:
- **Axel**
- **STARsolo**
- **Seqtk**

## Preparing Your Input Data

### Step 4: Download the SDRF File

Obtain the `sdrf.txt` file from ArrayExpress for the study you’re interested in. This file contains metadata for each sample in the study. You can locate the file on the ArrayExpress study page (see screenshot below for reference).

![Description of Image](Mapping-Scripts/array_express.png)

If you are only interested in specific samples, filter the `sdrf.txt` file to retain only those samples. You can load the `sdrf.txt` file as a DataFrame in Python to filter as needed.

### Demo Data

A sample `sdrf.txt` file with two sample entries is available in the `demo` folder for testing purposes. Note that running even this demo can take a while, so consider using `tmux` or a similar tool to manage your session.

## Running the Pipeline

Once your setup is complete, run the pipeline with:
```bash
nextflow run main.nf --metadata study.sdrf.txt --outdir results
```

This command will process the `fastq` files and automatically detect the necessary mapping parameters (e.g., UMI length, barcode length, whitelist).

## Repository Structure

- `preprocessing`: Contains example scripts for setting up STAR indexes.
- `demo`: Includes a demo `sdrf.txt` file with sample data for testing.

## Notes

- STAR indexes should be kept updated with the latest reference genome versions.
- It is recommended to run the pipeline on a server or with sufficient computing resources due to the time-intensive nature of the process.
