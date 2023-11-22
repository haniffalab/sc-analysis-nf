[![deploy-docs](https://github.com/haniffalab/sc-analysis-nf/actions/workflows/deploy-sphinx.yml/badge.svg)](https://github.com/haniffalab/sc-analysis-nf/actions/workflows/deploy-sphinx.yml)
# sc-analysis-nf Pipeline
Pipeline processing raw data to generate Anndata objects and a number of plots following ``PBMC3K ScanPy tutorial``.

[![docs](https://img.shields.io/badge/Documentation-online-purple)](https://github.com/haniffalab/sc-analysis-nf)

## Installation
All set-up required to utilise the pipeline. 

1. Clone the repository into local files.

        git clone git@github.com:haniffalab/sc-analysis-nf.git
        cd [path/of/project/repo]

2. Install Nextflow.

        # Make sure that Java v11+ is installed:
        java -version
 
        # Install Nextflow
        curl -fsSL get.nextflow.io | bash
 
        # Add Nextflow binary to your user's PATH:
        mv nextflow ~/bin/
        # OR system-wide installation:
        # sudo mv nextflow /usr/local/bin

3. Create conda environment with ScanPy
        conda create -n scanpy_env -f envs/environment.yml
        conda activate scanpy_env

## Execution
**NOTE: Ensure there are no spaces and each gene is separated by a comma**
1. Edit "marker_genes.csv" to insert genes of interest. 

        nano input/marker_genes.csv
        ## can also use VSCode or other TextEditor to edit file.

2. Execute the pipeline.

        nextflow run main.nf -c nextflow.config 
        ## View log of pipeline 
        nextflow log modest_mcclintock -t template.html > provenance.html