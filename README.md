# sc-analysis-nf
Pipeline to preprocess raw data from iRODS, developed using Nextflow (DSL=2)

## Installation
All setup required to utilise the pipeline. 

1. Clone the repository into local files

        git clone git@github.com:haniffalab/sc-analysis-nf.git
        cd [path/of/project/repo]

2. Create conda environment with Nextflow

        conda create -n nextflow_env -f environment.yml
        conda acitvate nextlfow_env 

## Execution
Execute the pipeline

        nextflow run main.nf -c nextflow.config
