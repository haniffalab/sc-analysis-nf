#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

input_exp = file(params.input)
output_exp = params.output
read_in = file(params.read_in)
read_out = params.read_out
pre_in = file(params.pre_in)
pre_out = params.pre_out
pre_vio_1 = params.vio1
pre_vio_2 = params.vio2
pre_vio_3 = params.vio3
pre_scat1 = params.scat1
pre_scat2 = params.scat2
pre_dispersion = params.dispersion
pre_adata_out = params.pre_adata_out
pca_plot1 = params.pca_plot1
pca_plot2 = params.pca_plot2
pca_plot3 = params.pca_plot3
pca_plot4 = params.pca_plot4
pca_plot5 = params.pca_plot5


process READ {

    conda 'scanpy_38'

    publishDir "results", mode: "copy"

    input: path(read_in)

    output: path("${read_out}")
    
    script:
    """
    1_read.py --input ${read_in} --output ${read_out}
    """
}

process PREPRO {

    conda 'scanpy_38'

    publishDir "results", mode: "copy"

    input: path(read_out)

    output: path("${pre_out}")
            path("${pre_vio_1}")
            path("${pre_vio_2}")
            path("${pre_vio_3}")
            path("${pre_scat1}")
            path("${pre_scat2}")
            path("${pre_dispersion}")
            path("${pre_adata_out}")
    
    script:
    """
    2_preprocess.py --input ${read_out} --output ${pre_out}
    """
}

process PCA {

    conda 'scanpy_38'

    publishDir "results", mode: "copy"

    input: path(pre_adata_out)

    output: path("${pca_plot1}")
            path("${pca_plot2}")
            path("${pca_plot3}")
            path("${pca_plot4}")
            path("${pca_plot5}")
    
    script:
    """
    3_PCA.py --input ${pre_adata_out} --output ${pca_plot1}
    """
}
workflow { 
    READ(read_in)
    PREPRO(READ.out)
    PCA(PREPRO.out[0])
}

