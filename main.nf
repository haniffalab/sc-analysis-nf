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
pca_in = params.pca_in
pca_plot1 = params.pca_plot1
pca_plot2 = params.pca_plot2
pca_plot3 = params.pca_plot3
pca_plot4 = params.pca_plot4
pca_plot5 = params.pca_plot5
pca_out = params.pca_out
mgenes_in = params.mgenes_in
mgenes_plot1 = params.mgenes_plot1
mgenes_plot2 = params.mgenes_plot2
mgenes_out = params.mgenes_out
mgenes_plot3 = params.mgenes_plot3
mgenes_plot4 = params.mgenes_plot4
mgenes_plot5 = params.mgenes_plot5
mgenes_plot6 = params.mgenes_plot6
mgenes_plot7 = params.mgenes_plot7
mgenes_plot8 = params.mgenes_plot8
mgenes_plot9 = params.mgenes_plot9
// mgenes_plot10 = params.mgenes_plot10
mgenes_out2 = params.mgenes_out2
mgenes_out3 = params.mgenes_out3


process READ {

    publishDir "output", mode: "copy"
    publishDir "output", pattern: "pbmc3k.h5ad"

    input: path(read_in)

    output: path("${read_out}")
    
    script:
    """
    1_read.py --path ${read_in} --outfile ${read_out}
    """
}

process PREPRO {

    publishDir "output", mode: "copy"

    input: path(read_out)

    output: path("figures/*")
            path("${pre_adata_out}"), emit: pre_adata_out
    
    script:
    """
    2_preprocess.py --anndata_path ${read_out} --results_file ${pre_adata_out}
    """
}

process PCA {

    conda 'scanpy_38'

    publishDir "output", mode: "copy"
    publishDir "output", pattern: ".h5ad"

    input: path(pre_adata_out)

    output: path("${pca_plot1}")
            path("${pca_plot2}")
            path("${pca_plot3}")
            path("${pca_plot4}")
            path("${pca_plot5}")
            path("${pca_out}")

    
    script:
    """
    3_PCA.py --pca_adata_input ${pre_adata_out} --pca_outfile ${pca_out}
    """
}

process M_GENES {

    conda 'scanpy_38'

    publishDir "output", mode: "copy"
    publishDir "output", pattern: ".h5ad"

    input: path(pca_out)

    output: path("${mgenes_plot1}")
            path("${mgenes_plot2}")
            path("${mgenes_out}")
            path("${mgenes_plot3}")
            path("${mgenes_plot4}")
            path("${mgenes_plot5}")
            path("${mgenes_plot6}")
            path("${mgenes_plot7}")
            path("${mgenes_plot8}")
            path("${mgenes_plot9}")
           // path("${mgenes_plot10}")
            path("${mgenes_out2}")
            path("${mgenes_out3}")
            
            
    script:
    """
    4_marker_genes.py --mg_adata_input ${pca_out} --mg_outfile_mid ${mgenes_out} --mg_outfile_end ${mgenes_out2} \
    --mg_outfile_withoutx ${mgenes_out3}
    """
}

workflow { 
    READ(read_in)
    PREPRO(READ.out[0])
    PCA(PREPRO.out.pre_adata_out)
    M_GENES(PCA.out[5])
}
