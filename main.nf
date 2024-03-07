#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2


params.irods_in = null
params.irods_in_metadata = null 
params.irods_out = null 
params.irods_out_metadata = null 

html_in = file(params.html_in)
html_template_dir = file(params.html_template_dir)
html_out = params.html_out

process IRODS {

    publishDir "output", mode: "copy"

    input: path(irods_in)
           path(irods_in_metadata)

    output: path("${irods_out}")
            path("${irods_out_metadata}")
    
    script:
    """
    irods2lustre.sh 
    """
}

process SOUP {

    publishDir "${params.outdir}/${params.sample}/${params.args.soupx_output_dir}", mode: "copy"

    input: path(input_dir)

    output: path("adjusted_counts/barcodes.tsv")
            path("adjusted_counts/genes.tsv")
            path("adjusted_counts/matrix.mtx")
            path("soup.channel.png")
    
    script:
    """
    soupx.R ${input_dir} 
    """
}

process LOAD_MTX {

    publishDir "${params.outdir}/${params.sample}", mode: "copy"

    input: path(barcodes)  
           path(genes)
           path(matrix)

    output: path(params.args.load_mtx.adata_filename) 
    
    script:
    """
    load_mtx.py \
        --adata_filename ${params.args.load_mtx.adata_filename} \
    """
}

process ADD_METADATA {

    publishDir "${params.outdir}/${params.sample}", mode: "copy"

    input: path(soup_barcodes_out)  
           file(params.metadata_file)
           file(params.sample)

    output: path("adata_metadata.h5ad")
    
    script:
    """
    add_metadata.py \
    --metadata_in ${params.metadata_file} \
    --samples ${params.sample}
    """
}


process SCRUBLET {

    publishDir "${params.outdir}/${params.sample}/adata", mode: "copy", pattern: ".h5ad"
    publishDir "${params.outdir}/${params.sample}/report/plots", mode: "copy", pattern: ".png"


    input: path(add_meta_out)
           file(params.sample)

    output: path("adata_scrublet.h5ad") 
            path("scrublet_histogram.png")
            path("scrublet_umap.png")
    
    script:
    """
    run_scrublet.py \
        --samples ${params.sample}
    """
}

process SCANPY {

    publishDir "${params.outdir}/${params.sample}/${params.args.scanpy_output_dir}", mode: "copy"

    input: path(scrub_out)
           file(params.sample)

    output: path("qcmetrics_totalcounts.png") //qcmetrics_totalcounts.png
            path("figures/violin_qcmetrics_pctcountsmt.png") //figures/violin_qcmetrics_pctcountsmt.png
            path("figures/scatter_qcmetrics.png") //figures/scatter_qcmetrics.png
            path("figures/violin_all_counts.png") //figures/violin_all_counts.png
            path("figures/violin_qc_all_counts.png") //figures/violin_qc_all_counts.png
            path("figures/highest_expr_genes_FCAImmP7241240.png") //figures/highest_expr_genes_FCAImmP7241240.png *
            path("figures/filter_genes_dispersion_scatter_FCAImmP7241240.png") //figures/filter_genes_dispersion_scatter_FCAImmP7241240.png *
            path("figures/pca_FCAImmP7241240.png") //figures/pca_FCAImmP7241240.png
            path("figures/pca_variance_ratio_FCAImmP7241240.png") //figures/pca_variance_ratio_FCAImmP7241240.png *
            path("figures/umap_neighbour_FCAImmP7241240.png") //figures/umap_neighbour_FCAImmP7241240.png *
            path("figures/umap_doubletscores.png") //figures/umap_doubletscores.png
            path("figures/umap_qc_doubletscores.png") //figures/umap_qc_doubletscores.png
            path("figures/umap_leiden.png") //figures/umap_leiden.png
            path("figures/umap_pct_counts.png") //figures/umap_pct_counts.png
            path("figures/violin_leiden_pct_counts.png") //figures/violin_leiden_pct_counts.png
            path("figures/umap_cellcycle.png") //figures/umap_cellcycle.png
            path("figures/violin_cellcycle.png") //figures/violin_cellcycle.png
            path("figures/violin_outlier.png") //figures/violin_outlier.png
            path("figures/umap_plot19.png") //figures/umap_plot19.png *
            path("figures/umap_plot20_leiden.png") //figures/umap_plot20_leiden.png *
            path("figures/violin_pctcounts_leiden.png") //figures/violin_pctcounts_leiden.png
            path("figures/violin_outlier_leiden.png") //figures/violin_outlier_leiden.png
            path("figures/umap_cellcycle_leiden.png") //figures/umap_cellcycle_leiden.png
            path("figures/pca_plot24_normalised.png") //figures/pca_plot24_normalised.png *
            path("figures/umap_metadata.png") //figures/umap_metadata.png
            path("figures/violin_n_genes_by_counts_sangerid.png") //figures/violin_n_genes_by_counts_sangerid.png *
            path("figures/violin_total_counts_sangerid.png") //figures/violin_total_counts_sangerid.png *
            path("figures/violin_pct_counts_mt_sangerid.png") //figures/violin_pct_counts_mt_sangerid.png *
            path("scanpy_adata.h5ad") //scanpy_adata.h5ad
    
    script:
    """
    scanpy_ftskin.py \
        --samples ${params.sample}
    """
}

process FILTER_READS {

    publishDir "${params.outdir}/${params.sample}/${params.args.filter_reads_output_dir}", mode: "copy"


    input: path(scrub_out)
           file(params.sample)

    output: path("adata_filter.h5ad") 
            path("qcmetrics_totalcounts.png") //qcmetrics_totalcounts.png
            path("figures/violin_qcmetrics_pctcountsmt.png") //figures/violin_qcmetrics_pctcountsmt.png
            path("figures/scatter_qcmetrics.png") //figures/scatter_qcmetrics.png
            path("figures/violin_all_counts.png") //figures/violin_all_counts.png
            path("figures/violin_qc_all_counts.png") //figures/violin_qc_all_counts.png
    
    script:
    """
    filter_reads.py \
        --samples ${params.sample}
    """
}

process FEATURE_SELECTION {

    publishDir "${params.outdir}/${params.sample}/${params.args.feature_selection_output_dir}", mode: "copy"


    input: path(filter_reads_out)

    output: path("feature_selection.h5ad") 
            path("figures/highest_expr_genes_FCAImmP7241240.png") //*
            path("figures/filter_genes_dispersion_scatter_FCAImmP7241240.png") //*
    
    script:
    """
    feature_selection.py 
    """
}

process DIMENSIONALITY {

    publishDir "${params.outdir}/${params.sample}/${params.args.dimensionality_output_dir}", mode: "copy"


    input: path(feature_selection_out)

    output: path("dimensionality.h5ad") 

    
    script:
    """
    dimensionality.py 
    """
}

process GENERATE_HTML {

    publishDir "output", mode: "copy" //"${params.outdir}/${params.sample}/${params.args.html_output_dir}"

    input: path(html_in)
           path(html_template_dir)

    output: html_out //eventually want whole directory 
    
    script:
    """
    generate_html.py \
    --html_table ${html_in} \
    --index_html ${html_out} \
    --html_template_dir ${html_template_dir}
    """
}

workflow { 
    //IRODS(irods_in,irods_in_metadata)
    SOUP(file(params.cellranger_dir))
    LOAD_MTX(SOUP.out[0],SOUP.out[1], SOUP.out[2])
    ADD_METADATA(LOAD_MTX.out[0], file(params.metadata_file), params.sample)
    SCRUBLET(ADD_METADATA.out[0], params.sample)
    SCANPY(SCRUBLET.out[0], params.sample)
    FILTER_READS(SCRUBLET.out[0], params.sample)
    FEATURE_SELECTION(FILTER_READS.out[0])
    GENERATE_HTML(html_in, html_template_dir)
}
