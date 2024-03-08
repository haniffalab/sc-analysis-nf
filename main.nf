#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2


params.irods_in = null
params.irods_in_metadata = null 
params.irods_out = null 
params.irods_out_metadata = null 

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

    publishDir "${params.outdir}/${params.sample}/report/plots/${params.args.soupx_output_dir}", mode: "copy", pattern: "*.png"
    publishDir "${params.outdir}/${params.sample}/${params.adata_output_dir}", mode: "copy", pattern: "*.tsv"
    publishDir "${params.outdir}/${params.sample}/${params.adata_output_dir}", mode: "copy", pattern: "*.mtx"

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

    publishDir "${params.outdir}/${params.sample}/${params.adata_output_dir}", mode: "copy", pattern: "*.h5ad"

    input: path(barcodes)  
           path(genes)
           path(matrix)

    output: path("${params.mtx_adata}")
    
    script:
    """
    load_mtx.py \
        --mtx_out ${params.mtx_adata}
    """
}

process ADD_METADATA {

    publishDir "${params.outdir}/${params.sample}/${params.adata_output_dir}", mode: "copy", pattern: "*.h5ad"

    input: path(soup_barcodes_out)  
           file(params.metadata_file)
           file(params.sample)

    output: path("${params.meta_adata}")
    
    script:
    """
    add_metadata.py \
        --add_meta_in ${soup_barcodes_out} \
        --metadata_in ${params.metadata_file} \
        --samples ${params.sample} \
        --add_metadata_out ${params.meta_adata}
    """
}


process SCRUBLET {

    publishDir "${params.outdir}/${params.sample}/report/plots/${params.args.scrublet_output_dir}", mode: "copy", pattern: "*.png"
    publishDir "${params.outdir}/${params.sample}/${params.adata_output_dir}", mode: "copy", pattern: "*.h5ad"


    input: path(add_meta_out)
           file(params.sample)

    output: path("${params.scrublet_adata}")
            path("scrublet_histogram.png")
            path("scrublet_umap.png")
    
    script:
    """
    scrublet_sc.py \
        --samples ${params.sample} \
        --scrub_in ${add_meta_out} \
        --scrub_out ${params.scrublet_adata}
    """
}

process FILTER_READS {

    publishDir "${params.outdir}/${params.sample}/report/plots/${params.args.filter_reads_output_dir}", mode: "copy", pattern: "figures/*.png"
    publishDir "${params.outdir}/${params.sample}/report/plots/${params.args.filter_reads_output_dir}", mode: "copy", pattern: "*.png"
    publishDir "${params.outdir}/${params.sample}/${params.adata_output_dir}", mode: "copy", pattern: "*.h5ad"


    input: path(scrub_out)
           file(params.sample)

    output: path("${params.filter_adata}")
            path("qcmetrics_totalcounts.png") 
            path("figures/violin_qcmetrics_pctcountsmt.png") 
            path("figures/scatter_qcmetrics.png")
            path("figures/violin_all_counts.png") 
            path("figures/violin_qc_all_counts.png")
    
    script:
    """
    filter_reads.py \
        --samples ${params.sample} \
        --filter_in ${scrub_out} \
        --filter_out ${params.filter_adata}
    """
}

process FEATURE_SELECTION {

    publishDir "${params.outdir}/${params.sample}/report/plots/${params.args.feature_selection_output_dir}", mode: "copy", pattern: "figures/*.png"
    publishDir "${params.outdir}/${params.sample}/${params.adata_output_dir}", mode: "copy", pattern: "*.h5ad"

    input: path(filter_reads_out)

    output: path("${params.feature_selection_adata}")
            path("figures/highest_expr_genes_FCAImmP7241240.png") //*
            path("figures/filter_genes_dispersion_scatter_FCAImmP7241240.png") //*
    
    script:
    """
    feature_selection.py \
        --feature_in ${filter_reads_out} \
        --feature_out ${params.feature_selection_adata}
    """
}

process DIMENSIONALITY {

    publishDir "${params.outdir}/${params.sample}/report/plots/${params.args.dimensionality_output_dir}", mode: "copy", pattern: "figures/*.png"
    publishDir "${params.outdir}/${params.sample}/${params.adata_output_dir}", mode: "copy", pattern: "*.h5ad"

    input: path(feature_selection_out)

    output: path("${params.dimensionality_adata}")
            path("figures/pca_FCAImmP7241240.png") 
            path("figures/pca_variance_ratio_FCAImmP7241240.png") // *
            path("figures/umap_neighbour_FCAImmP7241240.png") // *
            path("figures/umap_doubletscores.png") 
            path("figures/umap_qc_doubletscores.png") 
            path("figures/umap_leiden.png")
            path("figures/umap_pct_counts.png")
            path("figures/violin_leiden_pct_counts.png") 
            path("figures/umap_cellcycle.png") 
            path("figures/violin_cellcycle.png") 
            path("figures/violin_outlier.png") 
            path("figures/umap_plot19.png") // *
            path("figures/umap_plot20_leiden.png") // *
            path("figures/violin_pctcounts_leiden.png") 
            path("figures/violin_outlier_leiden.png") 
            path("figures/umap_cellcycle_leiden.png")


    
    script:
    """
    dimensionality_reduction.py \
        --dimensionality_in ${feature_selection_out} \
        --dimensionality_out ${params.dimensionality_adata}
    """
}

process COUNT_TRANSFORMATION {

    publishDir "${params.outdir}/${params.sample}/report/plots/${params.args.count_transformation_output_dir}", mode: "copy", pattern: "figures/*.png"
    publishDir "${params.outdir}/${params.sample}/${params.adata_output_dir}", mode: "copy", pattern: "*.h5ad"

    input: path(counts_in)

    output: path("${params.count_transformation}")
            path("figures/pca_plot24_normalised.png") // *
            path("figures/umap_metadata.png") //
            path("figures/violin_n_genes_by_counts_sangerid.png") // *
            path("figures/violin_total_counts_sangerid.png") // *
            path("figures/violin_pct_counts_mt_sangerid.png") // *
    
    script:
    """
    count_transformation.py \
        --counts_in ${counts_in} \
        --counts_out ${params.count_transformation}
    """
}

process GENERATE_HTML {

    publishDir "${params.outdir}/${params.sample}", mode: "copy"

    input: path(html_dist_dir)

    output: path("report") //eventually want whole directory 
    
    script:
    """
    generate_html.py \
    --html_table ${html_dist_dir}
    """
}

workflow { 
    html_dist_dir = projectDir/'html'/'dist'
    // IRODS(irods_in,irods_in_metadata)
    SOUP(file(params.cellranger_dir))
    LOAD_MTX(SOUP.out[0],SOUP.out[1], SOUP.out[2])
    ADD_METADATA(LOAD_MTX.out[0], file(params.metadata_file), params.sample)
    SCRUBLET(ADD_METADATA.out[0], params.sample)
    FILTER_READS(SCRUBLET.out[0], params.sample)
    FEATURE_SELECTION(FILTER_READS.out[0])
    DIMENSIONALITY(FEATURE_SELECTION.out[0])
    COUNT_TRANSFORMATION(DIMENSIONALITY.out[0])
    GENERATE_HTML(file(html_dist_dir))
}