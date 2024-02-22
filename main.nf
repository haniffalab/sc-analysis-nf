#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2


//irods_in = params.irods_in
params.irods_in = null
//irods_in_metadata = params.irods_in_metadata
params.irods_in_metadata = null 
//irods_out = params.irods_out
params.irods_out = null 
//irods_out_metadata = params.irods_out_metadata
params.irods_out_metadata = null 
soup_in = params.soup_in 
soup_raw_in = params.soup_raw_in 
soup_barcodes_out = params.soup_barcodes_out
soup_genes_out = params.soup_genes_out
soup_matrix_out = params.soup_matrix_out
soup_plot = params.soup_plot
load_in = params.load_in
load_raw_out = params.load_raw_out
load_metadata = params.load_metadata
load_meta_out = params.load_meta_out
scrub_in = params.scrub_in
scrub_out= params.scrub_out
scrub_histogram = params.scrub_histogram
scrub_umap = params.scrub_umap
scanpy_in = params.scanpy_in
scanpy_plot1 = params.scanpy_plot1
scanpy_plot2 = params.scanpy_plot2
scanpy_plot3 = params.scanpy_plot3
scanpy_plot4 = params.scanpy_plot4
scanpy_plot5 = params.scanpy_plot5
scanpy_plot6 = params.scanpy_plot6
scanpy_plot7 = params.scanpy_plot7
scanpy_plot8 = params.scanpy_plot8
scanpy_plot9 = params.scanpy_plot9
scanpy_plot10 = params.scanpy_plot10
scanpy_plot11 = params.scanpy_plot11
scanpy_plot12 = params.scanpy_plot12
scanpy_plot13 = params.scanpy_plot13
scanpy_plot14 = params.scanpy_plot14
scanpy_plot15 = params.scanpy_plot15
scanpy_plot16 = params.scanpy_plot16
scanpy_plot17 = params.scanpy_plot17
scanpy_plot18 = params.scanpy_plot18
scanpy_plot19 = params.scanpy_plot19
scanpy_plot20 = params.scanpy_plot20
scanpy_plot21 = params.scanpy_plot21
scanpy_plot22 = params.scanpy_plot22
scanpy_plot23 = params.scanpy_plot23
scanpy_plot24 = params.scanpy_plot24
scanpy_plot25 = params.scanpy_plot25
scanpy_plot26 = params.scanpy_plot26
scanpy_plot27 = params.scanpy_plot27
scanpy_plot28 = params.scanpy_plot28
scanpy_out = params.scanpy_out
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

process LOAD_DATA {

    publishDir "${params.outdir}/${params.sample}", mode: "copy"

    input: path(soup_barcodes_out)
           path(soup_genes_out)
           path(soup_matrix_out)
           path(load_metadata)

    output: path("${load_raw_out}")
            path("${load_meta_out}")
    
    script:
    """
    loaddata.py --load_in ${soup_matrix_out} --metadata_in ${load_metadata} --load_raw_out ${load_raw_out} --load_meta_out ${load_meta_out}
    """
}

process SCRUBLET {

    publishDir "output/FCAImmP7241240", mode: "copy"

    input: path(load_meta_out)

    output: path("${scrub_out}")
            path("${scrub_histogram}")
            path("${scrub_umap}")
    
    script:
    """
    scrublet_ftskin.py 
    """
}

process SCANPY {

    publishDir "output/FCAImmP7241240", mode: "copy"

    input: path(scrub_out)

    output: path("${scanpy_plot1}")
            path("${scanpy_plot2}")
            path("${scanpy_plot3}")
            path("${scanpy_plot4}")
            path("${scanpy_plot5}")
            path("${scanpy_plot6}")
            path("${scanpy_plot7}")
            path("${scanpy_plot8}")
            path("${scanpy_plot9}")
            path("${scanpy_plot10}")
            path("${scanpy_plot11}")
            path("${scanpy_plot12}")
            path("${scanpy_plot13}")
            path("${scanpy_plot14}")
            path("${scanpy_plot15}")
            path("${scanpy_plot16}")
            path("${scanpy_plot17}")
            path("${scanpy_plot18}")
            path("${scanpy_plot19}")
            path("${scanpy_plot20}")
            path("${scanpy_plot21}")
            path("${scanpy_plot22}")
            path("${scanpy_plot23}")
            path("${scanpy_plot24}")
            path("${scanpy_plot25}")
            path("${scanpy_plot26}")
            path("${scanpy_plot27}")
            path("${scanpy_plot28}")
            path("${scanpy_out}")
    
    script:
    """
    scanpy_ftskin.py 
    """
}

process GENERATE_HTML {

    publishDir "output", mode: "copy"

    input: path(html_in)
           path(html_template_dir)

    output: html_out //eventually want whole directory 
    
    script:
    """
    generate_html.py --html_table ${html_in} --index_html ${html_out} --html_template_dir ${html_template_dir}
    """
}

workflow { 
//    IRODS(irods_in,irods_in_metadata)
    SOUP(file(params.cellranger_dir)) //ensure relative 
    LOAD_DATA(SOUP.out[0], SOUP.out[1], SOUP.out[2], params.load_metadata)
    SCRUBLET(LOAD_DATA.out[1])
    SCANPY(SCRUBLET.out[0])
    GENERATE_HTML(html_in, html_template_dir)
}
