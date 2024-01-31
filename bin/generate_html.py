#!/usr/bin/env python3

import fire
import jinja2
from jinja2 import Environment, FileSystemLoader
import os
import random


def generate_html(html_table: str, index_html: str, html_template_dir: str):
    env = Environment(
        loader=FileSystemLoader(
            os.path.join(html_template_dir, "theme/jinja_templates")
        )
    )
    template_index = env.get_template("index.html")
    # template_soup = env.get_template("soup.html")

    ######### reading the csv #########
    import csv

    with open(html_table, newline="") as csvfile:
        reader = csv.reader(csvfile)
        build = list(reader)
    ###################################

    output_index = "html/theme/dist/index.html"
    with open(output_index, "w") as fh:
        fh.write(
            template_index.render(
                SAMPLE_ID="SC-analysis-nf Pipeline Output Page",
                show_preprocessing=True,
                img_1="/sc-analysis-nf/output/figures/highest_expr_genes_2.png",
                img_2="/sc-analysis-nf/output/figures/violin_pct_counts_mt_2.png",
                img_3="/sc-analysis-nf/output/figures/violin_n_genes_by_counts_2.png",
                img_4="/sc-analysis-nf/output/figures/violin_total_counts_2.png",
                img_5="/sc-analysis-nf/output/figures/scatter_pct_counts_mt_2.png",
                img_6="/sc-analysis-nf/output/figures/scatter_n_genes_by_counts_2.png",
                img_7="/sc-analysis-nf/output/figures/filter_genes_dispersionscatter_2.png",
                JOBID=random.randint(10, 10000),
                mytable=build,
            )
        )


if __name__ == "__main__":
    fire.Fire(generate_html)
