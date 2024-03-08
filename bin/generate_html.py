#!/usr/bin/env python3

import fire
import jinja2
from jinja2 import Environment, FileSystemLoader
import os
import random
import shutil


def generate_html(html_table: str):
    env = Environment(loader=FileSystemLoader("dist"))
    template_index = env.get_template("index.html")
    # template_soup = env.get_template("soup.html")

    ######### reading the csv #########
    # import csv

    # with open(html_table, newline="") as csvfile:
    #     reader = csv.reader(csvfile)
    #     build = list(reader)
    ###################################

    shutil.copytree("dist", "report")
    output_index = "report/index.html"
    with open(output_index, "w") as fh:
        fh.write(
            template_index.render(
                SAMPLE_ID="SC-analysis-nf Pipeline Output Page",
                show_preprocessing=True,
                PLOT_SOUPX_CHANNEL="/plots/soupx/soup.channel.png",
                img_2="/sc-analysis-nf/output/figures/violin_pct_counts_mt_2.png",
                img_3="/sc-analysis-nf/output/figures/violin_n_genes_by_counts_2.png",
                img_4="/sc-analysis-nf/output/figures/violin_total_counts_2.png",
                img_5="/sc-analysis-nf/output/figures/scatter_pct_counts_mt_2.png",
                img_6="/sc-analysis-nf/output/figures/scatter_n_genes_by_counts_2.png",
                img_7="/sc-analysis-nf/output/figures/filter_genes_dispersionscatter_2.png",
                JOBID=random.randint(10, 10000),
            )
        )


if __name__ == "__main__":
    fire.Fire(generate_html)
