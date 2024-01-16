#!/usr/bin/env python3

import fire
import jinja2
import scanpy as sc
import pandas as pd
import numpy as np
from jinja2 import Environment, FileSystemLoader
import os
import random

def generate_html(html_table:str, index_html:str):
    root = os.path.dirname(os.path.abspath(__file__))
    templates_dir = os.path.join(root, "/Users/nlg143/projects/sc-analysis-nf/jinja/templates")
    env = Environment( loader = FileSystemLoader(templates_dir) )
    template = env.get_template('index.html')

    ######### reading the csv ######### 
    import csv
    with open(html_table, newline='') as csvfile:
        reader = csv.reader(csvfile)
        build = list(reader)
    ###################################

    filename = os.path.join("/Users/nlg143/projects/sc-analysis-nf", index_html)
    with open(filename, 'w') as fh:
        fh.write(template.render(
            h1 = "SC-analysis-nf Pipeline Output Page",
            show_preprocessing = True,
            img_1 = "/Users/nlg143/projects/sc-analysis-nf/output/figures/highest_expr_genes_2.png",
            img_2 = "/Users/nlg143/projects/sc-analysis-nf/output/figures/violin_pct_counts_mt_2.png",
            img_3 = "/Users/nlg143/projects/sc-analysis-nf/output/figures/violin_n_genes_by_counts_2.png",
            img_4 = "/Users/nlg143/projects/sc-analysis-nf/output/figures/violin_total_counts_2.png",
            img_5 = "/Users/nlg143/projects/sc-analysis-nf/output/figures/scatter_pct_counts_mt_2.png",
            img_6 = "/Users/nlg143/projects/sc-analysis-nf/output/figures/scatter_n_genes_by_counts_2.png",
            img_7 = "/Users/nlg143/projects/sc-analysis-nf/output/figures/filter_genes_dispersionscatter_2.png",
            JOBID = random.randint(10, 10000),
            mytable = build,
        ))


if __name__ == "__main__":
    fire.Fire(generate_html)
