from jinja2 import Environment, FileSystemLoader
import os
import random



root = os.path.dirname(os.path.abspath(__file__))
templates_dir = os.path.join(root, '/Users/nlg143/projects/sc-analysis-nf/jinja/templates')
env = Environment( loader = FileSystemLoader(templates_dir) )
template = env.get_template('index.html')

######### reading the csv ######### 
import csv
with open('html_csv_test.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    build = list(reader)
###################################

filename = os.path.join(root, '/Users/nlg143/projects/sc-analysis-nf/jinja/html', 'index.html')
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



#####add this back into template index.html if you want a if statemnt in there 
#<ul>
#  {% for jobid in jobid %}
#     <p>{{ jobid }}</p>
#  {% endfor %}
#  </ul>