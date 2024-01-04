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
        img_1 = "/Users/nlg143/projects/sc-analysis-nf/output/FCAImmP7241240/soup.channel.png",
        img_2 = "/Users/nlg143/projects/sc-analysis-nf/output/FCAImmP7241240/scrublet_histogram.png",
        img_3 = "/Users/nlg143/projects/sc-analysis-nf/output/FCAImmP7241240/scrublet_umap.png",
        img_4 = "/Users/nlg143/projects/sc-analysis-nf/output/FCAImmP7241240/qcmetrics_totalcounts.png",
        img_5 = "/Users/nlg143/projects/sc-analysis-nf/output/FCAImmP7241240/figures/violin_qcmetrics_pctcountsmt.png",
        img_6 = "/Users/nlg143/projects/sc-analysis-nf/output/FCAImmP7241240/figures/scatter_qcmetrics.png",
        img_7 = "/Users/nlg143/projects/sc-analysis-nf/output/FCAImmP7241240/figures/violin_all_counts.png",
        JOBID = random.randint(10, 10000),
        mytable = build,
    ))



#####add this back into template index.html if you want a if statemnt in there 
#<ul>
#  {% for jobid in jobid %}
#     <p>{{ jobid }}</p>
#  {% endfor %}
#  </ul>