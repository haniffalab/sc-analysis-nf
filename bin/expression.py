import pandas as pd 
#obj ={'gene_name': ['TP53', 'TNF', 'EGFR', 'VEGFA', 'APOE', 
'IL6', 'TGFB1', 'MTHFR', 'ESR1', 'AKT1'], 
'expression':[1,2,3,4,5,6,7,8,9,10]}
#data = pd.DataFrame(data= obj)
#data.to_csv('data.csv', index=False) 
data = pd.read_csv('data.csv')

len(data)

sequencing = {'sequencing' 
:['visium','xenium','visium','visium','visium','visium','visium','xenium','visium','visium']}
seq = pd.DataFrame(sequencing)
data_out= pd.concat([data,seq], axis =1)
data_out
data_out.to_csv('data_out.csv', index=False) 
