
import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import pandas as pd
from pandas.plotting import scatter_matrix
import seaborn as sns


len_of_prot_seq = 584
#print("lenght of sequnce ", len_of_prot_seq)

#open output with numbers of interaction residues *.CSV file -> still buggy...
with open('ligads_outputs/ALL.csv','r') as f:
	reader = csv.reader(f)
	input_csv = list(reader)

result = {}

#with open('matiX_cons.csv','r') as f:
#	reader = csv.reader(f)
#	input_csv_con = list(reader)

#print(input_csv)
for y in input_csv:
	sequence = list(map(int,y[2:]))
	#print(y[0])
	temp = [0]*len_of_prot_seq
	for val in sequence:
		temp[val] = 1
		#print(temp)
	result[y[0]] = temp
#print(val, temp)


#print result

dataset = pd.DataFrame(result)
#dataset_con = pd.DataFrame(cons)
datasett= pd.DataFrame.transpose(dataset)
#dataset_con1 = pd.DataFrame.transpose(dataset_con)
#datasett_con= pd.DataFrame.transpose(dataset_con)
datasett.drop(datasett.columns[[0]], axis=1, inplace=True) #removing column "0", i could be not created before... ::
#binary matrix
#id_labels = datasett.columns[1:]
#print(datasett)

#print(dataset_con1)

#esult ["cons"]= [dataset_con1]
#print(datasett_con)
#ataset_merg = pd.concat(frames)#
#print(dataset_merg)


#datasett.to_csv('matiX.csv', float_format='%.2f', na_rep="NAN!")
#print(result)

plt.figure(figsize=(120, 50))

sns.set()
sns.set(font_scale=2)
#print(datasett[10:])
domain_I = pd.DataFrame(result)
domain_I=pd.DataFrame.transpose(domain_I)
domain_I.drop(domain_I.columns[[0]], axis=1, inplace=True)
#print (domain_I)
domain_II = pd.DataFrame(result)
domain_II=pd.DataFrame.transpose(domain_II)
# print(domain_I)
domain_III = pd.DataFrame(result)
domain_III=pd.DataFrame.transpose(domain_III)

#cutting matix on HSA domains DI - 1-196 -> split plot on domains by number in sequnce
domain_I.drop(domain_I.columns[[range(196,583)]], axis=1, inplace=True)
domain_II.drop(domain_II.columns[[range(0,197)]], axis=1, inplace=True)
domain_IIa=pd.DataFrame(domain_II) # - new column counting
domain_IIa.drop(domain_IIa.columns[[range(188,386)]], axis=1, inplace=True) #removing column "0", i could be not created before... ::
#print(domain_IIa)
domain_III.drop(domain_III.columns[[range(0,385)]], axis=1, inplace=True) #removing column "0", i could be not created before... ::

#ax = sns.heatmap(datasett,linewidths=0.05, linecolor='black', cbar=None, square=1, annot_kws={'size': 0}, cmap= 'coolwarm').get_figure().savefig('output.png')
ax = sns.heatmap(domain_I,linewidths=0.05, linecolor='black', cbar=None, square=1, annot_kws={'size': 0}, cmap= 'Greys').get_figure().savefig('output_DI.png')
ax = sns.heatmap(domain_IIa,linewidths=0.05, linecolor='black', cbar=None, square=1, annot_kws={'size': 0}, cmap= 'Greys').get_figure().savefig('output_DII.png')
ax = sns.heatmap(domain_III,linewidths=0.05, linecolor='black', cbar=None, square=1, annot_kws={'size': 0}, cmap= 'Greys').get_figure().savefig('output_DIII.png')


#sns.heatmap(datasett)
#plt.show()
