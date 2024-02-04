import os
import csv
import pandas as pd
import numpy as np
import statistics
import seaborn as sns
from scipy.stats import t
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import tikzplotlib

# set directory
os.chdir("C:/Users/ahadi/OneDrive/Desktop/ECV/final_results")

sns.set_theme(style="whitegrid")
plt.rcParams["figure.figsize"] = [12, 8]
plt.rcParams["figure.autolayout"] = True
plt.rcParams["axes.edgecolor"] = "black"
plt.rcParams["axes.linewidth"] = 1.5

#%% ################################# Single objective  ######################################
#################### Trait 1 _ Proportion of desirable alleles
output = pd.read_csv('output_final_SO.csv')

chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Proportion of Desirable Alleles", hue="Approach",
            data=output, palette="Set1", width=0.7,
            linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 15) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)
#legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Proportion of desirable alleles', fontsize=27)
# plt.savefig('SO_proportion.eps', dpi=300)
plt.savefig('Figure 1(a).pdf')
# tikzplotlib.save('SO_trait1_proportion.tex')
plt.show()

#%%##############  Trait 1 _ Phenotypic Value

chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Phenotypic Value", hue="Approach", data=output,
            palette="Set1", width=0.7,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 15) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)

# legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Phenotypic trait value', fontsize=27)
# plt.savefig('SO_pheno_values.eps', dpi=300)
plt.savefig('Figure 1(b).pdf')
# tikzplotlib.save('SO_trait1_pheno.tex')
plt.show()

#%%############### Inbreeding
output = pd.read_csv('Inbreeding_SO.csv')

chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Genomic Relatedness", hue="Approach",
            data=output, palette="Set1", width=0.7,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 15) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)

# legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Genomic relatedness', fontsize=27)
# plt.savefig('SO_inbreeding.eps', dpi=300)
plt.savefig('Figure 1(c).pdf')
# tikzplotlib.save('SO_inbreeding.tex')
plt.show()


#%%##### drawing plots of trait 2 and 3 for single objective case
###trait2 _proportion
output = pd.read_csv('output_SO_trait2.csv')

chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Proportion of Desirable Alleles", hue="Approach",
            data=output, palette={"lightseagreen"}, width=0.4,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
plt.xlabel('Generation',fontsize=27)
plt.ylabel('Proportion of desirable alleles', fontsize=27)
ax1.get_legend().remove()
# plt.savefig('SO_proportion_trait2_case.eps', dpi=300)
plt.savefig('Figure 2(a).pdf')
plt.show()


###trait2 _pheno
chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Phenotypic Value", hue="Approach",
            data=output, palette={"lightseagreen"}, width=0.4,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
plt.xlabel('Generation',fontsize=27)
plt.ylabel('Phenotypic trait value', fontsize=27)
ax1.get_legend().remove()
# plt.savefig('SO_pheno_trait2_case.eps', dpi=300)
plt.savefig('Figure 2(b).pdf')
plt.show()

#%%
###trait3_proportion
output = pd.read_csv('output_SO_trait3.csv')

chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Proportion of Desirable Alleles", hue="Approach",
            data=output, palette={"lightseagreen"}, width=0.4,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
plt.xlabel('Generation',fontsize=27)
plt.ylabel('Proportion of desirable alleles', fontsize=27)
ax1.get_legend().remove()
# plt.savefig('SO_proportion_trait3_case.eps', dpi=300)
plt.savefig('Figure 2(c).pdf')
plt.show()


###trait2 _pheno
chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Phenotypic Value", hue="Approach",
            data=output, palette={"lightseagreen"}, width=0.4,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
plt.xlabel('Generation',fontsize=27)
plt.ylabel('Phenotypic trait value', fontsize=27)
ax1.get_legend().remove()
# plt.savefig('SO_pheno_trait3_case.eps', dpi=300)
plt.savefig('Figure 2(d).pdf')
plt.show()

#%% ################################# Multi objective ###################################

############ Trait 3
output = pd.read_csv('output_final_MO_trait3.csv')

##prop
chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Proportion of Desirable Alleles", hue="Approach",
            data=output, palette="Set1", width=0.7,
            linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 12) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)
# legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Proportion of desirable alleles', fontsize=27)
# plt.savefig('MO_proportion_trait3.eps', dpi=300)
plt.savefig('Figure 3(a).pdf')
plt.show()

## pheno
chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Phenotypic Value", hue="Approach", data=output,
            palette="Set1", width=0.7,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 12) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)

# legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Phenotypic trait value', fontsize=27)
# plt.savefig('MO_pheno_values_trait3.eps', dpi=300)
plt.savefig('Figure 3(b).pdf')
plt.show()

#%%########## Trait 1

output = pd.read_csv('output_final_MO_trait1.csv')

##prop
chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Proportion of Desirable Alleles", hue="Approach",
            data=output, palette="Set1", width=0.7,
            linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 12) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)
# legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Proportion of desirable alleles', fontsize=27)
# plt.savefig('MO_proportion_trait1.eps', dpi=300)
plt.savefig('Figure 3(c).pdf')
plt.show()

## pheno
chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Phenotypic Value", hue="Approach", data=output,
            palette="Set1", width=0.7,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 12) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)

# legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Phenotypic trait value', fontsize=27)
# plt.savefig('MO_pheno_values_trait1.eps', dpi=300)
plt.savefig('Figure 3(d).pdf')
plt.show()

#%%########## Trait 2

output = pd.read_csv('output_final_MO_trait2.csv')

##prop
chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Proportion of Desirable Alleles", hue="Approach",
            data=output, palette="Set1", width=0.7,
            linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 12) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)
# legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Proportion of desirable alleles', fontsize=27)
# plt.savefig('MO_proportion_trait2.eps', dpi=300)
plt.savefig('Figure 3(e).pdf')
plt.show()

## pheno
chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Phenotypic Value", hue="Approach", data=output,
            palette="Set1", width=0.7,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 12) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)

# legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Phenotypic trait value', fontsize=27)
# plt.savefig('MO_pheno_values_trait2.eps', dpi=300)
plt.savefig('Figure 3(f).pdf')
plt.show()

#%% Inbreeding
output = pd.read_csv('Inbreeding_MO.csv')

chart1, ax1 = plt.subplots()
sns.boxplot(x="Generation", y="Genomic Relatedness", hue="Approach",
            data=output, palette="Set1", width=0.7,linewidth=2)

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
	label.set_fontsize(17)
    
legend = ax1.legend(loc='upper left' , borderpad=1.5,
                    labelspacing=2,ncol=1,
                    handlelength=3,fontsize = 15) #prop={'size':15}
legend.get_frame().set_linewidth(1.5)

# legend.get_frame().set_facecolor('lightgray')
legend.legendPatch.set_edgecolor("black")

plt.xlabel('Generation', fontsize=27)
plt.ylabel('Genomic relatedness', fontsize=27)
# plt.savefig('MO_inbreeding.eps', dpi=300)
plt.savefig('Figure 4.pdf')
plt.show()











