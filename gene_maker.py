### loading packages
import random
import sys
import csv

"""
Output of code:
    1. genes_output.txt: includes generated gene information required for simulating initial population
    2. markers_index.csv: index of genes corresponding to markers
    3. trait1_index.csv: index of genes corresponding to trait 1
    2. trait2_index.csv: index of genes corresponding to trait 2
    2. trait3_index.csv: index of genes corresponding to trait 3
"""


sys.stdout = open("genes_output.txt", "w")    

gene_pool=[] 
markers_index_list=[]
trait_1_index_list=[]
trait_2_index_list=[]
trait_3_index_list=[]
 
"""
 description of genes:
 1. marker
 
 21. QTL for trait 1 only-in epistasis
 22. QTL for trait 1 only-not in epistasis
 
 3. QTL for trait 2 only
 
 4. QTL for trait 3 only
 
 5. QTL for trait 1 and 3 in common
"""

for i in range(100):
    gene_pool.append(1)
for i in range(10):
    gene_pool.append(21)
for i in range(10):
    gene_pool.append(22)
for i in range(10):
    gene_pool.append(3)
for i in range(50):
    gene_pool.append(4)
for i in range(20):
    gene_pool.append(5)
    
random.shuffle(gene_pool)

for i in range(len(gene_pool)):
    if gene_pool[i]== 1:
        markers_index_list.append(i)
    elif gene_pool[i]== 21:
        trait_1_index_list.append(i)
    elif gene_pool[i]== 22:
        trait_1_index_list.append(i)
    elif gene_pool[i]== 3:
        trait_2_index_list.append(i) 
    elif gene_pool[i]== 4:
        trait_3_index_list.append(i)         
    elif gene_pool[i]== 5:
        trait_1_index_list.append(i)         
        trait_3_index_list.append(i)


marker_used = 0
QTL_used = 0  #change based on chromosome
  

number_of_genes=len(gene_pool)

b=(number_of_genes+1)/2
c=(number_of_genes-1)/2

for i in range(1,number_of_genes+1): #change based on chromosome
    
    if i <= b :
        x = (number_of_genes+1)- i
    else:
        x = i
    r1 = (x-b)**2/(2*c**2) #change base for different chromosome
    
    mu=0.1
    
    r2 = ((0.5-mu)/c)*x+(0.5-((0.5-mu)*number_of_genes/c))
    
    r = (r1+r2)/2
    r = round( r, 2)
    
    epsilon=round(random.gauss(0, 0.01), 3)
    
    gene_type = gene_pool.pop(0)
    
    if gene_type == 1 :
        marker_used = marker_used +1
        print(i) #change based on chromosome
        print("M01{}".format(marker_used)) # change third char to chromosome index
        print("1 {} 2 0".format(r)) # change first char to chromosome index
        print("    0  2  0")
        print("")
    elif gene_type == 21 :
        QTL_used = QTL_used+1
        print(i) #change based on chromosome
        print("L000{}".format(QTL_used))
        print("1 {} 2 1".format(r)) # change first char to chromosome index
        print("        1 1 2 1") # first char shows the trait index
        print("")
    elif gene_type == 22 :
        QTL_used = QTL_used+1
        print(i) #change based on chromosome
        print("L000{}".format(QTL_used))
        print("1 {} 2 1".format(r)) # change first char to chromosome index
        print("        1 1 1 0 1 {} 0".format(0.5+epsilon)) # first char shows the trait index
        print("")
    elif gene_type == 3 :
        QTL_used = QTL_used+1        
        print(i) #change based on chromosome
        print("L000{}".format(QTL_used))
        print("1 {} 2 1".format(r)) # change first char to chromosome index
        print("        2 1 1 0 1 {} 0".format(0.5+epsilon)) # first char shows the trait index
        print("")
    elif gene_type == 4:
        QTL_used = QTL_used+1        
        print(i) #change based on chromosome
        print("L000{}".format(QTL_used))
        print("1 {} 2 1".format(r)) # change first char to chromosome index
        print("        3 1 1 0 1 {} 0".format(0.5+epsilon)) # first char shows the trait index
        print("")
    elif gene_type == 5:
        QTL_used = QTL_used+1
        print(i) #change based on chromosome
        print("L000{}".format(QTL_used))
        print("1 {} 2 2".format(r)) # change first char to chromosome index
        print("        1 1 1 0 1 {} 0".format(0.5+epsilon)) # first char shows the  first trait index        
        print("        3 1 1 0 0 {} 1".format(0.5+epsilon)) # first char shows the second rait index
        print("")
        
        
with open("markers_index.csv",'w',newline='') as fp:
    a=csv.writer(fp,delimiter=',')
    a.writerow(markers_index_list)

with open("trait1_index.csv",'w',newline='') as fp:
    a=csv.writer(fp,delimiter=',')
    a.writerow(trait_1_index_list)

with open("trait2_index.csv",'w',newline='') as fp:
    a=csv.writer(fp,delimiter=',')
    a.writerow(trait_2_index_list)
    
with open("trait3_index.csv",'w',newline='') as fp:
    a=csv.writer(fp,delimiter=',')
    a.writerow(trait_3_index_list)

sys.stdout.close()

