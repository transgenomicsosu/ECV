import os
import pandas as pd
import csv


######################################> Single objective
##########################>trait 1
output=[]
output.append(['Approach', 'Generation', 'Replication', 'Proportion of Desirable Alleles', 'Phenotypic Value', 'Inbreeding'])

#reading ECV outputs
for i in range(30):
    ## path for the folder corresponding populations.
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results/pop{}".format(i))
    with open("output_ecv.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('ECV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[3].split(",")[j])
            row.append(data[7].split(",")[j])
            if j <= 3:
                row.append(data[1].split(",")[j])
            else:
                row.append("-")
            output.append(row)

#reading phenotypic selection outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results/pop{}".format(i))
    with open("output_pheno.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('Phenotypic Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[3].split(",")[j])
            row.append(data[7].split(",")[j])
            if j <= 3:
                row.append(data[1].split(",")[j])
            else:
                row.append("-")
            output.append(row)

#reading GEBV selection outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_GEBV_SO/pop{}".format(i))
    with open("output_GEBV_SO.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('GEBV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[3].split(",")[j])
            row.append(data[7].split(",")[j])
            if j <= 3:
                row.append(data[1].split(",")[j])
            else:
                row.append("-")
            output.append(row)

os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results") ## path to save final results
with open("output_final_SO.csv",'w',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerows(output)
    
##########################>trait 2
output=[]
output.append(['Approach', 'Generation', 'Replication', 'Proportion of Desirable Alleles', 'Phenotypic Value', 'Inbreeding'])

#reading ECV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_ecv.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('ECV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[4].split(",")[j])
            row.append(data[8].split(",")[j])
            if j <= 3:
                row.append(data[1].split(",")[j])
            else:
                row.append("-")
            output.append(row)
 
os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results")
with open("output_SO_trait2.csv",'w',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerows(output) 

##########################>trait 3
output=[]
output.append(['Approach', 'Generation', 'Replication', 'Proportion of Desirable Alleles', 'Phenotypic Value', 'Inbreeding'])

#reading ECV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_ecv.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('ECV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[5].split(",")[j])
            row.append(data[9].split(",")[j])
            if j <= 3:
                row.append(data[1].split(",")[j])
            else:
                row.append("-")
            output.append(row)
 
os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results")
with open("output_SO_trait3.csv",'w',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerows(output)    
    
    
######################################> Single objective
##########################> Inbreeding
output=[]
output.append(['Approach', 'Generation', 'Replication', 'Inbreeding'])

#reading ECV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_ecv.csv") as file:
        data=file.read().split("\n")
        for j in range(4):
            row=[]
            row.append('ECV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[1].split(",")[j])
            output.append(row)

#reading phenotypic selection outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_pheno.csv") as file:
        data=file.read().split("\n")
        for j in range(4):
            row=[]
            row.append('Phenotypic Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[1].split(",")[j])
            output.append(row)

#reading GEBV selection outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_GEBV_SO/pop{}".format(i))
    with open("output_GEBV_SO.csv") as file:
        data=file.read().split("\n")
        for j in range(4):
            row=[]
            row.append('GEBV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[1].split(",")[j])
            output.append(row)

os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/SO_results")
with open("Inbreeding.csv",'w',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerows(output)


######################################> Multi-objective
##########################>trait 3
output=[]
output.append(['Approach', 'Generation', 'Replication', 'Proportion of Desirable Alleles', 'Phenotypic Value', 'Inbreeding'])

#reading ECV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_ecv.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('ECV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[7].split(",")[j])
            row.append(data[11].split(",")[j])
            if j <= 3:
                row.append(data[3].split(",")[j])
            else:
                row.append("-")
            output.append(row)

#reading phenotypic selection outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_pheno.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('Phenotypic Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[7].split(",")[j])
            row.append(data[11].split(",")[j])
            if j <= 3:
                row.append(data[3].split(",")[j])
            else:
                row.append("-")
            output.append(row)


#reading GEBV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_GEBV/pop{}".format(i))
    with open("output_GEBV_generation0.csv") as file:
        row=[]
        row.append('GEBV Selection')
        row.append(0)
        row.append(i+1)        
        data=file.read().split("\n")
        row.append(data[3].split(",")[0])
        row.append(data[7].split(",")[0])        
    with open("Inbreeding_gen0.csv") as file:
        data=file.read().split("\n")
        row.append(data[0].split(",")[0])
    output.append(row)    
    row=[]

    with open("output_GEBV_gen1.csv") as file:
        row=[]
        row.append('GEBV Selection')
        row.append(1)
        row.append(i+1)        
        data=file.read().split("\n")
        row.append(data[6].split(",")[0])
        row.append(data[10].split(",")[0])        
    with open("Inbreeding_gen1.csv") as file:
        data=file.read().split("\n")
        row.append(data[0].split(",")[0])
    output.append(row)    
    row=[]            
        
    with open("output_GEBV_gen2-4.csv") as file:
        data=file.read().split("\n")
        for j in range(3):
            row=[]
            row.append('GEBV Selection')
            row.append(j+2)
            row.append(i+1)
            row.append(data[7].split(",")[j])
            row.append(data[11].split(",")[j])
            if j <= 1:
                row.append(data[3].split(",")[j])
            else:
                row.append("-")
            output.append(row)
####
os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results")
with open("output_final_MO_trait3.csv",'w',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerows(output)

######################################> Multi-objective
##########################>trait 1
output=[]
output.append(['Approach', 'Generation', 'Replication', 'Proportion of Desirable Alleles', 'Phenotypic Value', 'Inbreeding'])

#reading ECV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_ecv.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('ECV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[5].split(",")[j])
            row.append(data[9].split(",")[j])
            if j <= 3:
                row.append(data[3].split(",")[j])
            else:
                row.append("-")
            output.append(row)

#reading phenotypic selection outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_pheno.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('Phenotypic Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[5].split(",")[j])
            row.append(data[9].split(",")[j])
            if j <= 3:
                row.append(data[3].split(",")[j])
            else:
                row.append("-")
            output.append(row)


#reading GEBV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_GEBV/pop{}".format(i))
    with open("output_GEBV_generation0.csv") as file:
        row=[]
        row.append('GEBV Selection')
        row.append(0)
        row.append(i+1)        
        data=file.read().split("\n")
        row.append(data[1].split(",")[0])
        row.append(data[5].split(",")[0])        
    with open("Inbreeding_gen0.csv") as file:
        data=file.read().split("\n")
        row.append(data[0].split(",")[0])
    output.append(row)    
    row=[]

    with open("output_GEBV_gen1.csv") as file:
        row=[]
        row.append('GEBV Selection')
        row.append(1)
        row.append(i+1)        
        data=file.read().split("\n")
        row.append(data[4].split(",")[0])
        row.append(data[8].split(",")[0])        
    with open("Inbreeding_gen1.csv") as file:
        data=file.read().split("\n")
        row.append(data[0].split(",")[0])
    output.append(row)    
    row=[]            
        
    with open("output_GEBV_gen2-4.csv") as file:
        data=file.read().split("\n")
        for j in range(3):
            row=[]
            row.append('GEBV Selection')
            row.append(j+2)
            row.append(i+1)
            row.append(data[5].split(",")[j])
            row.append(data[9].split(",")[j])
            if j <= 1:
                row.append(data[3].split(",")[j])
            else:
                row.append("-")
            output.append(row)


#####
os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results")
with open("output_final_MO_trait1.csv",'w',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerows(output)

######################################> Multi-objective
##########################>trait 2
output=[]
output.append(['Approach', 'Generation', 'Replication', 'Proportion of Desirable Alleles', 'Phenotypic Value', 'Inbreeding'])

#reading ECV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_ecv.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('ECV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[6].split(",")[j])
            row.append(data[10].split(",")[j])
            if j <= 3:
                row.append(data[3].split(",")[j])
            else:
                row.append("-")
            output.append(row)

#reading phenotypic selection outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_pheno.csv") as file:
        data=file.read().split("\n")
        for j in range(5):
            row=[]
            row.append('Phenotypic Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[6].split(",")[j])
            row.append(data[10].split(",")[j])
            if j <= 3:
                row.append(data[3].split(",")[j])
            else:
                row.append("-")
            output.append(row)
            
#reading GEBV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_GEBV/pop{}".format(i))
    with open("output_GEBV_generation0.csv") as file:
        row=[]
        row.append('GEBV Selection')
        row.append(0)
        row.append(i+1)        
        data=file.read().split("\n")
        row.append(data[2].split(",")[0])
        row.append(data[6].split(",")[0])        
    with open("Inbreeding_gen0.csv") as file:
        data=file.read().split("\n")
        row.append(data[0].split(",")[0])
    output.append(row)    
    row=[]

    with open("output_GEBV_gen1.csv") as file:
        row=[]
        row.append('GEBV Selection')
        row.append(1)
        row.append(i+1)        
        data=file.read().split("\n")
        row.append(data[5].split(",")[0])
        row.append(data[9].split(",")[0])        
    with open("Inbreeding_gen1.csv") as file:
        data=file.read().split("\n")
        row.append(data[0].split(",")[0])
    output.append(row)    
    row=[]            
        
    with open("output_GEBV_gen2-4.csv") as file:
        data=file.read().split("\n")
        for j in range(3):
            row=[]
            row.append('GEBV Selection')
            row.append(j+2)
            row.append(i+1)
            row.append(data[6].split(",")[j])
            row.append(data[10].split(",")[j])
            if j <= 1:
                row.append(data[3].split(",")[j])
            else:
                row.append("-")
            output.append(row)

os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results")
with open("output_final_MO_trait2.csv",'w',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerows(output)

    
######################################> Multi-objective
##########################> Inbreeding
output=[]
output.append(['Approach', 'Generation', 'Replication', 'Inbreeding'])

#reading ECV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_ecv.csv") as file:
        data=file.read().split("\n")
        for j in range(4):
            row=[]
            row.append('ECV Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[3].split(",")[j])
            output.append(row)

#reading phenotypic selection outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results/pop - Copy ({}) - Copy".format(i))
    with open("output_pheno.csv") as file:
        data=file.read().split("\n")
        for j in range(4):
            row=[]
            row.append('Phenotypic Selection')
            row.append(j)
            row.append(i+1)
            row.append(data[3].split(",")[j])
            output.append(row)

#reading GEBV outputs
for i in range(30):
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_GEBV/pop{}".format(i))
    row=[]
    row.append('GEBV Selection')
    row.append(0)
    row.append(i+1)
    with open("Inbreeding_gen0.csv") as file:
        data=file.read().split("\n")
        row.append(data[0].split(",")[0])
    output.append(row)    
    row=[]    
    
    row.append('GEBV Selection')
    row.append(1)
    row.append(i+1)            
    with open("Inbreeding_gen1.csv") as file:
        data=file.read().split("\n")
        row.append(data[0].split(",")[0])
    output.append(row)    
    row=[]            
        
    with open("output_GEBV_gen2-4.csv") as file:
        data=file.read().split("\n")
        for j in range(2):
            row=[]
            row.append('GEBV Selection')
            row.append(j+2)
            row.append(i+1)
            row.append(data[3].split(",")[j])
            output.append(row)


os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_results")
with open("Inbreeding.csv",'w',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerows(output)
