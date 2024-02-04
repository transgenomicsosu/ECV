import csv
from gurobipy import *
import numpy as np 
import datetime
import os
import gurobipy as gp
import statistics

################################
##    
###
####  Initialization and reading
###
##
################################

n_traits=3
epsilon=0.25
generation_number=5 ### add one more to the generations you desire
cross_list=[50,10,5,5,5] ### last component is dummy
K=[10000, 5000, 1000, 500, 500, 500]  ### last component is dummy


tolerance_1=[0.17, 0.05 , 0.05 , 0.05]# tolerance for first objective
tolerance_2=[0   ,  0   , 0    , 0.05]# tolerance for second objective


"""
desired order of traits for optimization: T3>T1>T2
"""

markers_index=[]
trait1_index=[]
trait2_index=[]
trait3_index=[]

begin_time = datetime.datetime.now()

with open("markers_index.csv","r") as f:
    file = csv.reader(f, delimiter=',')
    for row in file:
        for a in row:
            markers_index.append(int(a))


with open("trait1_index.csv","r") as f:
    file = csv.reader(f, delimiter=',')
    for row in file:
        for a in row:
            trait1_index.append(int(a))

with open("trait2_index.csv","r") as f:
    file = csv.reader(f, delimiter=',')
    for row in file:
        for a in row:
            trait2_index.append(int(a))

with open("trait3_index.csv","r") as f:
    file = csv.reader(f, delimiter=',')
    for row in file:
        for a in row:
            trait3_index.append(int(a))


trait_1_3_index=[]

for a in trait1_index:
    if a in trait3_index:
        trait_1_3_index.append(a)

N1=len(trait1_index)
N2=len(trait2_index)
N3=len(trait3_index)
 
N4=len(trait_1_3_index)
   
number_of_genes= len(markers_index)+N1+N2+N3-N4

###DEfining call back
def call_back(model, where):
    if where == GRB.Callback.MIPSOL:
        # make a list of edges selected in the solution
        vals = model.cbGetSolution(model._t)
        indices = gp.tuplelist((m,k)  for m, k in model._t.keys() if vals[m, k] == 1)
        if G[indices[0][1]][indices[1][1]] >= epsilon:
            model.cbLazy(gp.quicksum(model._t[m, k]  for m, k in indices) <= 1)
            

Pheno_allele_freq_trait1=[]
Pheno_phenotype_trait1=[]
Pheno_allele_freq_trait2=[]
Pheno_phenotype_trait2=[]
Pheno_allele_freq_trait3=[]
Pheno_phenotype_trait3=[]

Pheno_Inbreeding=[]

ECV_allele_freq_trait1=[]
ECV_phenotype_trait1=[]
ECV_allele_freq_trait2=[]
ECV_phenotype_trait2=[]
ECV_allele_freq_trait3=[]
ECV_phenotype_trait3=[]

ECV_Inbreeding=[]

p=[] #vector of frequencies
P=[] #matrix of frequencies

p_ecv=[]
p_pheno=[]



################################
##    
###
####  ECV Section
###
##
################################  
    
for g in range(generation_number):
    if g == 0:
        with open('Error.pop','r') as file:
            data=file.read().split("\n")
            data=data[18:]
            m = len(data)
            m=int((m+1)/3)
    
            init1=[]
            init2=[]
            init3=[]
                    
            X=[]
            Xc=[] #centered X  
            for i in range(K[g]):
                
                list1=[]
                list2=[]
               
                list1=[2-int(a)  for a in data[3*i][:number_of_genes]]
            
                list2=[2-int(a)  for a in data[(3*i)+1][:number_of_genes]]
                    
                init1_row1=[]
                init1_row2=[]
                for j in trait1_index:
                  init1_row1.append(list1[j])
                  init1_row2.append(list2[j])     
                init1.append(init1_row1)
                init1.append(init1_row2)

                init2_row1=[]
                init2_row2=[]
                for j in trait2_index:
                  init2_row1.append(list1[j])
                  init2_row2.append(list2[j])
                init2.append(init2_row1)
                init2.append(init2_row2)
                               
                init3_row1=[]
                init3_row2=[]                
                for j in trait3_index:
                    if j in trait_1_3_index:
                      init3_row1.append(1-list1[j])
                      init3_row2.append(1-list2[j]) 
                    else:
                      init3_row1.append(list1[j])
                      init3_row2.append(list2[j])              
                init3.append(init3_row1)
                init3.append(init3_row2) 
                
                #creating direct Xc matrix
                Xc_row=[]
                for j in markers_index:               
                    if list1[j]+list2[j] == 0:
                        Xc_row.append(-1)
                    elif list1[j]+list2[j] == 1:
                        Xc_row.append(0)
                    elif list1[j]+list2[j] == 2:
                        Xc_row.append(1)
                        
                Xc.append(Xc_row)
        G=[]
        p=[ (sum(row[i] for row in Xc)+K[g])/(2*K[g]) for i in range(len(markers_index))]
        p_ecv.append(p)        
        P=[2*a-1 for a in p]        
        Z=np.array(Xc)-np.array(P)        
        Z_T = Z.T       
        G=np.dot(Z,Z_T)        
        d=2*sum(p[k]*(1-p[k])  for k in range(len(p)))        
        G = G/d

        p_ecv.append('g{}'.format(g))
        for e in range(cross_list[g]):                   
            p=[(sum(Xc[f][i] for f in range(100*e, 100*e+99)) +100)/(2*100) for i in range(len(markers_index))]
            p_ecv.append(p)
                
               
        ##trait1   
        P1=[]
        P1_table=[]
        P1_table_row=[]
        
        for i in range(N1):
            P1_table=[]
            for j in range(2):
                P1_table_row=[]
                for k in range(K[g]):
                    P1_table_row.append(init1[2*k+j][i])
                P1_table.append(P1_table_row)
            P1.append(P1_table)            
        H1=[sum(P1[i][j][k] for i in range(N1) for j in range(2))   for k in range(K[g])]                     
        ECV_allele_freq_trait1.append(statistics.mean(H1)/(2*N1)) #### reporting

        #trait2
        P2=[]
        P2_table=[]
        P2_table_row=[]
        
        for i in range(N2):
            P2_table=[]
            for j in range(2):
                P2_table_row=[]
                for k in range(K[g]):
                    P2_table_row.append(init2[2*k+j][i])
                P2_table.append(P2_table_row)
            P2.append(P2_table)            
        H2=[sum(P2[i][j][k] for i in range(N2) for j in range(2))   for k in range(K[g])] 
        ECV_allele_freq_trait2.append(statistics.mean(H2)/(2*N2)) #### reporting
        
        #trait3
        P3=[]
        P3_table=[]
        P3_table_row=[]
        
        for i in range(N3):
            P3_table=[]
            for j in range(2):
                P3_table_row=[]
                for k in range(K[g]):
                    P3_table_row.append(init3[2*k+j][i])
                P3_table.append(P3_table_row)
            P3.append(P3_table)            
        H3=[sum(P3[i][j][k] for i in range(N3) for j in range(2))   for k in range(K[g])]                     
        ECV_allele_freq_trait3.append(statistics.mean(H3)/(2*N3)) #### reporting   
        
        ###################### optimization #####################
        conflict_solutions=[]
        Female=[]
        Male=[]
        inbreeding_list=[]
        for q in range(cross_list[g]): 
            mod = Model()                           
            t=mod.addVars(2,K[g], vtype=GRB.BINARY, name="t")            
            mod.addConstrs(sum(t[m,k] for k in range(K[g])) == 1 for m in range(2))
            
            for row in conflict_solutions:
                mod.addConstr(t[0,row[0]]+t[1,row[1]] <= 1)
                mod.addConstr(t[1,row[0]]+t[0,row[1]] <= 1)

            obj1=0.25*quicksum(t[m,k]*H1[k] for k in range(K[g]) for m in range(2))
            obj2=0.25*quicksum(t[m,k]*H2[k] for k in range(K[g]) for m in range(2))
            obj3=0.25*quicksum(t[m,k]*H3[k] for k in range(K[g]) for m in range(2)) 
             
            mod.setObjectiveN(obj3, 0 ,  priority=2, reltol=tolerance_1[g]) #abstol=  , reltol=0.05 
            mod.setObjectiveN(obj1, 1 ,  priority=1, reltol=tolerance_2[g])
            mod.setObjectiveN(obj2, 2 ,  priority=0 )  
                            
            mod.ModelSense = GRB.MAXIMIZE
                
            mod.Params.OutputFlag=0
            mod.setParam(GRB.Param.Threads,32)            
            mod._t = t            
            mod.Params.lazyConstraints = 1                       
            mod.optimize(call_back)
            
            conflict_solutions_row=[]
            lis=[]
            for i in range(2):
                for j in range(K[g]):
                    if t[i,j].x == 1:
                        conflict_solutions_row.append(j)
                        #print("Inividual {} is selected as parent{}".format(j+1,i+1))
                        if i ==0:
                            Female.append(j+1)
                        else:
                            Male.append(j+1)
                        lis.append(j)
            conflict_solutions.append(conflict_solutions_row)
            inbreeding_list.append(G[lis[0]][lis[1]])
            
        print("parents from population {} are selected.".format(g))
        ECV_Inbreeding.append(statistics.mean(inbreeding_list))   ### Reporting
        
               
        #creating pox file
        List=[]
        List.append('Model CrossID Female Male BC TC DCFem DCMale PopFemale PopMale PopBC PopTC PopDCFem PopDCMale')
        for v in range(cross_list[0]):
            List.append('1    {}     {}    {}    {}    {}    {}    {}    1    1    1    1    1    1'.format(v+1,Female[v],Male[v],Female[v],1,Female[v],Male[v]))
        
        #writing pox file
        with open("Simulation.pox", "w") as output:
            for i in range(len(List)):
                if i == range(len(List)):
                    output.writelines(str(List[i]))
                else:
                    output.writelines(str(List[i])+'\n')
    
        #running simulation to obtain next progeny
        with open("QuLinePlus.mio", "w") as output:
            output.writelines(str('Simulation.ges')+'\n')
            output.writelines(str('Error.pop')+'\n')
            output.writelines(str('Simulation_0.qmp')+'\n')
            output.writelines(str('Simulation.pox')+'\n')
            output.writelines(str('Simulation_out'))
        os.system('cmd /c "C:/Users/Pouya Ahadi/Desktop/QL.exe"')
        print("generation {} is simulated.".format(g+1))
    if g > 0:
        with open('Simulation.pox_001_001_001_001.pop','r') as file:                        
            data=file.read().split("\n")
            data=data[15:]
            m = len(data)
            m=int((m+1)/3)
    
            init1=[]
            init2=[]
            init3=[]
                    
            X=[]
            Xc=[] #centered X  
            for i in range(K[g]):
                
                list1=[]
                list2=[]
               
                list1=[2-int(a)  for a in data[3*i][:number_of_genes]]
            
                list2=[2-int(a)  for a in data[(3*i)+1][:number_of_genes] ]
                    
                init1_row1=[]
                init1_row2=[]
                for j in trait1_index:
                  init1_row1.append(list1[j])
                  init1_row2.append(list2[j])     
                init1.append(init1_row1)
                init1.append(init1_row2)

                init2_row1=[]
                init2_row2=[]
                for j in trait2_index:
                  init2_row1.append(list1[j])
                  init2_row2.append(list2[j])
                init2.append(init2_row1)
                init2.append(init2_row2)
                               
                init3_row1=[]
                init3_row2=[]                
                for j in trait3_index:
                    if j in trait_1_3_index:
                      init3_row1.append(1-list1[j])
                      init3_row2.append(1-list2[j]) 
                    else:
                      init3_row1.append(list1[j])
                      init3_row2.append(list2[j])              
                init3.append(init3_row1)
                init3.append(init3_row2)
                
                #creating direct Xc matrix
                Xc_row=[]
                for j in markers_index:               
                    if list1[j]+list2[j] == 0:
                        Xc_row.append(-1)
                    elif list1[j]+list2[j] == 1:
                        Xc_row.append(0)
                    elif list1[j]+list2[j] == 2:
                        Xc_row.append(1)
                        
                Xc.append(Xc_row)
      

        G=[]
        Z=np.array(Xc)-np.array(P)        
        Z_T = Z.T       
        G=np.dot(Z,Z_T)        
        d=2*sum(p[k]*(1-p[k])  for k in range(len(p)))        
        G = G/d 
        #### scaling G
        row=np.reshape(np.array(G), len(G)*len(G))
        G=G/(row.max())

        p_ecv.append('g{}'.format(g))
        for e in range(cross_list[g]):                   
            p=[(sum(Xc[f][i] for f in range(100*e, 100*e+99)) +100)/(2*100) for i in range(len(markers_index))]
            p_ecv.append(p)        

        
        ##trait1   
        P1=[]
        P1_table=[]
        P1_table_row=[]
        
        for i in range(N1):
            P1_table=[]
            for j in range(2):
                P1_table_row=[]
                for k in range(K[g]):
                    P1_table_row.append(init1[2*k+j][i])
                P1_table.append(P1_table_row)
            P1.append(P1_table)            
        H1=[sum(P1[i][j][k] for i in range(N1) for j in range(2))   for k in range(K[g])]                     
        ECV_allele_freq_trait1.append(statistics.mean(H1)/(2*N1)) #### reporting

        #trait2
        P2=[]
        P2_table=[]
        P2_table_row=[]
        
        for i in range(N2):
            P2_table=[]
            for j in range(2):
                P2_table_row=[]
                for k in range(K[g]):
                    P2_table_row.append(init2[2*k+j][i])
                P2_table.append(P2_table_row)
            P2.append(P2_table)            
        H2=[sum(P2[i][j][k] for i in range(N2) for j in range(2))   for k in range(K[g])] 
        ECV_allele_freq_trait2.append(statistics.mean(H2)/(2*N2)) #### reporting
        
        #trait3
        P3=[]
        P3_table=[]
        P3_table_row=[]
        
        for i in range(N3):
            P3_table=[]
            for j in range(2):
                P3_table_row=[]
                for k in range(K[g]):
                    P3_table_row.append(init3[2*k+j][i])
                P3_table.append(P3_table_row)
            P3.append(P3_table)            
        H3=[sum(P3[i][j][k] for i in range(N3) for j in range(2))   for k in range(K[g])]                     
        ECV_allele_freq_trait3.append(statistics.mean(H3)/(2*N3)) #### reporting   
        
        # reading phenotypic trait values
        pheno1=[]
        pheno2=[]
        pheno3=[]        
        with open('Simulation.pox.pou','r') as file: # reading phenotypes
            data=file.read().split("\n")
            data=data[2:]
            if g==1:                
                for i in range(K[g-1]):
                    counter=0
                    for j in range(len(data[i].split(' '))):
                         if data[i].split(' ')[j] == '':
                              pass
                         else:
                             counter=counter+1                  
                             if counter == 8:
                                 pheno1.append(float(data[i].split(' ')[j]))
                             elif counter == 11:
                                 pheno2.append(float(data[i].split(' ')[j]))
                             elif counter == 14:
                                 pheno3.append(float(data[i].split(' ')[j]))                                 
                                 
                ECV_phenotype_trait1.append(statistics.mean(pheno1))   ###reporting
                ECV_phenotype_trait2.append(statistics.mean(pheno2))
                ECV_phenotype_trait3.append(statistics.mean(pheno3))                
                
                pheno1=[]
                pheno2=[]
                pheno3=[]
                for i in range(K[g-1], K[g-1]+K[g]):
                    counter=0
                    for j in range(len(data[i].split(' '))):
                         if data[i].split(' ')[j] == '':
                              pass
                         else:
                             counter=counter+1                  
                             if counter == 8:
                                 pheno1.append(float(data[i].split(' ')[j]))
                             elif counter == 11:
                                 pheno2.append(float(data[i].split(' ')[j]))
                             elif counter == 14:
                                 pheno3.append(float(data[i].split(' ')[j]))                                 
                                 
                ECV_phenotype_trait1.append(statistics.mean(pheno1))   ###reporting
                ECV_phenotype_trait2.append(statistics.mean(pheno2))
                ECV_phenotype_trait3.append(statistics.mean(pheno3))
                 
            else:
                for i in range(K[g-1], K[g-1]+K[g]):
                    counter=0
                    for j in range(len(data[i].split(' '))):
                         if data[i].split(' ')[j] == '':
                              pass
                         else:
                             counter=counter+1                  
                             if counter == 8:
                                 pheno1.append(float(data[i].split(' ')[j]))
                             elif counter == 11:
                                 pheno2.append(float(data[i].split(' ')[j]))
                             elif counter == 14:
                                 pheno3.append(float(data[i].split(' ')[j]))                                 
                                 
                ECV_phenotype_trait1.append(statistics.mean(pheno1))   ###reporting
                ECV_phenotype_trait2.append(statistics.mean(pheno2))
                ECV_phenotype_trait3.append(statistics.mean(pheno3))        

        ###################### optimization #####################        
        if g < generation_number-1 :
            conflict_solutions=[]
            Female=[]
            Male=[]
            inbreeding_list=[]
            for q in range(cross_list[g]): 
                mod = Model()                           
                t=mod.addVars(2,K[g], vtype=GRB.BINARY, name="t")            
                mod.addConstrs(sum(t[m,k] for k in range(K[g])) == 1 for m in range(2))
                
                for row in conflict_solutions:
                    mod.addConstr(t[0,row[0]]+t[1,row[1]] <= 1)
                    mod.addConstr(t[1,row[0]]+t[0,row[1]] <= 1)
                
                if g >= 2:
                    for i in range(K[g]):
                        for j in range(i, K[g]):
                            if G[i][j] >= epsilon:
                                mod.addConstr(t[0,i]+t[1,j] <= 1)
                                mod.addConstr(t[1,i]+t[0,j] <= 1)                            
                   
                    obj1=0.25*quicksum(t[m,k]*H1[k] for k in range(K[g]) for m in range(2))
                    obj2=0.25*quicksum(t[m,k]*H2[k] for k in range(K[g]) for m in range(2))
                    obj3=0.25*quicksum(t[m,k]*H3[k] for k in range(K[g]) for m in range(2)) 
                     
                    mod.setObjectiveN(obj3, 0 ,  priority=2, reltol=tolerance_1[g]) #abstol=  , reltol=0.05 
                    mod.setObjectiveN(obj1, 1 ,  priority=1, reltol=tolerance_2[g])
                    mod.setObjectiveN(obj2, 2 ,  priority=0 )  
                                    
                    mod.ModelSense = GRB.MAXIMIZE
                        
                    mod.Params.OutputFlag=0
                    mod.setParam(GRB.Param.Threads,32)                                              
                    mod.optimize()
                else:
                    obj1=0.25*quicksum(t[m,k]*H1[k] for k in range(K[g]) for m in range(2))
                    obj2=0.25*quicksum(t[m,k]*H2[k] for k in range(K[g]) for m in range(2))
                    obj3=0.25*quicksum(t[m,k]*H3[k] for k in range(K[g]) for m in range(2)) 
                    
                    mod.setObjectiveN(obj3, 0 ,  priority=3, reltol=tolerance_1[g]) #abstol=  , reltol=0.05 
                    mod.setObjectiveN(obj1, 1 ,  priority=2, reltol=tolerance_2[g])
                    mod.setObjectiveN(obj2, 2 ,  priority=1 )  
                                    
                    mod.ModelSense = GRB.MAXIMIZE
                        
                    mod.Params.OutputFlag=0
                    mod.setParam(GRB.Param.Threads,32)            
                    mod._t = t            
                    mod.Params.lazyConstraints = 1                       
                    mod.optimize(call_back)
                    
                conflict_solutions_row=[]
                lis=[]
                for i in range(2):
                    for j in range(K[g]):
                        if t[i,j].x == 1:
                            conflict_solutions_row.append(j)
                            #print("Inividual {} is selected as parent{}".format(j+1,i+1))
                            if i ==0:
                                Female.append(j+1)
                            else:
                                Male.append(j+1)
                            lis.append(j)
                conflict_solutions.append(conflict_solutions_row)
                inbreeding_list.append(G[lis[0]][lis[1]])
            
            print("parents from population {} are selected.".format(g))
            ECV_Inbreeding.append(statistics.mean(inbreeding_list))   ### Reporting
                                            
            #creating pox file
            List=[]
            List.append('Model CrossID Female Male BC TC DCFem DCMale PopFemale PopMale PopBC PopTC PopDCFem PopDCMale')
            for v in range(cross_list[g]):
                List.append('1    {}     {}    {}    {}    {}    {}    {}    1    1    1    1    1    1'.format(v+1,Female[v],Male[v],Female[v],Female[v],Female[v],Male[v]))
            #writing pox file                
            with open("Simulation.pox", "w") as output:
                for i in range(len(List)):
                    if i == range(len(List)):
                        output.writelines(str(List[i]))
                    else:
                        output.writelines(str(List[i])+'\n')
                        
            with open('Simulation.pox_001_001_001_001.pop','r') as file:                                
                with open("Simulation.pop", "w") as output:
                        output.writelines(file)                    
            ##### running simulation for next progeny
            with open("QuLinePlus.mio", "w") as output:
                output.writelines(str('Simulation.ges')+'\n')
                output.writelines(str('Simulation.pop')+'\n')
                output.writelines(str('Simulation_{}.qmp'.format(g))+'\n')
                output.writelines(str('Simulation.pox')+'\n')
                output.writelines(str('Simulation_out'))
            os.system('cmd /c "C:/Users/Pouya Ahadi/Desktop/QL.exe"')
            print("generation {} is simulated.".format(g+1))
           

## outputting summary
with open("output_ecv.csv",'a',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerow('r')
    wrt.writerow(tolerance_1)
    wrt.writerow(tolerance_2)
    wrt.writerow(ECV_Inbreeding)
    wrt.writerow('-')
    wrt.writerow(ECV_allele_freq_trait1)
    wrt.writerow(ECV_allele_freq_trait2)
    wrt.writerow(ECV_allele_freq_trait3)
    wrt.writerow('-')
    wrt.writerow(ECV_phenotype_trait1)
    wrt.writerow(ECV_phenotype_trait2)
    wrt.writerow(ECV_phenotype_trait3)
    
with open("p_ecv.csv",'a',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerows(p_ecv)



################################
##    
###
####  Phenotypic Section (without self cross)
###
##
################################
for g in range(generation_number):
    if g == 0:
        with open('Error.pop','r') as file:
            data=file.read().split("\n")
            data=data[18:]
            m = len(data)
            m=int((m+1)/3)
    
            init1=[]
            init2=[]
            init3=[]
                    
            X=[]
            Xc=[] #centered X  
            for i in range(K[g]):
                
                list1=[]
                list2=[]
               
                list1=[2-int(a)  for a in data[3*i][:number_of_genes]]
            
                list2=[2-int(a)  for a in data[(3*i)+1][:number_of_genes] ]
                    
                init1_row1=[]
                init1_row2=[]
                for j in trait1_index:
                  init1_row1.append(list1[j])
                  init1_row2.append(list2[j])     
                init1.append(init1_row1)
                init1.append(init1_row2)

                init2_row1=[]
                init2_row2=[]
                for j in trait2_index:
                  init2_row1.append(list1[j])
                  init2_row2.append(list2[j])
                init2.append(init2_row1)
                init2.append(init2_row2)
                               
                init3_row1=[]
                init3_row2=[]                
                for j in trait3_index:
                    if j in trait_1_3_index:
                      init3_row1.append(1-list1[j])
                      init3_row2.append(1-list2[j]) 
                    else:
                      init3_row1.append(list1[j])
                      init3_row2.append(list2[j])              
                init3.append(init3_row1)
                init3.append(init3_row2)
                
                #creating direct Xc matrix
                Xc_row=[]
                for j in markers_index:               
                    if list1[j]+list2[j] == 0:
                        Xc_row.append(-1)
                    elif list1[j]+list2[j] == 1:
                        Xc_row.append(0)
                    elif list1[j]+list2[j] == 2:
                        Xc_row.append(1)
                        
                Xc.append(Xc_row)
        p=[ (sum(row[i] for row in Xc)+K[0])/(2*K[0]) for i in range(len(markers_index))]
        p_pheno.append(p)        
        P=[2*a-1 for a in p]        
        Z=np.array(Xc)-np.array(P)        
        Z_T = Z.T       
        G=np.dot(Z,Z_T)        
        d=2*sum(p[k]*(1-p[k])  for k in range(len(p)))        
        G = G/d 

        p_pheno.append('g{}'.format(g))
        for e in range(cross_list[g]):                   
            p=[(sum(Xc[f][i] for f in range(100*e, 100*e+99)) +100)/(2*100) for i in range(len(markers_index))]
            p_pheno.append(p)
                
        ##trait1   
        P1=[]
        P1_table=[]
        P1_table_row=[]
        
        for i in range(N1):
            P1_table=[]
            for j in range(2):
                P1_table_row=[]
                for k in range(K[g]):
                    P1_table_row.append(init1[2*k+j][i])
                P1_table.append(P1_table_row)
            P1.append(P1_table)            
        H1=[sum(P1[i][j][k] for i in range(N1) for j in range(2))   for k in range(K[g])]                     
        Pheno_allele_freq_trait1.append(statistics.mean(H1)/(2*N1)) #### reporting

        #trait2
        P2=[]
        P2_table=[]
        P2_table_row=[]
        
        for i in range(N2):
            P2_table=[]
            for j in range(2):
                P2_table_row=[]
                for k in range(K[g]):
                    P2_table_row.append(init2[2*k+j][i])
                P2_table.append(P2_table_row)
            P2.append(P2_table)            
        H2=[sum(P2[i][j][k] for i in range(N2) for j in range(2))   for k in range(K[g])] 
        Pheno_allele_freq_trait2.append(statistics.mean(H2)/(2*N2)) #### reporting
        
        #trait3
        P3=[]
        P3_table=[]
        P3_table_row=[]
        
        for i in range(N3):
            P3_table=[]
            for j in range(2):
                P3_table_row=[]
                for k in range(K[g]):
                    P3_table_row.append(init3[2*k+j][i])
                P3_table.append(P3_table_row)
            P3.append(P3_table)            
        H3=[sum(P3[i][j][k] for i in range(N3) for j in range(2))   for k in range(K[g])]                     
        Pheno_allele_freq_trait3.append(statistics.mean(H3)/(2*N3)) #### reporting  
        
        #creating pox file
        List=[]
        List.append('Model CrossID Female Male BC TC DCFem DCMale PopFemale PopMale PopBC PopTC PopDCFem PopDCMale')
        for v in range(1):
            List.append('1    {}     {}    {}    {}    {}    {}    {}    1    1    1    1    1    1'.format(v+1,2,3,2,2,2,3))
            
        with open("Simulation.pox", "w") as output:
            for i in range(len(List)):
                if i == range(len(List)):
                    output.writelines(str(List[i]))
                else:
                    output.writelines(str(List[i])+'\n')        
        
        with open("QuLinePlus.mio", "w") as output: #dummy simulation for getting pheno1
            output.writelines(str('Simulation.ges')+'\n')
            output.writelines(str('Error.pop')+'\n')
            output.writelines(str('Simulation_dummy.qmp')+'\n')
            output.writelines(str('Simulation.pox')+'\n')
            output.writelines(str('Simulation_out'))
        os.system('cmd /c "C:/Users/Pouya Ahadi/Desktop/QL.exe"')        
        
        pheno1=[]
        pheno2=[]
        pheno3=[]
        with open('Simulation.pox.pou','r') as file: # reading phenotypes
            data=file.read().split("\n")
            data=data[2:]                
            for i in range(K[g]):
                counter=0
                for j in range(len(data[i].split(' '))):
                     if data[i].split(' ')[j] == '':
                          pass
                     else:
                         counter=counter+1                  
                         if counter == 8:
                             pheno1.append(float(data[i].split(' ')[j]))
                         elif counter == 11:
                             pheno2.append(float(data[i].split(' ')[j]))
                         elif counter == 14:
                             pheno3.append(float(data[i].split(' ')[j]))                                 
                                         
        Pheno_phenotype_trait1.append(statistics.mean(pheno1))   ### reporting
        Pheno_phenotype_trait2.append(statistics.mean(pheno2))   ### reporting
        Pheno_phenotype_trait3.append(statistics.mean(pheno3))   ### reporting
        
        ###################optimization################
        conflict_solutions=[]
        Female=[]
        Male=[]
        inbreeding_list=[]
        for q in range(cross_list[g]): 
            mod = Model()                            
            t=mod.addVars(2,K[g], vtype=GRB.BINARY, name="t")           
            mod.addConstrs(sum(t[m,k] for k in range(K[g])) == 1 for m in range(2))            
            mod.addConstrs(sum(t[m,k] for m in range(2)) <= 1 for k in range(K[g])) ##self cross constr
                       
            for row in conflict_solutions:
                mod.addConstr(t[0,row[0]]+t[1,row[1]] <= 1)
                mod.addConstr(t[1,row[0]]+t[0,row[1]] <= 1)
                     
            obj1=quicksum(t[m,k]*pheno1[k] for k in range(K[g]) for m in range(2))
            obj2=quicksum(t[m,k]*pheno2[k] for k in range(K[g]) for m in range(2))
            obj3=quicksum(t[m,k]*pheno3[k] for k in range(K[g]) for m in range(2)) 
             
            mod.setObjectiveN(obj3, 0 ,  priority=3, reltol=tolerance_1[g]) #abstol=  , reltol=0.05 
            mod.setObjectiveN(obj1, 1 ,  priority=2, reltol=tolerance_2[g])
            mod.setObjectiveN(obj2, 2 ,  priority=1 )  
                            
            mod.ModelSense = GRB.MAXIMIZE

            #mod.setParam(GRB.Param.Threads,32)                                      
            mod.Params.OutputFlag=0                        
            mod.optimize()
            
            conflict_solutions_row=[]
            lis=[]
            for i in range(2):
                for j in range(K[g]):
                    if t[i,j].x == 1:
                        conflict_solutions_row.append(j)
                        #print("Inividual {} is selected as parent{}".format(j+1,i+1))
                        if i ==0:
                            Female.append(j+1)
                        else:
                            Male.append(j+1)
                        lis.append(j)
            conflict_solutions.append(conflict_solutions_row)
            inbreeding_list.append(G[lis[0]][lis[1]])
            
        Pheno_Inbreeding.append(statistics.mean(inbreeding_list))   ### Reporting
        print("parents from population {} are selected.".format(g))

                   
        #creating pox file
        List=[]
        List.append('Model CrossID Female Male BC TC DCFem DCMale PopFemale PopMale PopBC PopTC PopDCFem PopDCMale')
        for v in range(cross_list[0]):
            List.append('1    {}     {}    {}    {}    {}    {}    {}    1    1    1    1    1    1'.format(v+1,Female[v],Male[v],Female[v],Female[v],Female[v],Male[v]))
            
        with open("Simulation.pox", "w") as output:
            for i in range(len(List)):
                if i == range(len(List)):
                    output.writelines(str(List[i]))
                else:
                    output.writelines(str(List[i])+'\n')
    
        #####
        with open("QuLinePlus.mio", "w") as output:
            output.writelines(str('Simulation.ges')+'\n')
            output.writelines(str('Error.pop')+'\n')
            output.writelines(str('Simulation_0.qmp')+'\n')
            output.writelines(str('Simulation.pox')+'\n')
            output.writelines(str('Simulation_out'))
        os.system('cmd /c "C:/Users/Pouya Ahadi/Desktop/QL.exe"')
        print("generation {} is simulated.".format(g+1))
    if g > 0:
        with open('Simulation.pox_001_001_001_001.pop','r') as file:                        
            data=file.read().split("\n")
            data=data[15:]
            m = len(data)
            m=int((m+1)/3)
    
            init1=[]
            init2=[]
            init3=[]
                    
            X=[]
            Xc=[] #centered X  
            for i in range(K[g]):                
                list1=[]
                list2=[]
               
                list1=[2-int(a)  for a in data[3*i][:number_of_genes]]
            
                list2=[2-int(a)  for a in data[(3*i)+1][:number_of_genes] ]
                    
                init1_row1=[]
                init1_row2=[]
                for j in trait1_index:
                  init1_row1.append(list1[j])
                  init1_row2.append(list2[j])     
                init1.append(init1_row1)
                init1.append(init1_row2)

                init2_row1=[]
                init2_row2=[]
                for j in trait2_index:
                  init2_row1.append(list1[j])
                  init2_row2.append(list2[j])
                init2.append(init2_row1)
                init2.append(init2_row2)
                               
                init3_row1=[]
                init3_row2=[]                
                for j in trait3_index:
                    if j in trait_1_3_index:
                      init3_row1.append(1-list1[j])
                      init3_row2.append(1-list2[j]) 
                    else:
                      init3_row1.append(list1[j])
                      init3_row2.append(list2[j])              
                init3.append(init3_row1)
                init3.append(init3_row2)
                
                #creating direct Xc matrix
                Xc_row=[]
                for j in markers_index:               
                    if list1[j]+list2[j] == 0:
                        Xc_row.append(-1)
                    elif list1[j]+list2[j] == 1:
                        Xc_row.append(0)
                    elif list1[j]+list2[j] == 2:
                        Xc_row.append(1)
                        
                Xc.append(Xc_row)
      
        Z=np.array(Xc)-np.array(P)        
        Z_T = Z.T       
        G=np.dot(Z,Z_T)        
        d=2*sum(p[k]*(1-p[k])  for k in range(len(p)))        
        G = G/d 
        ### scaling G
        row=np.reshape(np.array(G), len(G)*len(G))
        G=G/(row.max())        

        p_pheno.append('g{}'.format(g))
        for e in range(cross_list[g]):                   
            p=[(sum(Xc[f][i] for f in range(100*e, 100*e+99)) +100)/(2*100) for i in range(len(markers_index))]
            p_pheno.append(p)
                    
        ##trait1   
        P1=[]
        P1_table=[]
        P1_table_row=[]        
        for i in range(N1):
            P1_table=[]
            for j in range(2):
                P1_table_row=[]
                for k in range(K[g]):
                    P1_table_row.append(init1[2*k+j][i])
                P1_table.append(P1_table_row)
            P1.append(P1_table)            
        H1=[sum(P1[i][j][k] for i in range(N1) for j in range(2))   for k in range(K[g])]                     
        Pheno_allele_freq_trait1.append(statistics.mean(H1)/(2*N1)) #### reporting

        #trait2
        P2=[]
        P2_table=[]
        P2_table_row=[]        
        for i in range(N2):
            P2_table=[]
            for j in range(2):
                P2_table_row=[]
                for k in range(K[g]):
                    P2_table_row.append(init2[2*k+j][i])
                P2_table.append(P2_table_row)
            P2.append(P2_table)            
        H2=[sum(P2[i][j][k] for i in range(N2) for j in range(2))   for k in range(K[g])] 
        Pheno_allele_freq_trait2.append(statistics.mean(H2)/(2*N2)) #### reporting
        
        #trait3
        P3=[]
        P3_table=[]
        P3_table_row=[]        
        for i in range(N3):
            P3_table=[]
            for j in range(2):
                P3_table_row=[]
                for k in range(K[g]):
                    P3_table_row.append(init3[2*k+j][i])
                P3_table.append(P3_table_row)
            P3.append(P3_table)            
        H3=[sum(P3[i][j][k] for i in range(N3) for j in range(2))   for k in range(K[g])]                     
        Pheno_allele_freq_trait3.append(statistics.mean(H3)/(2*N3)) #### reporting  
        
        pheno1=[]
        pheno2=[]
        pheno3=[]
        with open('Simulation.pox.pou','r') as file: # reading phenotypes
            data=file.read().split("\n")
            data=data[2:]
            pheno1=[]
            for i in range(K[g-1], K[g-1]+K[g]):
                counter=0
                for j in range(len(data[i].split(' '))):
                     if data[i].split(' ')[j] == '':
                          pass
                     else:
                         counter=counter+1                  
                         if counter == 8:
                             pheno1.append(float(data[i].split(' ')[j]))
                         elif counter == 11:
                             pheno2.append(float(data[i].split(' ')[j]))
                         elif counter == 14:
                             pheno3.append(float(data[i].split(' ')[j]))                                 
                                         
        Pheno_phenotype_trait1.append(statistics.mean(pheno1))   ### reporting
        Pheno_phenotype_trait2.append(statistics.mean(pheno2))   ### reporting
        Pheno_phenotype_trait3.append(statistics.mean(pheno3))   ### reporting                                                          
                
        if g < generation_number-1:
            conflict_solutions=[]
            Female=[]
            Male=[]
            inbreeding_list=[]
            for q in range(cross_list[g]): 
                mod = Model()                            
                t=mod.addVars(2,K[g], vtype=GRB.BINARY, name="t")           
                mod.addConstrs(sum(t[m,k] for k in range(K[g])) == 1 for m in range(2))            
                mod.addConstrs(sum(t[m,k] for m in range(2)) <= 1 for k in range(K[g])) ##self cross constr
                           
                for row in conflict_solutions:
                    mod.addConstr(t[0,row[0]]+t[1,row[1]] <= 1)
                    mod.addConstr(t[1,row[0]]+t[0,row[1]] <= 1)
                         
                obj1=quicksum(t[m,k]*pheno1[k] for k in range(K[g]) for m in range(2))
                obj2=quicksum(t[m,k]*pheno2[k] for k in range(K[g]) for m in range(2))
                obj3=quicksum(t[m,k]*pheno3[k] for k in range(K[g]) for m in range(2)) 
                 
                mod.setObjectiveN(obj3, 0 ,  priority=3, reltol=tolerance_1[g]) #abstol=  , reltol=0.05 
                mod.setObjectiveN(obj1, 1 ,  priority=2, reltol=tolerance_2[g])
                mod.setObjectiveN(obj2, 2 ,  priority=1 )  
                                
                mod.ModelSense = GRB.MAXIMIZE
    
                #mod.setParam(GRB.Param.Threads,32)                                      
                mod.Params.OutputFlag=0                        
                mod.optimize()
                    
                conflict_solutions_row=[]
                lis=[]
                for i in range(2):
                    for j in range(K[g]):
                        if t[i,j].x == 1:
                            conflict_solutions_row.append(j)
                            #print("Inividual {} is selected as parent{}".format(j+1,i+1))
                            if i ==0:
                                Female.append(j+1)
                            else:
                                Male.append(j+1)
                            lis.append(j)
                conflict_solutions.append(conflict_solutions_row)
                inbreeding_list.append(G[lis[0]][lis[1]])
                
            Pheno_Inbreeding.append(statistics.mean(inbreeding_list))   ### Reporting
            print("parents from population {} are selected.".format(g))
    
            #creating pox file
            List=[]
            List.append('Model CrossID Female Male BC TC DCFem DCMale PopFemale PopMale PopBC PopTC PopDCFem PopDCMale')
            for v in range(cross_list[g]):
                List.append('1    {}     {}    {}    {}    {}    {}    {}    1    1    1    1    1    1'.format(v+1,Female[v],Male[v],Female[v],Female[v],Female[v],Male[v]))
                
            with open("Simulation.pox", "w") as output:
                for i in range(len(List)):
                    if i == range(len(List)):
                        output.writelines(str(List[i]))
                    else:
                        output.writelines(str(List[i])+'\n')
                        
            with open('Simulation.pox_001_001_001_001.pop','r') as file:                                
                with open("Simulation.pop", "w") as output:
                        output.writelines(file)                    
            #####
            with open("QuLinePlus.mio", "w") as output:
                output.writelines(str('Simulation.ges')+'\n')
                output.writelines(str('Simulation.pop')+'\n')
                output.writelines(str('Simulation_{}.qmp'.format(g))+'\n')
                output.writelines(str('Simulation.pox')+'\n')
                output.writelines(str('Simulation_out'))
            os.system('cmd /c "C:/Users/Pouya Ahadi/Desktop/QL.exe"')       
            print("generation {} is simulated.".format(g+1))

## outputting summary
with open("output_pheno.csv",'a',newline='') as fp:
    wrt=csv.writer(fp,delimiter=',')
    wrt.writerow('r')
    wrt.writerow(tolerance_1)
    wrt.writerow(tolerance_2)
    wrt.writerow(Pheno_Inbreeding)
    wrt.writerow('-')
    wrt.writerow(Pheno_allele_freq_trait1)
    wrt.writerow(Pheno_allele_freq_trait2)
    wrt.writerow(Pheno_allele_freq_trait3)
    wrt.writerow('-')
    wrt.writerow(Pheno_phenotype_trait1)
    wrt.writerow(Pheno_phenotype_trait2)
    wrt.writerow(Pheno_phenotype_trait3)
    
print("run time is: {}".format(datetime.datetime.now() - begin_time)) 

