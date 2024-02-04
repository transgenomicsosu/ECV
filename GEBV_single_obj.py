import subprocess
import os
import csv
from gurobipy import *
import numpy as np 
import gurobipy as gp
import statistics
############################


######### Subprocess for R
# Define command and arguments
command ='Rscript'
path2script ="F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/codes/GEBV.R" ## path to file GEVB.R 
# Build subprocess command
cmd = [command, path2script]
# check_output will run the command and store to result


############################
n_traits=3
epsilon=0.25
generation_number=4 ### add one more to the generations you desire
cross_list=[50,10,5,5,5] ### last component is dummy! mot important
K=[10000, 5000, 1000, 500, 500, 500]  ### last component is dummy! mot important

tolerance_1=[0.17, 0.05 , 0.05 , 0.05]# tolerance for first objective
tolerance_2=[0   ,  0   , 0    , 0.05]# tolerance for second objective

###DEfining call back
def call_back(model, where):
    if where == GRB.Callback.MIPSOL:
        # make a list of edges selected in the solution
        vals = model.cbGetSolution(model._t)
        indices = gp.tuplelist((m,k)  for m, k in model._t.keys() if vals[m, k] == 1)
        if G[indices[0][1]][indices[1][1]] >= epsilon:
            model.cbLazy(gp.quicksum(model._t[m, k]  for m, k in indices) <= 1)
################################
            
pop_num = 30  ## number of population we simulate

"""
We store input files for each population in a different folder. popi folder corresponds to folder for population i

In each folder we require following files:
    1. Error.pop
    2. markers_index.csv: index of genes corresponding to markers
    3. trait1_index.csv: index of genes corresponding to trait 1
    2. trait2_index.csv: index of genes corresponding to trait 2
    2. trait3_index.csv: index of genes corresponding to trait 3  
"""

for ww in  range(pop_num):
    g=0
    print("rep: {}".format(ww))    
    os.chdir("F:/Course Files/Ms.C/Rsearch/Final Simulations & Files/MO_GEBV_SO/pop{}".format(ww)) ## path to folder of each population

    markers_index=[]
    trait1_index=[]
    trait2_index=[]
    trait3_index=[]
    
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
    
    
    
    GEBV_allele_freq_trait1=[]
    GEBV_phenotype_trait1=[]
    GEBV_allele_freq_trait2=[]
    GEBV_phenotype_trait2=[]
    GEBV_allele_freq_trait3=[]
    GEBV_phenotype_trait3=[]
    
    GEBV_Inbreeding=[]
    
    
    p=[] #vector of frequencies
    P=[] #matrix of frequencies
    
    p_GEBV=[]
    
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
           
        for i in range(K[0]):      
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
                              
            Xc_row=[]
            for j in markers_index:               
                if list1[j]+list2[j] == 0:
                    Xc_row.append(-1)
                elif list1[j]+list2[j] == 1:
                    Xc_row.append(0)
                elif list1[j]+list2[j] == 2:
                    Xc_row.append(1)               
            Xc.append(Xc_row)


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
        GEBV_allele_freq_trait1.append(statistics.mean(H1)/(2*N1)) #### reporting
    
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
        GEBV_allele_freq_trait2.append(statistics.mean(H2)/(2*N2)) #### reporting
    
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
        GEBV_allele_freq_trait3.append(statistics.mean(H3)/(2*N3)) #### reporting  

                          
    p=[ (sum(row[i] for row in Xc)+K[0])/(2*K[0]) for i in range(len(markers_index))]        
    P=[2*a-1 for a in p] 
    
    
    Z=np.array(Xc)-np.array(P)        
    Z_T = Z.T       
    G=np.dot(Z,Z_T)        
    d=2*sum(p[k]*(1-p[k])  for k in range(len(p)))        
    G = G/d 

    p_GEBV.append('g{}'.format(g))
    for e in range(cross_list[g]):                   
        p=[(sum(Xc[f][i] for f in range(100*e, 100*e+99)) +100)/(2*100) for i in range(len(markers_index))]
        p_GEBV.append(p)
        
        
    #creating  dummy pox file
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
                                     
    GEBV_phenotype_trait1.append(statistics.mean(pheno1))   ### reporting
    GEBV_phenotype_trait2.append(statistics.mean(pheno2))   ### reporting
    GEBV_phenotype_trait3.append(statistics.mean(pheno3))   ### reporting
        
        
    GEBV1=[]
    with open("GEBV1_gen0.csv") as file:
        data=file.read().split("\n")
        for rr in range(K[g]):
            GEBV1.append(float(data[rr]))


    ###################optimization################
    conflict_solutions=[]
    Female=[]
    Male=[]
    inbreeding_list=[]
    for q in range(cross_list[g]): 
        mod = Model()
                        
        t=mod.addVars(2,K[0], vtype=GRB.BINARY, name="t")
        
        mod.addConstrs(sum(t[m,k] for k in range(K[0])) == 1 for m in range(2))
        
        mod.addConstrs(sum(t[m,k] for m in range(2)) <= 1 for k in range(K[0])) ##self cross constr
                   
        for row in conflict_solutions:
            mod.addConstr(t[0,row[0]]+t[1,row[1]] <= 1)
            mod.addConstr(t[1,row[0]]+t[0,row[1]] <= 1)
                 
        obj1=quicksum(t[m,k]*GEBV1[k] for k in range(K[0]) for m in range(2))

        mod.setObjective(obj1, GRB.MAXIMIZE)
        #mod.setParam(GRB.Param.Threads,32)                                      
        mod.Params.OutputFlag=0                        
        mod.optimize()
        
        conflict_solutions_row=[]
        lis=[]
        for i in range(2):
            for j in range(K[0]):
                if t[i,j].x >= 0.9:
                    conflict_solutions_row.append(j)
                    if i ==0:
                        Female.append(j+1)
                    else:
                        Male.append(j+1)
                    lis.append(j)
        conflict_solutions.append(conflict_solutions_row)
        inbreeding_list.append(G[lis[0]][lis[1]])           
    GEBV_Inbreeding.append(statistics.mean(inbreeding_list))   ### Reporting
                
    #creating random pox file
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

    with open("QuLinePlus.mio", "w") as output:
        output.writelines(str('Simulation.ges')+'\n')
        output.writelines(str('Error.pop')+'\n')
        output.writelines(str('Simulation_0.qmp')+'\n')
        output.writelines(str('Simulation.pox')+'\n')
        output.writelines(str('Simulation_out'))
    os.system('cmd /c "C:/Users/Pouya Ahadi/Desktop/QL.exe"')
    print("generation {} is simulated.".format(g+1))


    for g in range(1,5): 
        M=[]             
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
                M_row=[]
                M_row=[]
                list1=[]
                list2=[]
               
                list1=[2-int(a)  for a in data[3*i][:number_of_genes]]           
                list2=[2-int(a)  for a in data[(3*i)+1][:number_of_genes] ]

                for kk in range(number_of_genes):
                    if list1[kk]+list2[kk] == 2:
                        M_row.append(1)
                    elif list1[kk]+list2[kk] == 1:
                        M_row.append(0)
                    elif list1[kk]+list2[kk] == 0:
                        M_row.append(-1)
                M.append(M_row)
                  
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
    
        with open("M.csv".format(g),'w',newline='') as fp:
            wrt=csv.writer(fp,delimiter=',')
            wrt.writerows(M)
    
        p_GEBV.append('g{}'.format(g))
        for e in range(cross_list[g]):                   
            p=[(sum(Xc[f][i] for f in range(100*e, 100*e+99)) +100)/(2*100) for i in range(len(markers_index))]
            p_GEBV.append(p)
                
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
        GEBV_allele_freq_trait1.append(statistics.mean(H1)/(2*N1)) #### reporting
    
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
        GEBV_allele_freq_trait2.append(statistics.mean(H2)/(2*N2)) #### reporting
        
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
        GEBV_allele_freq_trait3.append(statistics.mean(H3)/(2*N3)) #### reporting 
        
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
                                         
        GEBV_phenotype_trait1.append(statistics.mean(pheno1))   ### reporting
        GEBV_phenotype_trait2.append(statistics.mean(pheno2))   ### reporting
        GEBV_phenotype_trait3.append(statistics.mean(pheno3))   ### reporting  
        ########################
        input_matrix=[]
        for tt in range(K[g]):
            input_matrix.append(pheno1[tt])
            
        with open("Input1.csv",'w',newline='') as fp:
            wrt=csv.writer(fp,delimiter=',')
            for row in input_matrix:
                wrt.writerow([row])
                
        #################
        GEBV1=[]   
        ## reading GEBV1
        Input=[]
        with open("Input1.csv","r") as f:
            file = csv.reader(f, delimiter=',')
            for row in file:
                Input.append(row)
            
        with open("Input.csv",'w',newline='') as fp:
            wrt=csv.writer(fp,delimiter=',')
            for row in Input:
                wrt.writerow(row)
            
        x = subprocess.check_output(cmd, universal_newlines=True)
        
        with open("GEBV.csv") as file:
            data=file.read().split("\n")
            for rr in range(K[g]):
                GEBV1.append(float(data[rr+1]))
        
        with open("GEBV1_gen{}.csv".format(g),'w',newline='') as fp:
            wrt=csv.writer(fp,delimiter=',')
            for i in range(len(GEBV1)):
                wrt.writerow([GEBV1[i]]) 
                
        ###################################
        if g <= 3:
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
                         
                obj1=quicksum(t[m,k]*GEBV1[k] for k in range(K[g]) for m in range(2))
                mod.setObjective(obj1) #abstol=  , reltol=0.05                                 
                mod.ModelSense = GRB.MAXIMIZE
                                                 
                mod.Params.OutputFlag=0                        
                mod.optimize()
                
                conflict_solutions_row=[]
                lis=[]
                for i in range(2):
                    for j in range(K[g]):
                        if t[i,j].x >= 0.9:
                            conflict_solutions_row.append(j)
                            if i ==0:
                                Female.append(j+1)
                            else:
                                Male.append(j+1)
                            lis.append(j)
                conflict_solutions.append(conflict_solutions_row)
                inbreeding_list.append(G[lis[0]][lis[1]])                
            GEBV_Inbreeding.append(statistics.mean(inbreeding_list))   ### Reporting            
                       
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
    with open("output_GEBV_SO.csv",'w',newline='') as fp:
        wrt=csv.writer(fp,delimiter=',')
        wrt.writerow('r')
        wrt.writerow(GEBV_Inbreeding)
        wrt.writerow('-')
        wrt.writerow(GEBV_allele_freq_trait1)
        wrt.writerow(GEBV_allele_freq_trait2)
        wrt.writerow(GEBV_allele_freq_trait3)
        wrt.writerow('-')
        wrt.writerow(GEBV_phenotype_trait1)
        wrt.writerow(GEBV_phenotype_trait2)
        wrt.writerow(GEBV_phenotype_trait3)

    print("replication {} is done!".format(ww))    
            