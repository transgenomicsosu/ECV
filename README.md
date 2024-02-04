# Simulating Initial Population
We use Qu-GENE proposed by Ali et al. (2020) to simulate populations. It includes two softwares Engine.exe for simulating the initial population and QuLinePlus (QL.exe) for simulating other generations. We use gene_maker.py to simulate initial populations. The output of gene_maker.py includes: (assuming 3 traits)
- gene_output which includes gene information required for the .qug input file for Engine.exe
- markers_index.csv: index of genes corresponding to markers
- trait1_index.csv: index of genes corresponding to trait 1
- trait2_index.csv: index of genes corresponding to trait 2
- trait3_index.csv: index of genes corresponding to trait 3

After making the .qug file required for Error.exe, we can simulate the initial population. The output Error.pop will include genotype information of the initial population.

# Optimization 
We have two types of optimization:

## Single Objective:
Where we optimize a single target trait.
- single_obj(ecv&pheno).py : This code returns results for two approaches: ECV selection approach and selection based on phenotypic trait value
This code reads, Error.pop file as well as markers_index.csv, trait1_index.csv, trait2_index.csv, trait3_index.csv.
By specifying all design parameters, these codes optimally choose parents and run for the desired number of cycles.
The simulation for a new generation based on selected parents is done by QL.exe file which will be called within the code repeatedly.
Output of code is allele frequency and phenotypic trait value corresponding to the target trait and other existing traits (not target).
- GEBV_single_obj.py: This code returns results for selection based on GEBV values. The general structure of the code is similar to the previous case. The only difference is that it requires to use GEBV.R file within compiling the code to generate GEBV trait values at each iteration.

Calculating GEBV values in GEBV.R is based on the rrBlup package proposed by Endelman2011.


## Multi-Objective:
where we optimize multiple objectives using the lexicographic approach. The order of importance for traits and the tolerances can be specified within the code.
- multi_obj(ecv&pheno).py: This code returns allele frequency and phenotypic trait value for all generations based on two approaches, ECV selection and selection based on phenotypic values.
- GEBV_multi_obj.py: This code returns allele frequency and phenotypic trait value for all of generations based on selection based on GEBV values.
 
The structure of codes is similar to the single trait case except that in the optimization part, we use the lexicographic approach to select based on multiple objectives.

# Results and Visualization
The results are saved in CSV files. The final step is to organize the results and visualize them.
- result_data_frame_maker.py: Given a correct path to the folder of each population, this code combines the results in a single data frame. The Output of the code includes results for single-objective cases and inbreeding of selected parents and results for multi-objective cases and inbreeding of selected parents.
- plot_maker.py: Given the correct path to the final CSV file created by result_data_frame_maker.py, this code returns corresponding plots.
