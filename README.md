# nCOP - Network Coverage of Patients

This is my version of an algorithm called "nCOP", the original version can be access by this link: https://github.com/Singh-Lab/nCOP.git

This version was create in my scientific initiation (09/21 - 09/22), with my professor Adenilso. The objective of this is otimizate the original algorithm and try add other functions can help the users. The resume of how this algorithm works can be access in the original link: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5997485/. Or my own explanation about how this algorithm works and which modification I implementate. All credits behind this algorithm are Borislav H. Hristov and Mona Singh. The perfil of both can be access through the link of original algorithm.

################################ \
How execute this algorithm\
################################

The pass by pass can be viewer in the original code, but how I modificate a little, I will explain again:

I. Input
There are two required inputs: 
1) a network file 
2) a mutational file containing genes with list of individuals that have variants in these genes. 

Additionally, the user may provide:\
3) a weight file specifying a weight for each node in the network;\
4) a value for alpha (in which case the program skips the step selecting alpha and uses the user specified value)\
5) output prefix which is used in the beginning of the name of the output file

If you want provide a value of alpha, but not a weight file, just write a "None" in the third parameter. The same logic can applied by other parameters.

II. Output

output_prefix_results.txt is written in the Outputs directory. The file contain a list of candidate genes ranked by how frequently they appear in the randomized runs. If you pass a maf file to algorithm, other output file is created, this file content the file of normal input of nCOP, is just the maf file but only with fundamental information to algorithm. Can be pass directly in next executions.

III. How to run

1. To run with basic inputs:\
  python3 network_file.txt mutational_file.txt

2. If you want to specify any additional parameter add in orden, pass "None" if you not can any especifically param:\
  python3 network_file.txt mutational_file.txt None None My_output_file (this example don't use weight file neither alpha value, but pass a name of output file)
  
IV. Input File Formats

1. Network file: each line specifies an edge, white space delimited:\
GENE_ID GENE_ID

2. Mutational file: each line is white space delimited:\
GENE_ID INDIVIDUAL_1 INDIVIDUAL_3\
GENE_ID INDIVIDUAL_5

3. Weights file:\
GENE_ID WEIGHT\
TP53 0.59
