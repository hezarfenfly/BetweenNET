# BetweenNet: Ranking Cancer Drivers via Betweenness-based Outlier Detection and Random Walks

### This is the original repository for the BetweenNet paper


**Installing the dependencies**

```
pip install -r requirements.txt
```
## **1 - Input** - Betweenness

1. The PPI network file:

The file is located at data/InfluenceMatrix

```
gene1 gene2 confidence value
MDM2  TP53  0.99
APP APP 0.99
MYC MAX 0.98
...
```
2. Patients data (Normal and Tumor):

Files are located at data/Betweenness
```
gene_id        TCGA-A7-A0CE-01A    gene_id        TCGA-A7-A0CE-11A
ACIN1|22985    2916.8574           ACIN1|22985    3377.1429
ACLY|47        4599.7118           ACLY|47        3405.7143
removed        removed             ACMSD|130013   26.3492
ACN9|57001     11.1874             ACN9|57001     176.1905
...
```
Note: "removed"  refers to a gene is not expressed or is a mutated gene

### **Run Betweenness calculation**



## **2 - Input** - Bipartite Graph
1. Mutation matrix
The file contains a matrix of mutated genes.
The file is located at data/MutationMatrix.csv

```
Genes   g1  g2  g3  g4 ... gn
p1      1   0   0   1  ... 1 
p2      0   1   1   0  ... 1 
p3      1   1   1   0  ... 0 
p4      1   0   0   1  ... 1 
....
```
2. Outliers matrix
The file contains a matrix of outlier genes.
The file is located at data/OutliersMatrix.csv

```
Genes   g1      g2      g3     ...  gn
p1      True    True    False  ...  True 
p2      True    False   False  ...  False 
p3      False   True    False  ...  False 
p4      True    False   True   ...  True 
....
```

### **Outliers matrix generation**
To generate the outliers matrix, run  generate_outliers.py script, as follow:

```
python generate_outliers.py [Input directory] [Betweenness Results Path] [Genes list]
```

## **3 - Input** - Random Walk
1 - The nodes to index file mapping:
The file is located at out/hint_index_file.txt




## **Run**

To run OLDRIM on the given input files:

```
cd src
python oldrim.py
```
`src/config.py` contains different parameters for the inout files path and the functions used for scoring

