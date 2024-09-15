# MECDMC

## Method Description
MECDMC is a novel approach to predict potential SM-miRNA associations based on matrix
enhancement and collaborative double matrix completion.
 <a href="https://www.example.com" title="论文链接(待添加)">论文链接</a>

## Operating Environmention
PyCharm == 2021.2.1\
Python == 3.11.2\
Windows == 10\
Processor == Intel(R) Core(TM) i5 8300H 2.30GHz CPU with 16G RAM

## Required Packages
numpy == 1.24.1\
matplotlib == 3.6.2\
scipy == 1.10.0\
scikit-learn == 1.2.0

## File Description

### Dataset:
This folder contains three datasets (dataset 1, dataset 2, and newTest). Taking dataset 1 as an example, it includes the following files:\
**SM number.xlsx:** It contains the CID of 831 SMs and the row number of each SM in the association matrix. For example, the first row of the association matrix represents SM 'CID:137'.\
**miRNA number.xlsx:** It contains the name of 541 human-related miRNAs and the column number of each miRNA in the association matrix. For example, the first column of the association matrix represents miRNA 'hsa-let-7a-1'.\
**miRNA similarity matrix.txt:** It contains the integrated miRNA similarity (the dimension is 541 $\times$ 541), where each row and its corresponding column represent a specific miRNA. The element in the i-th row and j-th column denotes the similarity score between 
$miRNA_{i}$ and $miRNA_{j}$.\
**SM-miRNA-confirmed associations.txt:** It contains 664 SM-miRNA known associations. For example, the first row in the file, (75,177), denotes that the SMs represented in the 76-th row of the association matrix, 'CID:2578', is associated with the miRNA represented in the 178-th column of the association matrix, 'hsa-mir-200c'.\
**SM-miRNA association matrix.txt:** It is the constructed association matrix (the dimension is 831 $\times$ 541), each row of which represents a specific SM, and each column represents a specific miRNA. The (i,j)-th element of the association matrix, 
$m_{ij}$, is set to 1 if $SM_{i}$ is associated with $miRNA_{j}$, otherwise it is set to 0. The matrix has a total of 664 ones and 448907 zeros.\
### GRBF.py:
This file contains Python code that processes the integrated SM/miRNA similarity using Gaussian radial basis functions.
### MECDMC.PY:
This file contains the Python code for our algorithm. 
### Running steps:
Step 1: we run “GRBF.py”. Its input are the SM similarity matrix and the miRNA similarity matrix. The output are the final SM similarity matrix and the final miRNA similarity matrix.\
Step 2: we execute “MECDMC.py”. Its input are the SM-miRNA association matrix, the final SM similarity matrix, the final miRNA similarity matrix. Its output is the SM-miRNA predictive score matrix.
### Contact
If you have any problems or find mistakes, please feel free to contact me: z22070050@s.upc.edu.cn
