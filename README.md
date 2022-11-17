DECEPTICON
-
DEconvolution for CEll Proportion esTImation by CONsensus

Installation
-
Most dependent packages can be installed by the code below:
```R
#Biobase
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")

#bseqsc
install.packages(pkgs = 'devtools')
devtools::install_github('shenorrLabTRDF/csSAM')
BiocManager::install('NMF')
devtools::install_github('shenorrlab/bseq-sc')
library(bseqsc)

#DeconRNAseq
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DeconRNASeq")
library(DeconRNASeq)

#EPIC
devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

#MCPcounter
devtools::install_github("ebecht/MCPcounter",ref="master", subdir="Source")

#MuSic(The benchmark was performed using MuSic version 0. 2. 0)
BiocManager::install("TOAST")
BiocManager::install("SingleCellExperiment")
devtools::install_github('xuranw/MuSiC')

#preprocessCore
BiocManager::install("preprocessCore")

#psych
install.packages('psych')

#SCDC
install.packages("remotes")
remotes::install_github("renozao/xbioc",force = TRUE)
devtools::install_github("meichendong/SCDC")

#vegan
install.packages("vegan")

#BATMAN(Batman installation package is downloaded from website https://r-forge.r-project.org/scm/viewvc.php/pkg/?root=batman.)
install.packages("plotrix")
install.packages("doSNOW")
install.packages("batman_1.2.1.13.tar.gz", repos = NULL, type = "source")
```
Install DECEPTICON
-
Now, we can install the DECEPTICON by downloading the installation file in this page (Decepticon_1.0.0.0.tar.gz), and install it:
```R
install.packages("Decepticon_1.0.0.0.tar.gz", repos = NULL, type = "source")
```
Preparatory Work
-
Download the required files from exdata in the current interface:
1. exdata contains five expression templates provided by DECEPTICON, which are stored in a folder and named as "signature_matrix".
2. exdata also contains other files required by DECEPTICON. After downloading, these files are stored in the same path as "signature_matrix".
   Note: the file type of the expression template is txt, with rows representing genes and columns representing cell types; the names of all files cannot be changed.
3. If the user wants to use the custom expression template,it needs to be stored in the "custom_signature_matrix" folder of the same path.
4. Create a folder named "res" at the same path, and all the results will be saved in this folder.
5.Create a folder named "batman" to store BATMAN's data.
6.Before running DECEPTICON, need to modify BATMAN's data. See "batman _ usage" for detailed use of BATMAN.
As shown in the figure

Run DECEPTICON 
-
```R
library(Decepticon)
```
The input requires a bulk sample (a m * n matrix with m genes and n samples) and expression template(s) for cell type (i * k matrix with i genes and k cell types)

run_DECEPTICON(bulk.samples, RUNpath, path, light=FALSE, custom.signature=FALSE, signature.matrix=NULL)

`bulk.samples` in demo data is HNSCC dataset, which contains 23 samples and 23,687 genes.

`RUNpath` is the path where the above files (“signature_matrix”, “cibersort.R”, “DOCKER_codes.R”, etc.) are stored.

`path` is the path where all results are stored, i.e. “./RUNpath/res”

`light` set TRUE for light version of DECEPTICON, running time reduced.

`custom.signature` set TRUE for custom signature matrix.

`signature.matrix` When custom.signature is set to TRUE, need to enter expression template(s).

Example:
```R
Run_DECEPTICON(bulk.samples, RUNpath, path, custom.signature=TRUE, signature.matrix=c(“./custom.signature1.txt”, “./custom.signature2.txt”)) #custom signature
Run_DECEPTICON(bulk.samples, RUNpath, path, light = TRUE) #light version
```
