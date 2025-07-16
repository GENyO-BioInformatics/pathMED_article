# pathMED_article

Code to reproduce the results from the article describing the pathMED R package.

PathMED may be installed from the current Bioconductor release:

``` R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("pathMED")
```

Other necessary dependencies are listed at the beginning of each script.

To reproduce the first use case (multi-omics prediction of lung adenocarcinoma), it is necessary to download the transcriptomics, proteomics and metadata from <https://kb.linkedomics.org/download#LUAD> and save into the *LinkedOmicsKB* folder. Then, the script may be executed:

``` Bash
Rscript useCase1_LUAD.R
```

The output figures and tables of this use case will be saved into the *figures* and *SupplTables* folders respectively,

To reproduce the second use case (prediction of the treatment response in breast cancer patients), the data must be obtained cloning the github repository:

``` Bash
git clone https://github.com/cclab-brca/neoadjuvant-therapy-response-predictor 
```

Then, the corresponding script must be launched:

``` Bash
Rscript useCase2_BRCA.R
```

Figures and tables will be also be saved in the *figures* and *SupplTables* folders.
