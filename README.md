# Conditional deletion of HIF-1a provides new insight regarding the murine response to gastrointestinal infection with Salmonella Typhimurium

## Abstract
The hypoxia-inducible transcription factor 1 (HIF-1) has been shown to ameliorate different bacterial infections through enhancement of microbial killing. While the impact of HIF-1 on inflammatory diseases of the gut has been studied intensively, its function in bacterial infections of the intestine remains largely elusive. With the help of a publicly available gene expression data set, we could infer significant activation of the HIF-1 transcription factor after oral infection of mice with Salmonella Typhimurium. This prompted us to apply lineage-restricted deletion of the Hif1a locus in mice to examine cell type-specific functions of HIF-1 in this model. We show hypoxia-independent induction of HIF-1 activity upon Salmonella infection in the intestinal epithelium as well as in macrophages. Surprisingly, Hif1a deletion in intestinal epithelial cells impacted neither disease outcome nor inflammatory activity. The conditional knockout of Hif1a in myeloid cells enhanced the mRNA expression of the largely pro-inflammatory chemokine Cxcl2, revealing a potentially inflammatory effect of HIF-1 deficiency in myeloid cells in the gut in vivo. Again, the disease outcome was not affected. In vitro HIF-1-deficient macrophages showed an overall impaired transcription of pro- inflammatory factors, however, Salmonella bypassed direct intracellular, bactericidal HIF-1-dependent mechanisms in a Salmonella pathogenicity island (SPI)-2 independent manner. Taken together, our data suggest that HIF-1 in intestinal epithelial and myeloid cells is either dispensable or compensable in the immune defense against Salmonella Typhimurium.

## About
This repository contains exclusively the [script](https://github.com/saezlab/HIF1A-activity-salmonella/blob/master/hif1a_salmonella_notebook.Rmd) for the transcriptome analysis of salmonella infected mice. The analysis comprises, differential gene expression analysis, inference of TF activities with special focus on the TF HIF1A/Hif1a (using [dorothea](http://saezlab.github.io/dorothea/)) and pathway analysis (using [progeny](http://saezlab.github.io/progeny/)).

## Results
A rendered version of the markdown document can be downloaded [here](https://github.com/saezlab/HIF1A-activity-salmonella/raw/master/hif1a_salmonella_notebook.html).

## Data access
The salmonella dataset was generated by third parties and is publicly available on GEO. It can be accessed via the accession ID [GSE19174](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19174). This step is included in the [analysis script](https://github.com/saezlab/HIF1A-activity-salmonella/blob/master/hif1a_salmonella_notebook.Rmd).

## How to cite
> Robrahn L, Dupont A, Jumpertz S, Zhang K, Holland CH, Guillaume J, Rappold S, Cerovic V, Saez-Rodriguez J, Hornef MW, Cramer T. "Conditional deletion of HIF-1a provides new insight regarding the murine response to gastrointestinal infection with Salmonella Typhimurium." _bioRxiv._ 2021. DOI: [10.1101/2021.01.16.426940](https://doi.org/10.1101/2021.01.16.426940).
