# Multi-angled transcriptome analyses reveal chameleon-like behavior of macrophages in the porcine respiratory tract
This is the R codes which were used to generate figures and conduct the downstream analysis in our manuscript, “Multi-angled transcriptome analyses reveal chameleon-like behavior of macrophages in the porcine respiratory tract”. _Under revision in Frontiers in Immunology_

Dayoung Oh, Nick Vereecke, Wim Trypsteen, Heesoo Song, Ward De Spiegelaere, Bert Devriendt, Sieglinde Coppens, Jo Vandesompele, Hans J. Nauwynck


<!-- ABOUT THE PROJECT -->
## Summary

Macrophages  are key players in the immune response against viral infections in pigs. Recent research has identified different subsets of porcine nasal macrophages with varying levels of susceptibility to different strains of porcine reproductive and respiratory syndrome virus (PRRSV). Here, total RNA sequencing was performed on nasal macrophage subsets (NaSn- and NaSn+) and lung macrophages (LuSn+), isolated by FACS or LCM  with immunofluorescence staining. The results revealed the transcriptomic profile of nasal macrophages was distinct from that of lung macrophages. Further analysis identified differentially expressed genes in the two subsets of nasal macrophages, and cell type signature analysis from LCM-RNAseq revealed microenvironment-related features, indicating their adaptation to neighboring cells and extracellular matrix and their different roles in interacting with pathogens and immune responses. Our study enhances the understanding of macrophage subsets and the nasal immunity in pigs and provides advanced strategies for targeting specific macrophage subsets critical in the immune response against viruses, improving treatment of viral infections in pigs.

<img src="https://github.com/HeesooSong/PorcineMOsubsets-FACS-LCM-RNASeq/blob/main/Source/Graphical_Abstract.png?raw=true" width=50%>


## Data Analysis & Visualization
All data analysis and visalizations can be reproduced from the scripts and data given above.

### Prerequisites

* Windows 10
* R (4.0.4)
* RStudio
* Python 3.8.5
* Spyder (Python IDE)

### Built With

* R packages
    ```{r}
    DESeq2
    dplyr
    edgeR
    extrafont
    ggplot2
    ggpubr
    ggrepel
    limma
    RColorBrewer
    remotes
    reshape2
    scales
    tibble
    tidyverse
    VennDiagram
    ```


<!-- CONTRIBUTING -->
## Author Contributions

D.O. and H.J.N. designed the project. D.O. designed and performed tissue isolation, FACS, LCM, RNA extraction, and immunofluorescence staining. B.D. supervised and performed FACS. W.D.S supervised LCM. W.T. and J.V. supervised sequencing material preparation and designed and performed library preparation and sequencing. N.V. and S.C. performed data processing. D.O. and H.S. analyzed and visualized the data. D.O wrote the manuscript with help from all the authors.


<!-- CONTACT -->
## Contact

DaYoung Oh - [LinkdIn](https://www.linkedin.com/in/dayoung-oh-6053b6132/) -  Dayoung.Oh@UGent.be    
Heesoo Song - [LinkdIn](https://www.linkedin.com/in/heesoosong/) - cleos0409@gmail.com    
Sieglinde Coppens - [LinkdIn](https://www.linkedin.com/in/sieglinde-coppens/) - sieglinde.coppens@pathosense.com


<!-- ACKNOWLEDGMENTS-->
## Acknowledgments

We thank the Center for Medical Genetics and CRIG for their technical help with the library preparations and sequencing and the PathoSense team for their assistance with data processing.


