# Alpha-PDI

Supplement to Montoya-Bustamante et al. *in prep.*

[Ecological Synthesis Lab](https://marcomellolab.wordpress.com) (SintECO).

Authors: Sebastián Montoya-Bustamante, Carsten F. Dormann, Boris R. Krasnov & Marco A. R. Mello.

E-mail: [s.montoyabustamante\@gmail.com](mailto:s.montoyabustamante@gmail.com).

Published originally on june 24th, 2022 (English version).

Run in R version 4.1.2 (2021-01-11) -- "Bird Hippie".

Disclaimer: You may use this script freely for commercial or non-commercial purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this script helps you produce any academic work (paper, book, chapter, dissertation, thesis, monograph, report, lecture, talk, etc.), please acknowledge the authors and cite this repo and the respective publication.

## List of folders and files

See further info in the respective sections.

1. Code (folder) -> Folder containing the scripts to calculate the indices of generality used in the manuscript.
   -AlphaPDIfun.R -> script of the `alpha_PDI` function for calculating the alpha PDI index of generality.
   -genfun.R -> script of the `genfun` function for calculating traditional generality indices.
   -wcfun -> script of the `wcfun` function for calculating Pierotti et al (2017) index of generality.

2. Figures (folder) -> Folder containing the scripts to reproduce the unedited versions of the figures in the manuscript and supplementary material.
2.1. Figure2.R -> script to reproduce the unedited version of Figure 2.  It requires the `alpha_PDI`, `genfun` and `wcfun` functions. Follow the sequence given in the script to create and export the figure.
2.2. Figure3.R -> script to reproduce the unedited version of Figure 3.  It requires the `alpha_PDI` function. Follow the sequence given in the script to create and export the figure.
2.3. FigureSup.R -> script to reproduce Figures S1 to S36 of the supplementary material.  It requires the `alpha_PDI` function. Follow the sequence given in the script to create and export the figure.

## Functionality and origin

R code provided in this repository can be used to quantify generality of species in consumer-resource interactions.

## Instructions

1.

## (1) AlphaPDIfun

Computes alpha PDI for an interaction matrix (or vectior) and its resource abundance vector.

### Arguments

1.  data -> matrix. The original interaction matrix with consumers in rows and resources in columns.

2.  abun -> vector. It contains the resource abundances of the columns of the interaction matrix.

## (2) genfun

Computes previousy published generality indices for an interaction matrix (or vector) and its resource abundance vector.

### Arguments

1.  data -> matrix. The original interaction matrix with consumers in rows and resources in columns.

2.  abun -> vector. It contains the resource abundances of the columns of the interaction matrix.

## Acknowledgements

We thank Baltazar González, Cristina A. Kita, Diego P. Vázquez, Francisco A. Rodrigues, Guillermo Flórez-Montero, José C. Motta Jr., Natalya Zapata-Mesa, Paulo R. Guimarães Jr., and Tiago B. Quental for the exciting discussions about ecological networks and niche indices that inspired us to carry out this study. SMB thanks Ministerio de Ciencia, Tecnología e Innovación de Colombia (MinCiencias, Convocatoria Doctorados en el Exterior, convocatoria 860) and Coordination for the Improvement of Higher Education Personnel (CAPES, 88887.388097/2019-00) for the doctoral scholarships. MARM was funded by the Alexander von Humboldt Foundation (AvH, 3.4-8151/15037 and 3.2-BRA/1134644), National Council for Scientific and Technological Development (CNPq, 304498/2019-0), São Paulo Research Foundation (FAPESP, 2018/20695-7), and Dean of Research of the University of São Paulo (PRP-USP, 18.1.660.41.7). We also thank the [Stack Overflow](https://stackoverflow.com) community, where we solve most of our coding dilemmas.

## Source studies

1.  Blüthgen, N., Menzel, F., & Blüthgen, N. (2006). Measuring specialization in species interaction networks. BMC Ecology, 6, 9. <https://doi.org/10.1186/1472-6785-6-9>
2.  Feinsinger, P., Spears, E., & Poole, R. (1981). A Simple Measure of Niche. Ecology, 62(1), 27–32.
3.  Fort, H., Vázquez, D. P., & Lan, B. L. (2016). Abundance and generalisation in mutualistic networks: Solving the chicken-and-egg dilemma. Ecology Letters, 19(1), 4–11. <https://doi.org/10.1111/ele.12535>
4.  Hurlbert, S. (1978). The Measurement of Niche Overlap and Some Relatives. Ecology, 59(1), 67–77. <https://www.jstor.org/stable/1936632>
5.  Montoya-Bustamante S., Dormann C. F., Krasnov B. R., Mello M. A. R. In prep. A new index to estimate ecological generalization in consumer-resource interactions.
6.  Petraitis, P. S. (1979). Likelihood Measures of Niche Breadth and Overlap. Ecology, 60(4), 703–710.
7.  Pierotti, M. E. R., Martín-Fernández, J. A., & Barceló-Vidal, C. (2017). The peril of proportions: robust niche indices for categorical data. Methods in Ecology and Evolution, 8(2), 223–231. <https://doi.org/10.1111/2041-210X.12656>
8. Poisot, T., Canard, E., Mouquet, N., & Hochberg, M. E. (2012). A comparative study of ecological specialization estimators. Methods in Ecology and Evolution, 3(3), 537–544. <https://doi.org/10.1111/j.2041-210X.2011.00174.x>
9. Schoener, T. W. (1974). Some Methods for Calculating Competition Coefficients from Resource-Utilization Spectra. The American Naturalist, 108(961), 332–340. <https://doi.org/10.1086/282911>
10. Smith, E. P. (1982). Niche breadth, resource availability, and inference. Ecology, 63(6), 1675–1681. <https://doi.org/10.2307/1940109>
