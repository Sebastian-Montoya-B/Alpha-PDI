# Alpha-PDI

Supplement to the paper Montoya-Bustamante et al. *in prep.*

[Ecological Synthesis Lab](https://marcomellolab.wordpress.com) (SintECO).

Authors: Sebastián Montoya-Bustamante, Carsten F. Dormann, Boris R. Krasnov & Marco A. R. Mello.

E-mail: [s.montoyabustamante\@gmail.com](mailto:s.montoyabustamante@gmail.com).

Published originally on June 24th, 2022 (English version).

Run in R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics".

Disclaimer: You may use this script freely for commercial or non-commercial purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this script helps you produce any academic work (paper, book, chapter, dissertation, thesis, monograph, report, lecture, talk, etc.), please acknowledge the authors and cite this repo and the respective publication.

## List of folders and files

See further info in the respective sections.

1.  Code (folder) -\> folder containing the scripts to calculate the indices of generality and other analyses used in our study.
    -   alpha_PDI.R -\> script of the `alpha_PDI` function for calculating the alpha PDI index of generality.
    -   genfun.R -\> script of the `genfun` function for calculating the other generality indices used in our study.
    -   wcfun.R -\> script of the `wcfun` function for calculating the Pierotti et al (2017) index of generality.
    -   Generality_correlations.R -\> script to reproduce the Spearman correlations used to test for a relationship between the specialization parameter and the indices of generality. It requires the `alpha_PDI` and `genfun` functions, and the vectors1.RDS data.
2.  Data (folder) -\> folder containing the theoretical data used in the analyses and figures. Data were generated using the "quantitative niche model" of Fründ et al. (2016).
    -   vectors1.RDS -\> list of vectors of resource use for 32 consumers and 51 potential resources. Vectors were generated with even resource abundance distributions, the specialization parameter varying from low (0.5) to high (60), and trait matching varying gradually from 0 to 1.
    -   matrices1.RDS -\> list of vectors of resource use for 7,200 consumers. Vectors were generated using a variable number of potential resources (5, 10, and 50), with the specialization parameter varying from low (0.5) to high (60), uneven resource abundance distributions, and random trait matching. Vectors are organized in matrices of five consumers.
    -   samp_vectors_even.RDS -\> list of vectors of resource use for 3,600 consumers. Vectors were generated using a variable number of potential resources (5, 10, and 50), with the specialization parameter varying from low (0.5) to high (60), even resource abundance distributions, and random trait matching.
    -   samp_vectors_uneven.RDS -\> list of vectors of resource use for 3,600 consumers. Vectors were generated using a variable number of potential resources (5, 10, and 50), with the specialization parameter varying from low (0.5) to high (60), uneven resource abundance distributions, and random trait matching.
    -   vectors432.RDS -\> list of vectors of resource use for 432 consumers. Vectors were generated using using a variable number of potential resources (5, 10, and 50), with the specialization parameter varying from low (0.5) to high (60), even resource abundance distributions, and random trait matching.
3.  Figures (folder) -\> folder containing the scripts to reproduce the unedited versions of the figures in the manuscript and supplementary material.
    -   Figure2.R -\> script to reproduce the unedited version of Figure 2. It requires the `alpha_PDI`, `genfun`, and `wcfun` functions. Follow the sequence given in the script to create and export the figure.
    -   Figure3.R -\> script to reproduce the unedited version of Figure 3, and its associated analysis. It requires the `alpha_PDI` function. Follow the sequence given in the script to create and export the figure.
    -   FigureSup.R -\> script to reproduce Figures S1 to S36 of the supplementary material. It requires the `alpha_PDI` function. Follow the sequence given in the script to create and export the figure.

## Functionality and origin

R code provided in this repository was designed to quantify the generality of species in consumer-resource interactions. It can also be used to quantify the generality of nodes in any consumer-resource network, given that the general assumptions explained in our study are met.

## Instructions

1.  Open the `AlphaPDI.Rproj` file.
2.  Run the `alpha_PDI`, `genfun`, and `wcfun` functions. Experiment with them creating hypothetical interaction matrices and resource abundance vectors, or analyze your own empirical data.
3.  Use the scripts Figure2.R, Figure3.R, and FigureSup.R to reproduce the figures and analyses of Montoya-Bustamante et al. *in prep.*

## (1) alpha_PDI

Computes *alpha PDI* for an interaction matrix (or vectior) and its resource abundance vector.

### Arguments

1.  data -\> matrix or vector. The original interaction matrix (or vector) with consumers in the rows and resources in the columns.

2.  abun -\> vector. It contains the resource abundances of the columns of the interaction matrix.

## (2) genfun

Computes the other generality indices calculated in our study for an interaction matrix (or vector) and its resource abundance vector.

### Arguments

1.  data -\> matrix or vector. The original interaction matrix (or vector) with consumers in rows and resources in columns.

2.  abun -\> vector. It contains the resource abundances of the columns of the interaction matrix.

## (3) wcfun

Computes the *Wc* index of generality proposed by Pierotti et al. (2017) for an interaction matrix (or vector) and its resource abundance vector.

### Arguments

1.  data -\> matrix. The original interaction matrix with consumers in rows and resources in columns.

2.  abun -\> vector. It contains the resource abundances of the columns of the interaction matrix.

## Acknowledgements

We thank Baltazar González, Cristina A. Kita, Diego P. Vázquez, Francisco A. Rodrigues, Guillermo Flórez-Montero, José C. Motta Jr., Natalya Zapata-Mesa, Nico Blüthgen, Paulo R. Guimarães Jr., and Tiago B. Quental for the exciting discussions about ecological networks and niche indices that inspired us to carry out this study. SMB thanks Ministerio de Ciencia, Tecnología e Innovación de Colombia (MinCiencias, Convocatoria Doctorados en el Exterior, convocatoria 860) and Coordination for the Improvement of Higher Education Personnel (CAPES, 88887.388097/2019-00) for the doctoral scholarships. MARM was funded by the Alexander von Humboldt Foundation (AvH, 3.4-8151/15037 and 3.2-BRA/1134644), National Council for Scientific and Technological Development (CNPq, 304498/2019-0), São Paulo Research Foundation (FAPESP, 2018/20695-7), and Dean of Research of the University of São Paulo (PRP-USP, 18.1.660.41.7). We also thank the [Stack Overflow](https://stackoverflow.com) community, where we solve most of our coding dilemmas.

## Reference

-   Montoya-Bustamante S., Dormann C. F., Krasnov B. R., Mello M. A. R. In prep. A new index to estimate ecological generalization in consumer-resource interactions. *In prep*.

## Source repos

[ihsmodel](https://github.com/pinheirorbp/ihsmodel)

[Restricted-Null-Model](https://github.com/gabrielmfelix/Restricted-Null-Model)

[nestedness](https://github.com/pinheirorbp/nestedness)

## Source studies

-   Blüthgen, N., Menzel, F., & Blüthgen, N. (2006). Measuring specialization in species interaction networks. BMC Ecology, 6, 9. <https://doi.org/10.1186/1472-6785-6-9>
-   Feinsinger, P., Spears, E., & Poole, R. (1981). A Simple Measure of Niche. Ecology, 62(1), 27–32.
-   Fort, H., Vázquez, D. P., & Lan, B. L. (2016). Abundance and generalisation in mutualistic networks: Solving the chicken-and-egg dilemma. Ecology Letters, 19(1), 4–11. <https://doi.org/10.1111/ele.12535>
-   Fründ, J., Mccann, K. S., & Williams, N. M. (2016). Sampling bias is a challenge for quantifying specialization and network structure : lessons from a quantitative niche model. Oikos, 502–513. <https://doi.org/10.1111/oik.02256>
-   Hurlbert, S. (1978). The Measurement of Niche Overlap and Some Relatives. Ecology, 59(1), 67–77. <https://www.jstor.org/stable/1936632>
-   Manly, B. F. J., McDondald, L. L., Thomas, D. L., McDonald, T. L., & Erickson, W. P. (2002). Resource Selection by Animals (Second). Springer Netherlands. <https://doi.org/10.1007/0-306-48151-0>
-   Petraitis, P. S. (1979). Likelihood Measures of Niche Breadth and Overlap. Ecology, 60(4), 703–710.
-   Pierotti, M. E. R., Martín-Fernández, J. A., & Barceló-Vidal, C. (2017). The peril of proportions: robust niche indices for categorical data. Methods in Ecology and Evolution, 8(2), 223–231. <https://doi.org/10.1111/2041-210X.12656>
-   Poisot, T., Canard, E., Mouquet, N., & Hochberg, M. E. (2012). A comparative study of ecological specialization estimators. Methods in Ecology and Evolution, 3(3), 537–544. <https://doi.org/10.1111/j.2041-210X.2011.00174.x>
-   Schoener, T. W. (1974). Some Methods for Calculating Competition Coefficients from Resource-Utilization Spectra. The American Naturalist, 108(961), 332–340. <https://doi.org/10.1086/282911>
-   Smith, E. P. (1982). Niche breadth, resource availability, and inference. Ecology, 63(6), 1675–1681. <https://doi.org/10.2307/1940109>
