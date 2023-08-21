# Alpha-PDI

Supplement to the paper Montoya-Bustamante et al. *under review.*

[Ecological Synthesis Lab](https://marcomellolab.wordpress.com) (SintECO).

Authors: Sebastián Montoya-Bustamante, Carsten F. Dormann, Boris R. Krasnov & Marco A. R. Mello.

E-mail: [s.montoyabustamante\@gmail.com](mailto:s.montoyabustamante@gmail.com).

Published originally on June 24th, 2022 (English version).

Run in R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics".

Disclaimer: You may use this script freely for commercial or non-commercial purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this script helps you produce any academic work (software, paper, book, chapter, dissertation, thesis, monograph, report, lecture, keynote, talk, poster etc.), please acknowledge the authors and cite this repo and the respective publication.

## List of folders and files

See further info in each section.

1.  Code (folder) -\> folder containing the scripts used to calculate the generalization indices and reproduce other analyses performed in our paper.
    -   alpha_PDI.R -\> script of the `alpha_PDI` function for calculating the **αPDI** index of generalization.
    -   genfun.R -\> script of the `genfun` function for calculating the other generalization indices used in our study.
    -   wcfun.R -\> script of the `wcfun` function for calculating the Wc index of generalization (Pierotti et al. 2017).
    -   correlations.R -\> script to reproduce the Spearman correlations used to test for a relationship between the specialization parameter and the indices of generalization. It requires the `alpha_PDI` and `genfun` functions, and the vectors1.RDS data.
    -   example.R -\> script with a detailed description on how to use `alpha_PDI`, `genfun`, and `wcfun`.
    -	  QNM.R -\> script of the quantitative niche model by Fründ et al. (2016) and functions to generate the analized theoretical vectors.
    -   Data generator.R -\> script to generate all the simulated data in our analysis by using the quantitative niche model


2.  Data (folder) -\> folder containing the empirical data used in the analyses and figures.
    -   Fleas.RDS -\> list of 74 matrices of flea-mammal interactions.
    -   resource_abundances.RDS -\> list of 74 vectors of resource abundance distribution (mammal abundance) for each matrix in Fleas.RDS.
    -   sim_data.RDS -\> Data generated with the Data generator script.
    -   sim_data2.RDS -\> Data generated with the Data generator script to be used exclusively with the `wcfun` function.


3.  Figures (folder) -\> folder containing the scripts used to reproduce the unedited versions of the figures presented in the manuscript and supplementary material.
    -   `Figure 2.R` -\> script used to reproduce Figure 2 and its associated analysis. It requires the `alpha_PDI`, `genfun`, and `wcfun` functions. Follow the sequence given in the script to plot and export the figure.
    -   `Figure 3.R` -\> script used to reproduce Figure 3 and its associated analysis. It requires the `alpha_PDI` and `wcfun` functions. Follow the sequence given in the script to plot and export the figure.
    -   `Figures 4 S9-S12.R` -\> script used to reproduce Figure 4 and Figures S9 to S12. It requires the `alpha_PDI` function. Follow the sequence given in the script to plot and export the figures.
    -   `Figures S1 S2 S3.R` -\> script used to reproduce Figure S1 to S3. It requires the It requires the `alpha_PDI`, `genfun`, and `wcfun` functions. Follow the sequence given in the script to plot and export the figures.
    -   `Figures S4 S5.R` -\> script used to reproduce Figure S4 and S5. It requires the It requires the `genfun` function. Follow the sequence given in the script to plot and export the figures.
    -   `Figure S6.R` -\> script used to reproduce Figure S6 and its associated analysis. It requires the `alpha_PDI` function. Follow the sequence given in the script to plot and export the figure.
    -   `Figure S7.R` -\> script used to reproduce Figure S7 and its associated analysis. It requires the `alpha_PDI` function. Follow the sequence given in the script to plot and export the figure.
    -   `Figure S8.R` -\> script used to reproduce Figure S8 and its associated analysis. It requires the `alpha_PDI` function. Follow the sequence given in the script to plot and export the figure.
    -   Exported (folder) -\> folder containing the exported figures.

## Functionality and origin

The R code provided in this repository was designed to quantify the degree of generalization of species in consumer-resource interactions. Our analysis can also be used to quantify the degree of generalization of nodes in any consumer-resource network, given that the main assumptions explained in our paper are met. Therefore, read our paper carefully before using the functions provided here.

## Instructions

### 1. If you want to study our analysis in detail:

1.  Open the file `AlphaPDI.Rproj`.
2.  Run the `alpha_PDI`, `genfun`, and `wcfun` functions. Experiment with them by creating hypothetical interaction matrices and resource abundance vectors, or analyze your own empirical data.
3.  Use the files `correlations.R`, `Figure 2.R`, `Figure 3.R`, `Figures 4 S9-S12`, `Figures S1 S2 S3.R`, `Figures S4 S5.R`, `Figure S6.R`, `Figure S7.R`, and `Figure S8.R` to reproduce the respective figures and analyses presented in our paper.

### 2. If you want to make a quick test of our analysis or apply it to your own data:

1.  Open the file `example.R`, which contains a **tutorial**.
2.  Follow the instructions given in the tutorial.

## (1) alpha_PDI

Computes **αPDI** or **αPDI'** for an interaction matrix (or vector) and its resource abundance vector.

### Arguments

1.  data -\> matrix or vector. The original interaction matrix (or vector) with consumers in the rows and resources in the columns.

2.  abun -\> vector. It contains the resource abundances of the columns of the interaction matrix.

3.  corrected -\> logical. If "TRUE" it calculates alpha PDI corrected by the maximum possible value given the total number of interactions of the consumer.

4.  m -\> numeric. Symmetry parameter. Default is 1.

## (2) genfun

Computes the other indices of generalization calculated in our study for an interaction matrix (or vector) and its resource abundance vector.

### Arguments

1.  data -\> matrix or vector. The original interaction matrix (or vector) with consumers in rows and resources in columns.

2.  abun -\> vector. It contains the resource abundances of the columns of the interaction matrix.

## (3) wcfun

Computes the *Wc* index of generalization proposed by Pierotti et al. (2017) for an interaction matrix (or vector) and its resource abundance vector.

### Arguments

1.  data -\> matrix. The original interaction matrix with consumers in the rows and resources in the columns.

2.  abun -\> vector. It contains the resource abundances of the columns of the interaction matrix.

## Acknowledgements

We thank Baltazar González, Cristina A. Kita, Diego P. Vázquez, Francisco A. Rodrigues, Guillermo Flórez-Montero, José C. Motta Jr., Mariana Bender, Natalya Zapata-Mesa, Nico Blüthgen, Paula Lemos, Paulo R. Guimarães Jr., and Tiago B. Quental for the exciting discussions about generalization indices that inspired us to carry out this study. Special thanks go to Jochen Fründ for his recommendations on how to use his quantitative niche model, and Daniela Arenas for helping us with mammalian taxonomy. SMB thanks Ministerio de Ciencia, Tecnología e Innovación de Colombia (MinCiencias, Convocatoria Doctorados en el Exterior 860) and Coordination for the Improvement of Higher Education Personnel (CAPES, 88887.388097/2019-00) for the doctoral scholarships. MARM was funded by the Alexander von Humboldt Foundation (AvH, 1134644), National Council for Scientific and Technological Development (CNPq, 304498/2019-0), São Paulo Research Foundation (FAPESP, 2018/20695-7 and 2023/02881-6), and Dean of Research of the University of São Paulo (PRP-USP, 18.1.660.41.7). We also thank the Stack Overflow community (https://stackoverflow.com/), where we solve most of our coding dilemmas.


## Feedback

If you have any questions, corrections, or suggestions, please feel free to open an [issue](https://github.com/Sebastian-Montoya-B/Alpha-PDI/issues) or make a [pull request](https://github.com/Sebastian-Montoya-B/Alpha-PDI/pulls).


## Reference

-   Montoya-Bustamante S., Dormann C. F., Krasnov B. R., Mello M. A. R. A new index to estimate ecological generalization in consumer-resource interactions. *In prep*.


## Source repos

[Quantitative niche model](https://github.com/JochenFruend/Defaunation_SeedDispersalWebs)


## Source studies

-   Blüthgen, N., Menzel, F., & Blüthgen, N. (2006). Measuring specialization in species interaction networks. BMC Ecology, 6, 9. <https://doi.org/10.1186/1472-6785-6-9>
-   Feinsinger, P., Spears, E., & Poole, R. (1981). A Simple Measure of Niche. Ecology, 62(1), 27–32.
-   Felix, G. M., Pinheiro, R. B. P., Poulin, R., Krasnov, B. R., & Mello, M. A. R. (2022). The compound topology of host–parasite networks is explained by the integrative hypothesis of specialization. Oikos, 2022(1). <https://doi.org/10.1111/oik.08462%3C/div>
-   Fort, H., Vázquez, D. P., & Lan, B. L. (2016). Abundance and generalisation in mutualistic networks: Solving the chicken-and-egg dilemma. Ecology Letters, 19(1), 4–11. <https://doi.org/10.1111/ele.12535>
-   Fründ, J., Mccann, K. S., & Williams, N. M. (2016). Sampling bias is a challenge for quantifying specialization and network structure: lessons from a quantitative niche model. Oikos, 502–513. <https://doi.org/10.1111/oik.02256>
-   Hurlbert, S. (1978). The Measurement of Niche Overlap and Some Relatives. Ecology, 59(1), 67–77. <https://www.jstor.org/stable/1936632>
-   Manly, B. F. J., McDondald, L. L., Thomas, D. L., McDonald, T. L., & Erickson, W. P. (2002). Resource Selection by Animals (Second). Springer Netherlands. <https://doi.org/10.1007/0-306-48151-0>
-   Petraitis, P. S. (1979). Likelihood Measures of Niche Breadth and Overlap. Ecology, 60(4), 703–710.
-   Pierotti, M. E. R., Martín-Fernández, J. A., & Barceló-Vidal, C. (2017). The peril of proportions: robust niche indices for categorical data. Methods in Ecology and Evolution, 8(2), 223–231. <https://doi.org/10.1111/2041-210X.12656>
-   Pinheiro, R. B. P., Felix, G. M. F., Dormann, C. F., & Mello, M. A. R. (2019). A new model explaining the origin of different topologies in interaction networks. Ecology,100(9), e02796. <https://doi.org/10.1002/ecy.2796%3C/div>
-   Poisot, T., Canard, E., Mouquet, N., & Hochberg, M. E. (2012). A comparative study of ecological specialization estimators. Methods in Ecology and Evolution, 3(3), 537–544. <https://doi.org/10.1111/j.2041-210X.2011.00174.x>
-   Schoener, T. W. (1974). Some Methods for Calculating Competition Coefficients from Resource-Utilization Spectra. The American Naturalist, 108(961), 332–340. <https://doi.org/10.1086/282911>
-   Smith, E. P. (1982). Niche breadth, resource availability, and inference. Ecology, 63(6), 1675–1681. <https://doi.org/10.2307/1940109>
