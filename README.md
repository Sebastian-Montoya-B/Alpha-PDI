# Alpha-PDI

Supplement to Montoya-Bustamante et al. *in prep.*

[Ecological Synthesis Lab](https://marcomellolab.wordpress.com) (SintECO).

Authors: Sebastián Montoya-Bustamante, Carsten F. Dormann, Boris R. Krasnov & Marco A. R. Mello.

E-mail: [s.montoyabustamante\@gmail.com](mailto:s.montoyabustamante@gmail.com).

Published originally on june 24th, 2022 (English version).

Run in R version 4.1.2 (2021-01-11) -- "Bird Hippie".

Disclaimer: You may use this script freely for commercial or non-commercial purposes at your own risk. We assume no responsibility or liability for the use of this software, convey no license or title under any patent, copyright, or mask work right to the product. We reserve the right to make changes in the software without notification. We also make no representation or warranty that such application will be suitable for the specified use without further testing or modification. If this script helps you produce any academic work (paper, book, chapter, dissertation, thesis, monograph, report, lecture, talk, etc.), please acknowledge the authors and cite this repo and the respective publication.

## List of files

See further info in the respective sections.

1. AlphaPDIfun.R -> script of the `alpha_PDI` function for calculating the alpha PDI index of generality.
2. genfun.R -> script of the `genfun` function for calculating traditional generality indices.

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
