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

## Functionality and origin

R code provided in this repository can be used to quantify generality of species in consumer-resource interactions.

## Instructions

1.

## (1) AlphaPDIfun

Computes alpha PDI for an interaction matrix (or vectior) and its resource abundance vector.

### Arguments

1.  data -> matrix. The original interaction matrix with consumers in rows and resources in columns.

2.  abun -> vector. It contains the resource abundances of the columns of the interaction matrix.
