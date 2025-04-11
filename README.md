# Chili-Pepper-F1-Fruit-Shape-Prediction

Prediction of fruit shapes in F1 progenies based on parental image data using elliptic Fourier analysis: an initial case study with chili peppers (Capsicum annuum)

Fumiya Kondoa, Yui Kumanomido, Mariasilvia D’Andrea, Valentino Palombo, Nahed Ahmed, Shino Futatsuyama, Kazuhiro Nemoto, Kenichi Matsushima
(20205-04-11 Computer and Electronics in Agriculture Elsevier, Under review)


# NOTE
・This repositry shares all datasets and Rscripts to perform phenomic and genomic prediction methods to estimate fruit shapes of F1 progenies in chili peppers (Capsicum annuum).

・In these prediction, Elliptic Fourier descriptors (EFDs) extracted fruit image data are used and fruit contours of F1 progenies will be predicted based on parental EFDs for phenomic prediction with (phenomic prediction) or without (for genomic prediction) genomic data. 

・Using shared datasets and R scripts, four kinds of predictions for F1 progenies can pe performed as below.

1. PPmid: Prediction methods based on parental EFDs data with constant values (constant ratio between dominance and additive effects), considering dominance effects observed in F1 accession.

2. PPδ: Phenomic prediction based solely on parental EFDs data, non-considering dominance effects observed in F1 accession.

3. GP[20]: Genomic prediction based on parental EFDs and genotypic data, almost non-considering dominance effects observed in F1 accession.

4. GP[132]: Genomic prediction based on EFDs and genotypic data of 132 inbred accessions including F1 parents., almost non-considering dominance effects observed in F1 accession.

・The detailed concepts were written in our papper (Kondo et al. unpublised).

・This repositry consists four directories to perform F1 fruit shape predictions from raw EFD data. The details were described below.

#Explanations for four directories
#Part1.Averaged_EFD_Calculation
・This directory is prepared for calculating averaged EFDs of C. annuum 291 accessions (132 inbred accession + 159 F1 accessions derived from 20 parents in the 132 inbred accessions). For later tasks, these averaged EFDs will be used.

・Dataset directory in this directory contain two datasets below.

1. Accession_list.csv: Lists of plant materials in our study

2. Raw_EFD_data.csv: Raw EFD data of the 291 accessions. In our study, EFDs were obtained from a total of 8,730 images (2 angles(Direction1 and Direction2) X 291 accesssions X 3 cultivation year X 5 fruits) 
