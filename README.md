# Chili-Pepper-F1-Fruit-Shape-Prediction

Prediction of fruit shapes in F1 progenies based on parental image data using elliptic Fourier analysis: an initial case study with chili peppers (Capsicum annuum)

Fumiya Kondoa, Yui Kumanomido, Mariasilvia D’Andrea, Valentino Palombo, Nahed Ahmed, Shino Futatsuyama, Kazuhiro Nemoto, Kenichi Matsushima
(20205-04-11 Computer and Electronics in Agriculture Elsevier, Under review)


# NOTE
・This repositry shares all datasets and Rscripts to perform phenomic and genomic prediction methods to estimate fruit shapes of F1 progenies in chili peppers (Capsicum annuum).

・In these prediction, Elliptic Fourier descriptors (EFDs) extracted fruit image data are used and fruit contours of F1 progenies will be predicted based on parental EFDs for phenomic prediction with (phenomic prediction) or without (for genomic prediction) genomic data. 

・Using shared datasets and R scripts, four kinds of predictions for F1 progenies can pe performed as below.

1. GP[20]: Genomic prediction based on parental EFDs and genotypic data, almost non-considering dominance effects observed in F1 accession.

2. GP[132]: Genomic prediction based on EFDs and genotypic data of 132 inbred accessions including F1 parents., almost non-considering dominance effects observed in F1 accession.

3. PPmid: Prediction methods based on parental EFDs data with constant values (constant ratio between dominance and additive effects), considering dominance effects observed in F1 accession.

4. PPδ: Phenomic prediction based solely on parental EFDs data, non-considering dominance effects observed in F1 accession.

・The detailed concepts were written in our papper (Kondo et al. unpublised).

・This repositry consists four directories to perform F1 fruit shape predictions from raw EFD data. We can perform each tasks by downloading each directory. The details were described below.


# Explanations for four directories

# Part1.Averaged_EFD_Calculation

・This directory is prepared for calculating averaged EFDs of C. annuum 291 accessions (132 inbred accession + 159 F1 accessions derived from 20 parents in the 132 inbred accessions). For later tasks, these averaged EFDs will be used.

・Dataset directory in this directory contain two datasets below.

1. Accession_list.csv: Lists of plant materials in our study
2. Raw_EFD_data.csv: Raw EFD data of the 291 accessions. In our study, EFDs were obtained from a total of 8,730 images (2 angles(Direction1 and Direction2) X 291 accesssions X 3 cultivation year X 5 fruits)


・R script in this directory enables calculation averaged EFDs for each accession in each angle from raw EFD

・This task doesn't need specific R packages.

# Part2.Genomic_and_Phenomic_prediction_of_Fruit_shape_using_Averaged_EFDs

・This directory is prepared for performing genomic and phenomic prediction methods (GP[20], GP[132], PPmid, PPδ).


・In each prediction methods, fruit contours in th 159 F1 accessions are predicted based on the averaged EFDs data of 20 F1 parents and other inbred accessions.


・Dataset directory in this directory contains two datasets below.

1. Accession_list.csv: Lists of plant materials in our study
2. Averaged_EFD_data.csv: Averaged EFD data to perform all prediction methods. This data is obtained in the task "Part1.Averaged_EFD_Calculation" as described above
3. Genotypic_data.csv: Numeric genotypic data (-1, 0, 1) of the 132 inbred accessions to perform GP[20] and GP[132] This data is derived from MIG-seq, resulting in 3,149 SNPs.
4. Parental_combinations_of_F1.csv: Lists for 159 F1 accessions with their parental accessions.

・R script in this directory enables prediction of EFDs in the 159 F1 accessions using GP[20], GP[132], PPmid, PPδ, respectively, and the illustration of fruit contours based on the real and predcted EFDs.

・This task needs three R packages "stringr", "rrBLUP", "RAINBOWR". 

# Part2.Genomic_and_Phenomic_prediction_of_Fruit_shape_using_Averaged_EFDs

・This directory is prepared for performing genomic and phenomic prediction methods (GP[20], GP[132], PPmid, PPδ).


・In each prediction methods, fruit contours in th 159 F1 accessions are predicted based on the averaged EFDs data of 20 F1 parents and other inbred accessions.


・Dataset directory in this directory contains two datasets below.

1. Accession_list.csv: Lists of plant materials in our study
2. Averaged_EFD_data.csv: Averaged EFD data to perform all prediction methods. This data is obtained in the task "Part1.Averaged_EFD_Calculation" as described above
3. Genotypic_data.csv: Numeric genotypic data (-1, 0, 1) of the 132 inbred accessions to perform GP[20] and GP[132] This data is derived from MIG-seq, resulting in 3,149 SNPs.
4. Parental_combinations_of_F1.csv: Lists for 159 F1 accessions with their parental accessions.

・R script in this directory enables prediction of EFDs in the 159 F1 accessions using GP[20], GP[132], PPmid, PPδ, respectively, and the illustration of fruit contours based on the real and predcted EFDs.

・This task needs three R packages "stringr", "rrBLUP", "RAINBOWR". 




