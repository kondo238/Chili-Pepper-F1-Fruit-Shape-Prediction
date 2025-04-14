# Chili-Pepper-F1-Fruit-Shape-Prediction
**Research paper**：*Prediction of fruit shapes in F1 progenies of chili peppers (Capsicum annuum) based on parental image data using elliptic Fourier analysis*

**Authors**：Fumiya Kondo, Yui Kumanomido, Mariasilvia D’Andrea, Valentino Palombo, Nahed Ahmed, Shino Futatsuyama, Kazuhiro Nemoto, Kenichi Matsushima

**Status**: Under review (Computer and Electronics in Agriculture, Elsevier, submitted 2025-04-11)

# Overview
This repository contains all datasets and R scripts used for implementing **phenomic and genomic prediction models** to estimate the fruit shapes of F1 chili pepper progenies (Capsicum annuum).

・Predictions are based on **Elliptic Fourier Descriptors (EFDs)** extracted from fruit images.
・F1 fruit shapes are predicted from parental EFDs, either:
    ・using **only image-derived EFDs** (*phenomic prediction*, PPδ), or
    ・combining EFDs with **genotypic data** (*genomic prediction*, GP).

Using the provided datasets and R scripts, the following four prediction strategies can be applied:

・**GP[20]**: Genomic prediction using EFDs and genotypic data from 20 F1 parents. Dominance effects in F1 progenies are mostly not considered.

・**GP[132]**: Genomic prediction using EFDs and genotypic data from 132 inbred accessions (including the 20 F1 parents). Dominance effects are mostly not considered.

・**PPmid**: Phenomic prediction using parental EFDs and a fixed dominance-to-additive effect ratio, accounting for dominance effects.

・**PPδ**: Pure phenomic prediction based solely on parental EFDs, without considering dominance.

Full methodological details are described in the manuscript (Kondo et al., under review).

# Repository Structure
The repository is organized into four main directories, each corresponding to a specific step in the F1 fruit shape prediction pipeline. Each directory is self-contained and can be used independently.

# Part1.Averaged_EFD_Calculation
**Purpose:** Calculate the average EFDs for 291 accessions (132 inbreds + 159 F1s) from raw image data. These averaged EFDs are used in subsequent prediction steps.

# Datasets:
・Accession_list.csv: List of plant materials in our study

・Raw_EFD_data.csv: Raw EFD data of the 291 accessions. In our study, EFDs were obtained from a total of 8,730 images (2 angles [Direction1 and Direction2] × 291 accessions × 3 cultivation years × 5 fruits)

# Script:
・Computes the average EFDs per accession per view direction.
・No R packages required.

# Part2.Genomic_and_Phenomic_prediction_of_Fruit_shape_using_Averaged_EFDs
・This directory is prepared for performing genomic and phenomic prediction methods (GP[20], GP[132], PPmid, PPδ).

・In each prediction method, fruit contours in the 159 F1 accessions are predicted based on the averaged EFDs data of 20 F1 parents and other inbred accessions.

・The dataset directory in this folder contains four datasets below:

Accession_list.csv: List of plant materials in our study

Averaged_EFD_data.csv: Averaged EFD data to perform all prediction methods. This data is obtained in the task "Part1.Averaged_EFD_Calculation" as described above

Genotypic_data.csv: Numeric genotypic data (-1, 0, 1) of the 132 inbred accessions to perform GP[20] and GP[132]. This data is derived from MIG-seq, resulting in 3,149 SNPs.

Parental_combinations_of_F1.csv: List of 159 F1 accessions with their parental accessions.

・The R script in this directory enables prediction of EFDs in the 159 F1 accessions using GP[20], GP[132], PPmid, and PPδ, and the illustration of fruit contours based on the real and predicted EFDs.

・This task needs three R packages: "stringr", "rrBLUP", and "RAINBOWR".

# Part3.Raw_EFD-based_phenomic_prediction_of_F1_fruit
・This directory is prepared for performing raw EFD data-based phenomic prediction using PPδ (showing the highest prediction accuracy compared to the other three methods [GP[20], GP[132], PPmid]).

・In each prediction method, fruit contours in the 159 F1 accessions are predicted based on the averaged EFDs data of 20 F1 parents and other inbred accessions.

・The dataset directory in this folder contains four datasets below:

Accession_list.csv: List of plant materials in our study

Parental_combinations_of_F1.csv: List of 159 F1 accessions with their parental accessions

Raw_EFD_data.csv: Raw EFD data of the 291 accessions. Same data used in the first task "Part1.Averaged_EFD_Calculation"

Representative_ratio_between_dominance_and_additive_effect.csv: Constant values to perform PPδ. These constant values are calculated in the former task "Part3.Raw_EFD-based_phenomic_prediction_of_F1_fruit"

・The R script in this directory enables raw EFD data-based PPδ to predict possible fruit contours in the 159 F1 accessions, and illustrates fruit contours based on the real and predicted EFDs.

・This task needs two R packages: "stringr" and "progress".

# Trial_for_F1_fruit_shape_prediction
・This directory is prepared for a test trial of PPδ using initial EFD data (".nef" format file) obtained from SAPE program software (This software extracts EFD data from raw image data [.jpg])

・Using raw EFD data (".nef") of two inbred accessions (Parent1 and Parent2) as examples, we can perform PPδ using the averaged and raw EFD.

・The detailed explanation is written in the R script in this directory.

・The dataset directory in this folder contains five datasets below:

Parent1_a.nef: Raw EFD data in Direction 1 (one of the two fruit angles) for Parent1

Parent1_b.nef: Raw EFD data in Direction 2 (the other angle) for Parent1

Parent2_a.nef: Raw EFD data in Direction 1 (one of the two fruit angles) for Parent2

Parent2_b.nef: Raw EFD data in Direction 2 (the other angle) for Parent2

Representative_ratio_between_dominance_and_additive_effect.csv: Constant values to perform PPδ. These constant values are calculated in the former task "Part3.Raw_EFD-based_phenomic_prediction_of_F1_fruit"

・The R script in this directory enables raw EFD data-based PPδ to predict possible fruit contours in the 159 F1 accessions, and illustrates fruit contours based on the real and predicted EFDs.

・This task doesn't need any R packages.
