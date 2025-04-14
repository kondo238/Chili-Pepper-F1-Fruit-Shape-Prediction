# Chili-Pepper-F1-Fruit-Shape-Prediction
**Research paper**：*Prediction of fruit shapes in F1 progenies of chili peppers (Capsicum annuum) based on parental image data using elliptic Fourier analysis*

**Authors**：Fumiya Kondo, Yui Kumanomido, Mariasilvia D’Andrea, Valentino Palombo, Nahed Ahmed, Shino Futatsuyama, Kazuhiro Nemoto, Kenichi Matsushima

**Status**: Under review (Computer and Electronics in Agriculture, Elsevier, submitted 2025-04-11)

# Overview
This repository contains all datasets and R scripts used for implementing **phenomic and genomic prediction models** to estimate the fruit shapes of F1 chili pepper progenies (Capsicum annuum).

* Predictions are based on **Elliptic Fourier Descriptors (EFDs)** extracted from fruit images.

* F1 fruit shapes are predicted from parental EFDs, either:
    * using **only image-derived EFDs** (*phenomic prediction*, PPδ), or
    * combining EFDs with **genotypic data** (*genomic prediction*, GP).

Using the provided datasets and R scripts, the following four prediction strategies can be applied:

* **GP[20]**: Genomic prediction using EFDs and genotypic data from 20 F1 parents. Dominance effects in F1 progenies are mostly not considered.

* **GP[132]**: Genomic prediction using EFDs and genotypic data from 132 inbred accessions (including the 20 F1 parents). Dominance effects are mostly not considered.

* **PPmid**: Phenomic prediction using parental EFDs and a fixed dominance-to-additive effect ratio, accounting for dominance effects.

* **PPδ**: Pure phenomic prediction based solely on parental EFDs, without considering dominance.

**Full methodological details are described in the manuscript (Kondo et al., under review).**

# Repository Structure
The repository is organized into four main directories, each corresponding to a specific step in the F1 fruit shape prediction pipeline. Each directory is self-contained and can be used independently.

# Part 1. Averaged EFD Calculation
**Purpose:** 
Calculate the average EFDs for 291 accessions (132 inbreds + 159 F1s) from raw image data. These averaged EFDs are used in subsequent prediction steps.

**Datasets:**
* Accession_list.csv: List of plant materials in our study
* Raw_EFD_data.csv: EFD data extracted from 8,730 fruit images (2 views × 291 accessions × 3 years × 5 fruits).

**Script:**
* Computes the average EFDs per accession per view direction.
* **No R packages required**.

# Part 2. Genomic and Phenomic prediction of Fruit shape using Averaged EFDs
**Purpose:** 
Perform genomic (GP[20], GP[132]) and phenomic (PPmid, PPδ) prediction using the averaged EFDs calculated in Part 1.
In each prediction method, fruit contours in the 159 F1 accessions are predicted based on the averaged EFDs data of 20 F1 parents and other inbred accessions.

**Datasets:**
* Accession_list.csv: List of plant materials in our study
* Averaged_EFD_data.csv: Generated in Part 1.
* Genotypic_data.csv: Genotype matrix (−1, 0, 1) from MIG-seq (3,149 SNPs).
* Parental_combinations_of_F1.csv: Links each F1 to its parental accessions.

**Script:**
・Predicts EFDs for the 159 F1 accessions.
・Visualizes fruit contours using both observed and predicted EFDs.
・**Required R packages**: stringr, rrBLUP, RAINBOWR.

# Part 3. Raw EFD-based phenomic prediction of F1 fruit
**Purpose:**
・Perform PPδ using raw EFD data, which showed the highest accuracy among all methods.
In each prediction method, fruit contours in the 159 F1 accessions are predicted based on the averaged EFDs data of 20 F1 parents and other inbred accessions.

**Datasets:**
・Accession_list.csv: List of plant materials in our study
・Parental_combinations_of_F1.csv: List of 159 F1 accessions with their parental accessions
・Raw_EFD_data.csv: Same as in Part 1.
・Representative_ratio_between_dominance_and_additive_effect.csv: Constant values to perform PPδ. These constant values are calculated in the former task "Part3.Raw_EFD-based_phenomic_prediction_of_F1_fruit"

**Script:**
・Performs PPδ using raw EFD data.
・Visualizes predicted fruit contours.
・**Required R packages**: stringr, progress. 

# Trial for F1 fruit shape prediction
**Purpose:**
Trial run of the PPδ method using .nef files (raw EFD data exported from SHAPE software)
[This software extracts EFD data from raw image data ".jpg"]

The detailed explanation is written in the R script in this directory.

**Example Data:**
・Parent1_a.nef, Parent1_b.nef: Two views of Parent1.
・Parent2_a.nef, Parent2_b.nef: Two views of Parent2.
・Representative_ratio_between_dominance_and_additive_effect.csv: Constant values to perform PPδ. These constant values are calculated in the former task "Part3.Raw_EFD-based_phenomic_prediction_of_F1_fruit"

**Script:**
・Demonstrates prediction and visualization of F1 contours based on .nef data.
・**No R packages required**.
