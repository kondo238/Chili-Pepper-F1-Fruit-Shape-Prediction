###################### Part 2. Genomic and phenomic prediction of F1 fruit contours ###################################################################

# Set working directory
setwd(getwd())

# Load required libraries 
library(stringr)
library(rrBLUP)
library(RAINBOWR)

# Load datasets
Acclist <- as.data.frame(read.csv("./Dataset/Accession_list.csv")) # Load accession list
Ave_data <- as.data.frame(read.csv("./Dataset/Averaged_EFD_data.csv")) # Load averaged EFD data
Inbred_ID <- Acclist[Acclist$PopulationType == "Inbred",3]
combi <- as.data.frame(read.csv("./Dataset/Parental_combinations_of_F1.csv")) # Load parental combinations
gt.score <- as.data.frame(read.csv("./Dataset/Genotypic_data.csv", row.names = 1)) # Load genotypic data
# Filter genotypic data for inbred accessions
gt.score <- gt.score[Inbred_ID,]

############################################## 1.Genomic Prediction GP[132] & GP[20] #############################################################

# 1-1: Estimate F1 genotypic data based on parental genotypic data
sim_mt <- as.data.frame(matrix(NA,
                               ncol = ncol(gt.score),
                               nrow = nrow(combi)))
colnames(sim_mt) <- colnames(gt.score)
rownames(sim_mt) <- combi[,3]

F1_ID <- as.character(combi[,3])
g <- gt.score
g_list <- rownames(g)

# Compute F1 genotypes by averaging parental genotypic data
for(i in 1:length(F1_ID)){
  p1 <- as.numeric(g[which(g_list == combi[i,1]),]) # Mother parent
  p2 <- as.numeric(g[which(g_list == combi[i,2]),]) # Father parent

  # Before_Impute
  f1 <- (p1+p2)/2
  f1[which(f1 == 0.5)] <- 1
  f1[which(f1 == -0.5)] <- -1
  sim_mt[i,] <- f1
}
rm(g, f1, g_list, i, p1, p2)  # Clean up workspace

# 1-2: GP[132] - Genomic prediction for 159 F1 accessions based on 132 inbred parents
Inbred_ave <- Ave_data[c(which(Acclist$PopulationType == "Inbred"),  #for direction 1
                         (length(Acclist$ID)+which(Acclist$PopulationType == "Inbred"))) #for direction 2
                         ,]

Pre <- Inbred_ave[0,] # Matrix for storing predicted values
all_ID <- c(Inbred_ID, F1_ID)

# Perform genomic prediction
for(x in 1:length(Direction)){
  #Preparation_of_training_EFD_data(from_132_inbred_accessions)
  pheno <- Inbred_ave[Inbred_ave$Direction == Direction[x],] 
  row.names(pheno) <- pheno$ID
  pheno <- pheno[Inbred_ID,-c(1,2)]

  # Initialize prediction matrix
  a <- as.data.frame(matrix(NA,
                            ncol = ncol(pheno),
                            nrow = nrow(sim_mt)))
  
  colnames(a) <- colnames(pheno)
  pheno <- rbind(pheno, a)
  rownames(pheno) <- all_ID
  rm(a)
  
  #Prepare genotypic data
  gt.score2 <- rbind(gt.score, sim_mt)
  X <- gt.score2
  poly <- apply(X, 2, var) > 0
  X <- X[, poly]
  
  A <- A.mat(X,
             min.MAF = 0.05,
             max.missing = NULL,
             return.imputed = TRUE)
  X <- A$imputed
  gt.score2 <- X
  rm(X, A, poly)
  
  # Compute Gaussian relationship matrix
  A <- as.matrix(gt.score2[all_ID,])
  A <- calcGRM(A,
               methodGRM = "gaussian",
               kernel.h = "tuned")

  # Initialize prediction result matrix
  df <- as.data.frame(matrix(0,
                             ncol = ncol(pheno),
                             nrow = nrow(pheno)))
  
  # Genomic prediction using GBLUP-GAUSS
  for(h in 1:ncol(pheno)){
  y_kin <- pheno[,h]
  geno <- as.character(rownames(A))
  X_kin <- as.data.frame(cbind(geno, y_kin)) #Kernel_setting
  kin <- kin.blup(X_kin, geno="geno", pheno="y_kin", GAUSS=FALSE,K=A,fixed=NULL,covariate=NULL,
                  PEV=FALSE,n.core=1,theta.seq=NULL)
    df[,h] <- as.numeric(kin$pred)
  }
  
  colnames(df) <- colnames(pheno)
  rownames(df) <- all_ID
  df <- cbind(Direction = rep(Direction[x], length(F1_ID)),
              ID = F1_ID,
              df[-c(1:length(Inbred_ID)),])
  Pre <- rbind(Pre, df)
}

Pre[,(1+2)] <- rep(1, nrow(Pre)) # Fill 1 as constant values for a1
Pre[,(21+2)] <- rep(0, nrow(Pre)) #Fill 0 as constant values for a1
Pre[,(41+2)] <- rep(0, nrow(Pre)) #Fill 0 as constant values for a1

GP132_Pre <- Pre

# Save predictions
write.csv(GP132_Pre, "GP132_predicted_EFDs.csv", row.names = T) 
rm(A, df, gt.score2, kin, pheno, Pre, X_kin, Direction, geno, h, x, y_kin, all_ID) # Clean up workspace

# 1-3: GP[20] - Genomic prediction for 159 F1 accessions based on 20 F1 parents in inbred accessions
# Data preparation
F1parent_No <- which(str_detect(Acclist$Note, pattern = "F1 parent"))
F1parent_ID <- Acclist[F1parent_No,3]
F1parent_ave <- Ave_data[c(F1parent_No,  # for direction 1
                         (length(Acclist$ID)+F1parent_No)) # for direction 2
                       ,]

Pre <- F1parent_ave[0,] #Matrix_for_saving_predicted_values_in_GP[20](Genomic prediction for 159 F1 accessions based on 20 F1 paretnts in inbred accessions)
Direction <- c("a", "b") # a=direction1 & b=direction2
all_ID <- c(F1parent_ID, F1_ID)

# Genomic prediction
for(x in 1:length(Direction)){
  #Preparation_of_training_EFD_data(from_20_inbred_accessions)
  pheno <- F1parent_ave[F1parent_ave$Direction == Direction[x],] 
  row.names(pheno) <- pheno$ID
  pheno <- pheno[F1parent_ID,-c(1,2)]
  a <- as.data.frame(matrix(NA,
                            ncol = ncol(pheno),
                            nrow = nrow(sim_mt)))
  
  colnames(a) <- colnames(pheno)
  pheno <- rbind(pheno, a)
  rownames(pheno) <- all_ID
  rm(a)
  
  # Preparation of training genotypic data (from 132 inbred accessions) and test genotypic data (from 159 F1 accessions)
  gt.score2 <- rbind(gt.score[F1parent_ID,], sim_mt)
  
  X <- gt.score2
  poly <- apply(X, 2, var) > 0
  X <- X[, poly]
  
  A <- A.mat(X,
             min.MAF = 0.05,
             max.missing = NULL,
             return.imputed = TRUE)
  X <- A$imputed
  gt.score2 <- X
  rm(X, A, poly)
  
  # Compute gaussian relation matrix
  A <- as.matrix(gt.score2[all_ID,])
  A <- calcGRM(A,
               methodGRM = "gaussian",
               kernel.h = "tuned")
  
  # Initialize prediction result matrix
  df <- as.data.frame(matrix(0,
                             ncol = ncol(pheno),
                             nrow = nrow(pheno))
  )
  # Genomic prediction using GBLUP-GAUSS
  for(h in 1:ncol(pheno)){
    y_kin <- pheno[,h]
    # Calculation GBLUP as Fixed effect
    geno <- as.character(rownames(A))
    
    X_kin <- as.data.frame(cbind(geno, y_kin)) # Kernel_setting
    kin <- kin.blup(X_kin, geno="geno", pheno="y_kin", GAUSS=FALSE,K=A,fixed=NULL,covariate=NULL,
                    PEV=FALSE,n.core=1,theta.seq=NULL)
    # Predict y based on X.test
    df[,h] <- as.numeric(kin$pred)
  }
  
  colnames(df) <- colnames(pheno)
  rownames(df) <- all_ID
  
  df <- cbind(Direction = rep(Direction[x], length(F1_ID)),
              ID = F1_ID,
              df[-c(1:length(F1parent_ID)),])
  
  Pre <- rbind(Pre, df)
  
}

Pre[,(1+2)] <- rep(1, nrow(Pre)) # Fill 1 as constant values for a1
Pre[,(21+2)] <- rep(0, nrow(Pre)) # Fill 0 as constant values for a1
Pre[,(41+2)] <- rep(0, nrow(Pre)) # Fill 0 as constant values for a1

GP20_Pre <- Pre

# Save predictions
write.csv(GP20_Pre, "GP20_predicted_EFDs.csv", row.names = T) #Save_the_predicted_EFDs_as_csv_format
rm(A, df, gt.score2, kin, pheno, Pre, X_kin, Direction, geno, h, x, y_kin, all_ID, F1parent_No, gt.score) # Clean up workspace


############################################## 2. Phenomic Prediction PPmid & PPδ #############################################################

# 2-1: PPmid - phenomic prediction for 159 F1 accessions based on averaged EFDs of 20 F1 paretnts in inbred accessions
# Data preparation
F1parent_No <- which(str_detect(Acclist$Note, pattern = "F1 parent"))
F1parent_ID <- Acclist[F1parent_No,3]
F1parent_ave <- Ave_data[c(F1parent_No,  #for direction 1
                           (length(Acclist$ID)+F1parent_No)) # for direction 2
                         ,]
Pre <- F1parent_ave[0,] # Matrix for saving predicted values in PPmid (phenomic prediction for 159 F1 accessions based on averaged EFDs of 20 F1 paretnts in inbred accessions)
Direction <- c("a", "b") # a=direction1 & b=direction2

# Perform PPmid
for(x in 1:length(Direction)){
  # Preparation of EFD of F1 parent data (from 20 inbred accessions)
  pheno <- F1parent_ave[F1parent_ave$Direction == Direction[x],] 
  row.names(pheno) <- pheno$ID
  pheno <- pheno[F1parent_ID,-c(1,2)]
  
  # Compute midpoint EFD as predicted EFDs for F1
  Pre1 <- F1parent_ave[0,]
  for(i in 1:length(F1_ID)){
    p1 <- as.character(combi[i,1]) # Mother parent ID
    p2 <- as.character(combi[i,2]) # Father parent ID
    p1h <- pheno[p1,] # Mother parent EFD
    p2h <- pheno[p2,] # Father parent EFD
      midpoint <- (p1h + p2h)/2 # Midpoint calculation
      pre_F1 <- cbind(Direction=(Direction[x]),
                      ID=F1_ID[i],
                      midpoint)
      Pre1 <- rbind(Pre1, pre_F1)
  }
  Pre <- rbind(Pre, Pre1)
  rm(Pre1)
}

Pre[,(1+2)] <- rep(1, nrow(Pre)) # Fill 1 as constant values for a1
Pre[,(21+2)] <- rep(0, nrow(Pre)) # Fill 0 as constant values for a1
Pre[,(41+2)] <- rep(0, nrow(Pre)) # Fill 0 as constant values for a1

PPmid_Pre <- Pre
# Save predictions
write.csv(PPmid_Pre, "PPmid_predicted_EFDs.csv", row.names = T) #Save_the_predicted_EFDs_as_csv_format
rm(F1parent_ave, midpoint, p1h, p2h, Pre, pre_F1, F1parent_No, i, p1, p2, x) # Clean up workspace


# 2-2: PPδ - phenomic prediction for 159 F1 accessions based on averaged EFDs of 20 F1 paretnts in inbred accessions and representative ratio between dominance and additive effects
# Data preparation
F1parent_No <- which(str_detect(Acclist$Note, pattern = "F1 parent"))
F1parent_ID <- Acclist[F1parent_No,3]
F1parent_ave <- Ave_data[c(F1parent_No,  # for direction 1
                           (length(Acclist$ID)+F1parent_No)) # for direction 2
                         ,]
F1_ave <- Ave_data[c(which(Acclist$PopulationType == "F1"),  # for direction 1
                     (length(Acclist$ID)+which(Acclist$PopulationType == "F1"))) # for direction 2
                   ,]
Pre <- F1parent_ave[0,] # Matrix for saving predicted values in PPmid (phenomic prediction for 159 F1 accessions based on averaged EFDs of 20 F1 paretnts in inbred accessions)
Direction <- c("a", "b") # a=direction1 & b=direction2

# Compute (I) additive effect, (II) dominance effects and (III) their ratio 159 F1 accessions
A <- D <- R <- as.data.frame(F1_ave[,])
A[,-c(1:2)] <- D[,-c(1:2)] <- R[,-c(1:2)] <- NA # Prepare matrix for saving additive effect, dominance effect, and their ratio

for(x in 1:length(Direction)){
  pheno <- F1parent_ave[F1parent_ave$Direction == Direction[x],] 
  row.names(pheno) <- pheno$ID
  pheno <- pheno[F1parent_ID,-c(1,2)]
  
  f1_pheno <- F1_ave[F1_ave$Direction == Direction[x],]
  row.names(f1_pheno) <- f1_pheno$ID
  f1_pheno <- f1_pheno[F1_ID,-c(1,2)]
  
  for(i in 1:length(F1_ID)){
    p1 <- as.character(combi[i,1]) # Mother parent ID
    p2 <- as.character(combi[i,2]) # Father parent ID
    f1 <- as.character(combi[i,3]) # F1 ID
    p1h <- pheno[p1,] # Mother parent EFD
    p2h <- pheno[p2,] # Father parent EFD
    f1h <- f1_pheno[f1,] # F1 EFD
    midpoint <- (p1h + p2h)/2 # Midpoint calculation
    a <- abs(p1h-p2h)/2 # Additive effect calculation
    d <- f1h - midpoint # Dominance effect calculation
    r <- d/a # Ratio between additive and dominance effects calculation
    
    A[A$Direction == Direction[x] &
        A$ID == F1_ID[i],c(3:82)] <- a
    D[D$Direction == Direction[x] &
        D$ID == F1_ID[i],c(3:82)] <- d
    R[R$Direction == Direction[x] &
        R$ID == F1_ID[i],c(3:82)] <- r
    }
}

rm(a, d, f1_pheno, f1h, p1h, p2h, pheno, r, f1, i, p1, p2) # Clean up workspace

# Determine the representative ratio between additive and dominance effect among 156 F1 accessions
Rep_R <- as.data.frame(F1_ave[1:2,])
Rep_R[1:2,1] <- c("a", "b")
Rep_R[1:2,2] <- rep("Representative", 2)
Rep_R[,-c(1:2)] <- NA

for(i in 1:80){ #calculate_median_of_the_ratio
  Rep_R[1,(i+2)] <- median(R[R$Direction == "a",(i+2)])
  Rep_R[2,(i+2)] <- median(R[R$Direction == "b",(i+2)])
}

# Save results
write.csv(Rep_R, "Representative_ratio_between_dominance_and_additive_effect.csv", row.names = F)  
rm(i, A, D, R) # Clean up workspace

# Perform PPδ
for(x in 1:length(Direction)){
  # Preparation of EFD of F1 parent data (from 20 inbred accessions)
  pheno <- F1parent_ave[F1parent_ave$Direction == Direction[x],] 
  row.names(pheno) <- pheno$ID
  pheno <- pheno[F1parent_ID,-c(1,2)]
  
  # Prediction of F1 EFDs
  Pre1 <- F1parent_ave[0,]
  for(i in 1:length(F1_ID)){
    p1 <- as.character(combi[i,1]) # Mother parent ID
    p2 <- as.character(combi[i,2]) # Father parent ID
    p1h <- pheno[p1,] # Mother parent EFD
    p2h <- pheno[p2,] # Father parent EFD
    midpoint <- (p1h + p2h)/2 # Midpoint calculation
    a <- abs(p1h-p2h)/2 # Additive effect calculation
    d_est <- a*Rep_R[Rep_R$Direction == Direction[x],-c(1:2)] # Estimate dominance effects in the crossing combination
    pre_pheno <- midpoint + d_est # Calculate predicted EFDs for the F1 by adding estimated dominance effects on the midpoint EFDs
    pre_F1 <- cbind(Direction=(Direction[x]),
                    ID=F1_ID[i],
                    pre_pheno)
  Pre1 <- rbind(Pre1, pre_F1)
  }
  Pre <- rbind(Pre, Pre1)
}

Pre[,(1+2)] <- rep(1, nrow(Pre)) #Fill_1_as_constant_values_for_a1
Pre[,(21+2)] <- rep(0, nrow(Pre)) #Fill_0_as_constant_values_for_a1
Pre[,(41+2)] <- rep(0, nrow(Pre)) #Fill_0_as_constant_values_for_a1

PPdelta_Pre <- Pre
# Save predictions
write.csv(PPdelta_Pre, "PPdelta_predicted_EFDs.csv", row.names = T) #Save_the_predicted_EFDs_as_csv_format
rm(a, d_est, midpoint, p1h, p2h, pheno, Pre, pre_F1, Pre1, Direction, F1parent_No, p1, p2, x, i, pre_pheno) # Clean up workspace


################################# 3. Draw fruit contours based on predicted and real EFDs ######################################################

# Define the Elliptic Fourier function
ef2coord <- function(ef, theta = seq(0, 2*pi, 0.01)) {
  x <- 0; y <- 0
  z <- length(ef)/4
  for(i in 1:(length(ef)/4)) {
    x <- x + ef[i] * cos(i * theta) + ef[i+z] * sin(i * theta)
    y <- y + ef[i+z*2] * cos(i * theta) + ef[i + z*3] * sin(i * theta)
  }
  coord <- list(x = x, y = y)
  return(coord)
}

# Setting for illustration
lw <- 2 # Line thickness
Direction <- "b" # Determine the direction (a=Direction1 & b=Direction2)
par(mfrow = c(10, 16), mar = c(0.0, 0.0, 0.0, 0.0)) # Determine the number of row (left) and column (right)

# Draw real and predicted EFD-based fruit contours
for(i in 1:length(F1_ID)){
  # Draw Real EFD-based fruit contours (Black)
  x <- F1_ave[F1_ave$Direction == Direction &
                F1_ave$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "black", lwd = lw)
  
  # Draw predicted EFD (GP[132]) - based fruit contours (Red)
  x <- GP132_Pre[GP132_Pre$Direction == Direction &
                   GP132_Pre$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  lines(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "red", lwd = lw)
  
  # Draw predicted EFD (GP[20]) - based fruit contours (Green)
  x <- GP20_Pre[GP20_Pre$Direction == Direction &
                   GP20_Pre$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  lines(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
        ann = F, axes = F, col = "green4", lwd = lw)
  
  # Draw predicted EFD (PPmid) - based fruit contours (Gold)
  x <- PPmid_Pre[PPmid_Pre$Direction == Direction &
                  PPmid_Pre$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  lines(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
        ann = F, axes = F, col = "gold2", lwd = lw)
  
  # Draw predicted EFD (PPdelta) - based fruit contours (Blue)
  x <- PPdelta_Pre[PPdelta_Pre$Direction == Direction &
                   PPdelta_Pre$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  lines(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
        ann = F, axes = F, col = "blue", lwd = lw)
  }

rm(x, ef, coord, lw, i, ef2coord, Direction)
