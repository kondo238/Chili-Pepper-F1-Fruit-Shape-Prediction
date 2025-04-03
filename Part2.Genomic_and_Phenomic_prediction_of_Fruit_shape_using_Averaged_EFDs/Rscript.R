######################Genomic and phenomic prediction of F1 fruit contours###################################################################
setwd(getwd())
library(stringr)
#library(SKM)
library(rrBLUP)
library(RAINBOWR)


#Data_loading
Acclist <- as.data.frame(read.csv("./Dataset/Accession_list.csv")) #Load_accession_list_data
Ave_data <- as.data.frame(read.csv("./Dataset/Averaged_EFD_data.csv"))# Load_raw_EFD_data_(80 variables)
Inbred_ID <- Acclist[Acclist$PopulationType == "Inbred",3]
combi <- as.data.frame(read.csv("./Dataset/Parental_combinations_of_F1.csv")) #Load_list_for_parental_combinations_of_F1
gt.score <- as.data.frame(read.csv("./Dataset/Genotypic_data.csv", row.names = 1)) #load_numeric_genotypic_data_of_132_inbred_accessions_(-1, 0, 1)_by_MIG-seq
gt.score <- gt.score[Inbred_ID,]








##############################################1.Genomic Prediction_GP[132]&GP[20]#############################################################
#1-1:Estimation_of_F1_genotypic_data_based_on_their_parental_genotypic_data
sim_mt <- as.data.frame(matrix(NA,
                               ncol = ncol(gt.score),
                               nrow = nrow(combi)))
colnames(sim_mt) <- colnames(gt.score)
rownames(sim_mt) <- combi[,3]


F1_ID <- as.character(combi[,3])
g <- gt.score
g_list <- rownames(g)

#Calculate_F1_genotypes(Calculating_Mean_of_parental_numeric_genotypic_data)
for(i in 1:length(F1_ID)){
  p1 <- as.numeric(g[which(g_list == combi[i,1]),]) #For mother parent
  p2 <- as.numeric(g[which(g_list == combi[i,2]),]) #For father parent

  #Before_Impute
  f1 <- (p1+p2)/2
  f1[which(f1 == 0.5)] <- 1
  f1[which(f1 == -0.5)] <- -1
  sim_mt[i,] <- f1
}

rm(g, f1, g_list, i, p1, p2)  








#1-2:GP[132](Genomic prediction for 159 F1 accessions based on 132 inbred accessions including the F1 parents)
#Data_preparation
Inbred_ave <- Ave_data[c(which(Acclist$PopulationType == "Inbred"),  #for direction 1
                         (length(Acclist$ID)+which(Acclist$PopulationType == "Inbred"))) #for direction 2
                         ,]

Pre <- Inbred_ave[0,] #Matrix_for_saving_predicted_values_in_GP[132](Genomic prediction for 159 F1 accessions based on 132 inbred accessions including the F1 parents)
Direction <- c("a", "b") #a=direction1 & b=direction2
all_ID <- c(Inbred_ID, F1_ID)
#genomic_prediction
for(x in 1:length(Direction)){
  #Preparation_of_training_EFD_data(from_132_inbred_accessions)
  pheno <- Inbred_ave[Inbred_ave$Direction == Direction[x],] 
  row.names(pheno) <- pheno$ID
  pheno <- pheno[Inbred_ID,-c(1,2)]
  a <- as.data.frame(matrix(NA,
                            ncol = ncol(pheno),
                            nrow = nrow(sim_mt)))
  
  colnames(a) <- colnames(pheno)
  pheno <- rbind(pheno, a)
  rownames(pheno) <- all_ID
  rm(a)
  
  #Preparation_of_training_genotypic_data(from 132 inbred accessions) and_test_genotypic_data(from_159_F1_accessions)
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
  
  #Calculate_gaussian_relation_matrix
  A <- as.matrix(gt.score2[all_ID,])
  A <- calcGRM(A,
               methodGRM = "gaussian",
               kernel.h = "tuned")
  
  
  #Preparation_of_matrix_for_saving_predicted_EFDs
  df <- as.data.frame(matrix(0,
                             ncol = ncol(pheno),
                             nrow = nrow(pheno)))
  

  #Genomic_prediction_using_GBLUP-GAUSS
  for(h in 1:ncol(pheno)){
  y_kin <- pheno[,h]
  #Calculation_GBLUP_as_Fixed_effect
  geno <- as.character(rownames(A))

  X_kin <- as.data.frame(cbind(geno, y_kin)) #Kernel_setting
  kin <- kin.blup(X_kin, geno="geno", pheno="y_kin", GAUSS=FALSE,K=A,fixed=NULL,covariate=NULL,
                  PEV=FALSE,n.core=1,theta.seq=NULL)
  
  #predict y based on X.test
  df[,h] <- as.numeric(kin$pred)
  }
  
  colnames(df) <- colnames(pheno)
  rownames(df) <- all_ID
  
  df <- cbind(Direction = rep(Direction[x], length(F1_ID)),
              ID = F1_ID,
              df[-c(1:length(Inbred_ID)),])
  
  Pre <- rbind(Pre, df)
  
}

Pre[,(1+2)] <- rep(1, nrow(Pre)) #Fill_1_as_constant_values_for_a1
Pre[,(21+2)] <- rep(0, nrow(Pre)) #Fill_0_as_constant_values_for_a1
Pre[,(41+2)] <- rep(0, nrow(Pre)) #Fill_0_as_constant_values_for_a1

GP132_Pre <- Pre

write.csv(GP132_Pre, "GP132_predicted_EFDs.csv", row.names = T) #Save_the_predicted_EFDs_as_csv_format

rm(A, df, gt.score2, kin, pheno, Pre, X_kin, Direction, geno, h, x, y_kin, all_ID)













#1-3:GP[20](Genomic prediction for 159 F1 accessions based on 20 F1 paretnts in inbred accessions)
#Data_preparation
F1parent_No <- which(str_detect(Acclist$Note, pattern = "F1 parent"))
F1parent_ID <- Acclist[F1parent_No,3]
F1parent_ave <- Ave_data[c(F1parent_No,  #for direction 1
                         (length(Acclist$ID)+F1parent_No)) #for direction 2
                       ,]

Pre <- F1parent_ave[0,] #Matrix_for_saving_predicted_values_in_GP[20](Genomic prediction for 159 F1 accessions based on 20 F1 paretnts in inbred accessions)
Direction <- c("a", "b") #a=direction1 & b=direction2
all_ID <- c(F1parent_ID, F1_ID)

#genomic_prediction
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
  
  #Preparation_of_training_genotypic_data(from 132 inbred accessions) and_test_genotypic_data(from_159_F1_accessions)
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
  
  #Calculate_gaussian_relation_matrix
  A <- as.matrix(gt.score2[all_ID,])
  A <- calcGRM(A,
               methodGRM = "gaussian",
               kernel.h = "tuned")
  
  
  #Preparation_of_matrix_for_saving_predicted_EFDs
  df <- as.data.frame(matrix(0,
                             ncol = ncol(pheno),
                             nrow = nrow(pheno))
  )
  
  
  #Genomic_prediction_using_GBLUP-GAUSS
  for(h in 1:ncol(pheno)){
    y_kin <- pheno[,h]
    #Calculation_GBLUP_as_Fixed_effect
    geno <- as.character(rownames(A))
    
    X_kin <- as.data.frame(cbind(geno, y_kin)) #Kernel_setting
    kin <- kin.blup(X_kin, geno="geno", pheno="y_kin", GAUSS=FALSE,K=A,fixed=NULL,covariate=NULL,
                    PEV=FALSE,n.core=1,theta.seq=NULL)
    
    #predict y based on X.test
    df[,h] <- as.numeric(kin$pred)
  }
  
  colnames(df) <- colnames(pheno)
  rownames(df) <- all_ID
  
  df <- cbind(Direction = rep(Direction[x], length(F1_ID)),
              ID = F1_ID,
              df[-c(1:length(F1parent_ID)),])
  
  Pre <- rbind(Pre, df)
  
}

Pre[,(1+2)] <- rep(1, nrow(Pre)) #Fill_1_as_constant_values_for_a1
Pre[,(21+2)] <- rep(0, nrow(Pre)) #Fill_0_as_constant_values_for_a1
Pre[,(41+2)] <- rep(0, nrow(Pre)) #Fill_0_as_constant_values_for_a1

GP20_Pre <- Pre

write.csv(GP20_Pre, "GP20_predicted_EFDs.csv", row.names = T) #Save_the_predicted_EFDs_as_csv_format

rm(A, df, gt.score2, kin, pheno, Pre, X_kin, Direction, geno, h, x, y_kin, all_ID, F1parent_No, gt.score)

















##############################################2.Phenomic Prediction_PPmid&PPδ#############################################################
#2-1:PPmid(phenomic prediction for 159 F1 accessions based on averaged EFDs of 20 F1 paretnts in inbred accessions)
#Data_preparation
F1parent_No <- which(str_detect(Acclist$Note, pattern = "F1 parent"))
F1parent_ID <- Acclist[F1parent_No,3]
F1parent_ave <- Ave_data[c(F1parent_No,  #for direction 1
                           (length(Acclist$ID)+F1parent_No)) #for direction 2
                         ,]
Pre <- F1parent_ave[0,] #Matrix_for_saving_predicted_values_in_PPmid(phenomic prediction for 159 F1 accessions based on averaged EFDs of 20 F1 paretnts in inbred accessions)
Direction <- c("a", "b") #a=direction1 & b=direction2

#Perform_PPmid
for(x in 1:length(Direction)){
  #Preparation_of_EFD_of_F1_parent_data(from_20_inbred_accessions)
  pheno <- F1parent_ave[F1parent_ave$Direction == Direction[x],] 
  row.names(pheno) <- pheno$ID
  pheno <- pheno[F1parent_ID,-c(1,2)]
  
  #Calculate_midpoint_EFD_as_predicted_EFDs_for_F1
  Pre1 <- F1parent_ave[0,]
  for(i in 1:length(F1_ID)){
    p1 <- as.character(combi[i,1]) #Mother_parent_ID
    p2 <- as.character(combi[i,2]) #Father_parent_ID
    p1h <- pheno[p1,] #Mother_parent_EFD
    p2h <- pheno[p2,] #Father_parent_EFD
      midpoint <- (p1h + p2h)/2 #Midpoint_calculation
      pre_F1 <- cbind(Direction=(Direction[x]),
                      ID=F1_ID[i],
                      midpoint)
      Pre1 <- rbind(Pre1, pre_F1)
  }
  Pre <- rbind(Pre, Pre1)
  rm(Pre1)
}

Pre[,(1+2)] <- rep(1, nrow(Pre)) #Fill_1_as_constant_values_for_a1
Pre[,(21+2)] <- rep(0, nrow(Pre)) #Fill_0_as_constant_values_for_a1
Pre[,(41+2)] <- rep(0, nrow(Pre)) #Fill_0_as_constant_values_for_a1

PPmid_Pre <- Pre
write.csv(PPmid_Pre, "PPmid_predicted_EFDs.csv", row.names = T) #Save_the_predicted_EFDs_as_csv_format

rm(F1parent_ave, midpoint, p1h, p2h, Pre, pre_F1, F1parent_No, i, p1, p2, x)







#2-2:PPδ(phenomic prediction for 159 F1 accessions based on averaged EFDs of 20 F1 paretnts in inbred accessions and representative ratio between dominance and additive effects)
#Data_preparation
F1parent_No <- which(str_detect(Acclist$Note, pattern = "F1 parent"))
F1parent_ID <- Acclist[F1parent_No,3]
F1parent_ave <- Ave_data[c(F1parent_No,  #for direction 1
                           (length(Acclist$ID)+F1parent_No)) #for direction 2
                         ,]
F1_ave <- Ave_data[c(which(Acclist$PopulationType == "F1"),  #for direction 1
                     (length(Acclist$ID)+which(Acclist$PopulationType == "F1"))) #for direction 2
                   ,]
Pre <- F1parent_ave[0,] #Matrix_for_saving_predicted_values_in_PPmid(phenomic prediction for 159 F1 accessions based on averaged EFDs of 20 F1 paretnts in inbred accessions)
Direction <- c("a", "b") #a=direction1 & b=direction2

#Calculation_of_1.additive_effect and 2.dominance effects and 3. their ratio 159 F1 accessions.
A <- D <- R <- as.data.frame(F1_ave[,])
A[,-c(1:2)] <- D[,-c(1:2)] <- R[,-c(1:2)] <- NA #Prepare matrix for saving additive effect, dominance effect, and their ratio


for(x in 1:length(Direction)){
  pheno <- F1parent_ave[F1parent_ave$Direction == Direction[x],] 
  row.names(pheno) <- pheno$ID
  pheno <- pheno[F1parent_ID,-c(1,2)]
  
  f1_pheno <- F1_ave[F1_ave$Direction == Direction[x],]
  row.names(f1_pheno) <- f1_pheno$ID
  f1_pheno <- f1_pheno[F1_ID,-c(1,2)]
  
  for(i in 1:length(F1_ID)){
    p1 <- as.character(combi[i,1]) #Mother_parent_ID
    p2 <- as.character(combi[i,2]) #Father_parent_ID
    f1 <- as.character(combi[i,3]) #F1_ID
    p1h <- pheno[p1,] #Mother_parent_EFD
    p2h <- pheno[p2,] #Father_parent_EFD
    f1h <- f1_pheno[f1,] #F1_EFD
    midpoint <- (p1h + p2h)/2 #Midpoint_calculation
    a <- abs(p1h-p2h)/2 #Additive_effect_calculation
    d <- f1h - midpoint #Dominance_effect_calculation
    r <- d/a #Ratio_between_additive_and_dominance_effects_calculation
    
    A[A$Direction == Direction[x] &
        A$ID == F1_ID[i],c(3:82)] <- a
    D[D$Direction == Direction[x] &
        D$ID == F1_ID[i],c(3:82)] <- d
    R[R$Direction == Direction[x] &
        R$ID == F1_ID[i],c(3:82)] <- r
    }
}

rm(a, d, f1_pheno, f1h, p1h, p2h, pheno, r, f1, i, p1, p2)

#Determine_the_representative_ratio_between_additive_and_dominance_effect_among_156_F1_accessions
Rep_R <- as.data.frame(F1_ave[1:2,])
Rep_R[1:2,1] <- c("a", "b")
Rep_R[1:2,2] <- rep("Representative", 2)
Rep_R[,-c(1:2)] <- NA

for(i in 1:80){ #calculate_median_of_the_ratio
  Rep_R[1,(i+2)] <- median(R[R$Direction == "a",(i+2)])
  Rep_R[2,(i+2)] <- median(R[R$Direction == "b",(i+2)])
}
  
write.csv(Rep_R, "Representative_ratio_between_dominance_and_additive_effect.csv", row.names = F)  
rm(i, A, D, R)

#Perform_PPδ
for(x in 1:length(Direction)){
  #Preparation_of_EFD_of_F1_parent_data(from_20_inbred_accessions)
  pheno <- F1parent_ave[F1parent_ave$Direction == Direction[x],] 
  row.names(pheno) <- pheno$ID
  pheno <- pheno[F1parent_ID,-c(1,2)]
  
  #Prediction_of_F1_EFDs
  Pre1 <- F1parent_ave[0,]
  for(i in 1:length(F1_ID)){
    p1 <- as.character(combi[i,1]) #Mother_parent_ID
    p2 <- as.character(combi[i,2]) #Father_parent_ID
    p1h <- pheno[p1,] #Mother_parent_EFD
    p2h <- pheno[p2,] #Father_parent_EFD
    midpoint <- (p1h + p2h)/2 #Midpoint_calculation
    a <- abs(p1h-p2h)/2 #Additive_effect_calculation
    d_est <- a*Rep_R[Rep_R$Direction == "a",-c(1:2)] #Estimate_dominance_effects_in_the_crossing_combination
    pre_pheno <- midpoint + d_est #Calculate_predicted_EFDs_for_the_F1_by_adding_estimated_dominance_effects_on_the_midpoint_EFDs
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
write.csv(PPdelta_Pre, "PPdelta_predicted_EFDs.csv", row.names = T) #Save_the_predicted_EFDs_as_csv_format


rm(a, d_est, midpoint, p1h, p2h, pheno, Pre, pre_F1, Pre1, Direction, F1parent_No, p1, p2, x, i, pre_pheno)














#################################3.Draw_fruit_contours_based_on_predicted_and_real_EFDs######################################################
#Define_the_Elliptic Fourier function
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


#Setting_for_illustration
lw <- 2 #line thickness
Direction <- "b" #Determine_the_direction (a=Direction1 & b=Direction2)
par(mfrow = c(10, 16), mar = c(0.0, 0.0, 0.0, 0.0)) #Determine_the_number_of_row(reft)_and_column(right)


#Draw_real_and_predicted_EFD-based fruit contours
for(i in 1:length(F1_ID)){
  #Draw_Real_EFD-based_fruit_contours(Black)
  x <- F1_ave[F1_ave$Direction == Direction &
                F1_ave$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "black", lwd = lw)
  
  #Draw_predicted_EFD(GP[132])-based_fruit_contours(Red)
  x <- GP132_Pre[GP132_Pre$Direction == Direction &
                   GP132_Pre$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  lines(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "red", lwd = lw)
  
  #Draw_predicted_EFD(GP[20])-based_fruit_contours(Green)
  x <- GP20_Pre[GP20_Pre$Direction == Direction &
                   GP20_Pre$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  lines(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
        ann = F, axes = F, col = "green4", lwd = lw)
  
  #Draw_predicted_EFD(PPmid)-based_fruit_contours(Gold)
  x <- PPmid_Pre[PPmid_Pre$Direction == Direction &
                  PPmid_Pre$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  lines(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
        ann = F, axes = F, col = "gold2", lwd = lw)
  
  #Draw_predicted_EFD(PPdelta)-based_fruit_contours(Blue)
  x <- PPdelta_Pre[PPdelta_Pre$Direction == Direction &
                   PPdelta_Pre$ID == F1_ID[i],-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  lines(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
        ann = F, axes = F, col = "blue", lwd = lw)
  }

rm(x, ef, coord, lw, i, ef2coord, Direction)



