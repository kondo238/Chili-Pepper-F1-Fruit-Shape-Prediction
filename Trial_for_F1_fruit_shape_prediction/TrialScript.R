##########Trial_script:Perform_fruit_shape_prediction_using_PPδ######################################################################
######################NOTE!!!################################################################################################################################
#1.This trial uses EFD data of the two accessions (Parent1 and Parent2) as parental examples and predicts their F1 fruit shape contours
#2.The EFD data was obtained by SHAPE program saved as ".nef" format (The number of harmonics was 20)
#3.In the shape program, EFDs for one direction (Direction 1 or Direction 2) were extracted from five fruits per one accession and saved as same ".nef" file
#4.For implementation of PPδ, representative ratio between dominance and additive_effect derived from the 156 F1 accession in author's research.
#5. fruit shape prediction was performed based on parental raw EFD and averaged EFD data

#################################Step.1_Data_loading_from_".nef"_format_file_and_Data_formatting###################################################################################################################
#Parental_raw_EFD_data_was_saved_in_"df"
#Parental_Averaged_EFD_data_was_saved_in_"Ave_df"
setwd(getwd())
ID <- c("Parent1", "Parent2") #Define parental ID
Direction <- c("a", "b") #Define direction (a=Direction1, b=Direction2)
Replicate <- c("P1", "P2", "P3", "P4", "P5") #Define_fruit_replicates_(In this trial, we have five data for one direction per one accession)

#Load_EFD_data
Parent1_a <- as.data.frame(read.delim("./Dataset/Parent1_a.nef", sep = "")) #Load_EFD_data(Direction 1)_from_five_fruits_of_Parent1_using_SHAPE_program
Parent1_b <- as.data.frame(read.delim("./Dataset/Parent1_b.nef", sep = "")) #Load_EFD_data(Direction 2)_from_five_fruits_of_Parent1_using_SHAPE_program

Parent2_a <- as.data.frame(read.delim("./Dataset/Parent2_a.nef", sep = "")) #Load_EFD_data(Direction 1)_from_five_fruits_of_Parent2_using_SHAPE_program
Parent2_b <- as.data.frame(read.delim("./Dataset/Parent2_b.nef", sep = "")) #Load_EFD_data(Direction 2)_from_five_fruits_of_Parent2_using_SHAPE_program

#Data_collection_formatting_for_raw_EFD_data
df <- as.data.frame(matrix(NA,
                           ncol = 83,
                           nrow = length(ID)*length(Direction)*length(Replicate)))
colnames(df) <- c("Direction", "ID", "Replicate",
                  c(paste(rep("a", 20), 1:20, sep = ""),
                    paste(rep("b", 20), 1:20, sep = ""),
                    paste(rep("c", 20), 1:20, sep = ""),
                    paste(rep("d", 20), 1:20, sep = "")))
df[,1] <- c(rep("a", length(ID)*length(Replicate)),
            rep("b", length(ID)*length(Replicate)))
df[,2] <- rep(c(rep("Parent1", length(Replicate)),
                rep("Parent2", length(Replicate))),
              length(Direction))
df[,3] <- rep(Replicate, length(Direction)*length(ID))

for(p in 1:length(Replicate)){
  no <- 2 + (p - 1)*21
  mt_p1_a <- Parent1_a[(no:(no+20)),-5]
  mt_p1_b <- Parent1_b[(no:(no+20)),-5] 
  mt_p2_a <- Parent2_a[(no:(no+20)),-5] 
  mt_p2_b <- Parent2_b[(no:(no+20)),-5] 
  rm(no)
  for(i in 1:20){
  #For_Parent1
    no <- 4 + (i -1)*1
  df[df$Direction == "a" &
       df$ID == "Parent1" &
       df$Replicate == Replicate[p],c(no, no+20, no+40, no+60)] <- mt_p1_a[(i+1),]
  
  df[df$Direction == "b" &
       df$ID == "Parent1" &
       df$Replicate == Replicate[p],c(no, no+20, no+40, no+60)] <- mt_p1_b[(i+1),]
  
  df[df$Direction == "a" &
       df$ID == "Parent2" &
       df$Replicate == Replicate[p],c(no, no+20, no+40, no+60)] <- mt_p2_a[(i+1),]
  
  df[df$Direction == "b" &
       df$ID == "Parent2" &
       df$Replicate == Replicate[p],c(no, no+20, no+40, no+60)] <- mt_p2_b[(i+1),]
  }
}

df[,(3+1)] <- 1
df[,(3+21)] <- 0
df[,(3+41)] <- 0

for(i in 1:80){
  df[,(i+3)] <- as.numeric(df[,(i+3)])
}


#Data_collection_of_formatting_for_averaged_EFD_data
Ave_df <- as.data.frame(df[0,-3])
for(j in 1:length(ID)){
    for(i in 1:length(Direction)){
  df2 <- df[df$Direction == Direction[i] & df$ID == ID[j],-c(1:3)]
  m <- (df2[1,] + df2[2,] + df2[3,] + df2[4,] + df2[5,])/5
  ave <- df[0,-3]
  ave[1,] <- c(Direction[i], ID[j], m)
  Ave_df <- rbind(Ave_df, ave)
  }
}



rm(mt_p1_a, mt_p1_b, mt_p2_a, mt_p2_b, Parent1_a, Parent1_b, Parent2_a, Parent2_b, i, no, p, df2, m, j, ave, ID, Replicate, Direction)



#################################2.Draw_parental_fruit_contours_based_on_the_raw_and_averaged_EFD#####################################################
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

ID <- c("Parent1", "Parent2") #Define parental ID
Direction <- c("a", "b") #Define direction (a=Direction1, b=Direction2)
Replicate <- c("P1", "P2", "P3", "P4", "P5") #Define_fruit_replicates_(In this trial, we have five data for one direction per one accession)

#Setting_for_illustration
lw <- 2 #line thickness
par(mfrow = c(length(ID), 12), mar = c(0.0, 0.0, 0.0, 0.0)) #Determine_the_number_of_row(reft)_and_column(right)

#Draw_raw_and_averaged_EFD-based fruit contours of the two parents(Parent1 and Parent2)
#MEMO_in_plotted_figure:Averaged EFD-based contours were drawn as red lines & Raw EFD-based contours were drawn as black lines
#MEMO_in_plotted_figure:The contours in the first row showed parent1 and those in second row showed parent2
#MEMO_in_plotted_figure:For_raw_EFD-based_contours,_left_five_contours_showed_direction1_and_right_five_contours_showed_direction2.

for(i in 1:length(ID)){
  mt <- df[df$ID == ID[i],]
  ave_mt <- Ave_df[Ave_df$ID == ID[i],]
  
  #For Averaged EFD (Direction 1)
  x <- ave_mt[ave_mt$Direction =="a",-c(1:2)]
  x <- as.numeric(x[1,])
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "red", lwd = lw)
  
  #For raw EFD (Direction 1)
  for(r in 1:length(Replicate)){
  x <- mt[mt$Direction =="a",-c(1:3)]
  x <- as.numeric(x[r,])
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "black", lwd = lw)
  }
  
  #For Averaged EFD (Direction 2)
  x <- ave_mt[ave_mt$Direction =="b",-c(1:2)]
  x <- as.numeric(x[1,])
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "red", lwd = lw)
  
  #For raw EFD (Direction 1)
  for(r in 1:length(Replicate)){
    x <- mt[mt$Direction =="b",-c(1:3)]
    x <- as.numeric(x[r,])
    ef <- as.matrix(x)
    coord <- ef2coord(ef)
    plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
         ann = F, axes = F, col = "black", lwd = lw)
  }
}

rm(x, ef, coord, lw, i, Direction, r, ef2coord, mt, ave_mt, ID, Replicate)




################################3.Prediction_of_EFDs_of_F1_by_PPdelta_using_parental_raw_and_averaged_EFD_data#############################################################################################
#Predicted_F1_EFD_based_on_parental_raw_EFD_data_was_saved_in_"Pre_df"
#Predicted_F1_EFD_based_on_parental_averaged_EFD_data_was_saved_in_"Ave_Pre_df"

ID <- c("Parent1", "Parent2") #Define parental ID
Direction <- c("a", "b") #Define direction (a=Direction1, b=Direction2)
Replicate <- c("P1", "P2", "P3", "P4", "P5") #Define_fruit_replicates_(In this trial, we have five data for one direction per one accession)
Rep_R <- as.data.frame(read.csv("./Dataset/Representative_ratio_between_dominance_and_additive_effect.csv")) #Load_representative_ratio_between_dominance_and_additive_effect_derived_from_156_F1_accessions_in_author's_research

#Prepare_matrix_for_saving_predicted_EFDs
Pre_df <- as.data.frame(matrix(NA,
                           ncol = 82,
                           nrow = 0))
colnames(Pre_df) <- c("Direction", "Combination", colnames(df)[-c(1:3)])
Ave_Pre_df <- Pre_df
Pre_df2 <- Ave_Pre_df2 <- Pre_df[1,]
Pre_df2[1,] <- Ave_Pre_df2[1,] <- NA

#Perform_PPδ_for_Averaged_EFD
for(x in 1:length(Direction)){
  pheno <- Ave_df[Ave_df$Direction == Direction[x],]
      p1 <- ID[1] #Mother_parent_ID
      p2 <- ID[2] #Father_parent_ID
      p1h <- as.numeric(pheno[pheno$ID == p1,-c(1:2)]) #Mother_parent_EFD
      p2h <- as.numeric(pheno[pheno$ID == p2,-c(1:2)]) #Mother_parent_EFD
      
      midpoint <- (p1h + p2h)/2 #Midpoint_calculation
      a <- abs(p1h - p2h)/2 #Additive_effect_calculation
      d_est <- a*Rep_R[Rep_R$Direction == Direction[x],-c(1:2)]#Estimate_dominance_effects_in_the_crossing_combination
      d_est[,c(1,21,41)] <- 0
      pre_pheno <- as.numeric(midpoint + d_est) #Calculate_predicted_EFDs_for_the_F1_by_adding_estimated_dominance_effects_on_the_midpoint_EFDs
      Ave_Pre_df3 <- Ave_Pre_df2
      Ave_Pre_df3[1,] <- c(Direction[x],
                           "Parent1_X_Parent2",
                       pre_pheno)
      Ave_Pre_df <- rbind(Ave_Pre_df, Ave_Pre_df3)
}

Ave_Pre_df[,(1+2)] <- 1 #Put_1_on_a1 which should be constant values
Ave_Pre_df[,(21+2)] <- 0 #Put_0_on_b1 which should be constant values
Ave_Pre_df[,(41+2)] <- 0 #Put_0_on_c1 which should be constant values

for(i in 1:80){
  Ave_Pre_df[,c(i+2)] <- as.numeric(Ave_Pre_df[,c(i+2)])
}

rm(x, a, pheno, p1, p2, p1h, p2h, midpoint, d_est, pre_pheno, Ave_Pre_df2, Ave_Pre_df3)

#Perform_PPδ_for_raw_EFD
for(x in 1:length(Direction)){
  pheno_P1 <- df[df$Direction == Direction[x] & df$ID == "Parent1",]
  pheno_P2 <- df[df$Direction == Direction[x] & df$ID == "Parent2",]
  for(i in 1:nrow(pheno_P1)){
    for(j in i:nrow(pheno_P2)){
    p1 <- Replicate[i] #Mother_parent_ID
    p2 <- Replicate[j] #Father_parent_ID
    p1h <- as.numeric(pheno_P1[pheno_P1$Replicate == p1,-c(1:3)]) #Mother_parent_EFD
    p2h <- as.numeric(pheno_P2[pheno_P2$Replicate == p2,-c(1:3)]) #Mother_parent_EFD
    
    midpoint <- (p1h + p2h)/2 #Midpoint_calculation
    a <- abs(p1h - p2h)/2 #Additive_effect_calculation
    d_est <- a*Rep_R[Rep_R$Direction == Direction[x],-c(1:2)]#Estimate_dominance_effects_in_the_crossing_combination
    d_est[,c(1,21,41)] <- 0
    pre_pheno <- as.numeric(midpoint + d_est) #Calculate_predicted_EFDs_for_the_F1_by_adding_estimated_dominance_effects_on_the_midpoint_EFDs
    Pre_df3 <- Pre_df2
    Pre_df3[1,] <- c(Direction[x],
                     paste("Parent1_",
                           p1,
                           "X",
                           "Parent2_",
                           p2, sep = ""),
                            pre_pheno)
    Pre_df <- rbind(Pre_df, Pre_df3)
    }
  }
}



Pre_df[,(1+2)] <- 1 #Put_1_on_a1 which should be constant values
Pre_df[,(21+2)] <- 0 #Put_0_on_b1 which should be constant values
Pre_df[,(41+2)] <- 0 #Put_0_on_c1 which should be constant values

for(i in 1:80){
  Pre_df[,c(i+2)] <- as.numeric(Pre_df[,c(i+2)])
}


rm(d_est, pheno_P1, pheno_P2, Pre_df2, Pre_df3, a, Direction, i, ID, j, midpoint, p1, p1h, p2, p2h, pre_pheno, Replicate, x)




#################################3.Draw_fruit_contours_based_on_predicted_EFDs_of_F1######################################################
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
ID <- c("Parent1", "Parent2") #Define parental ID
Direction <- c("a", "b") #Define direction (a=Direction1, b=Direction2)
Replicate <- c("P1", "P2", "P3", "P4", "P5") #Define_fruit_replicates_(In this trial, we have five data for one direction per one accession)
lw <- 2 #line thickness
par(mfrow = c(2, 16), mar = c(0.0, 0.0, 0.0, 0.0)) #Determine_the_number_of_row(reft)_and_column(right)

#Draw_real_and_predicted_EFD-based fruit contours of the F1
#MEMO_in_plotted_figure:Averaged EFD-based contours were drawn as red lines & Raw EFD-based contours were drawn as black lines
#MEMO_in_plotted_figure:The contours in the first row showed Direction 1 and those in second row showed Direction 2
for(x in 1:length(Direction)){
  mt <- Pre_df[Pre_df$Direction == Direction[x],]
  ave_mt <- Ave_Pre_df[Ave_Pre_df$Direction == Direction[x],]
  
  x <- as.numeric(ave_mt[1,-c(1:2)])
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.0, 1.0), asp = 1,
       ann = F, axes = F, col = "red", lwd = lw)
  
  for(i in 1:nrow(mt)){
    x <- as.numeric(mt[i,-c(1:2)])
    ef <- as.matrix(x)
    coord <- ef2coord(ef)
    plot(coord$y, -coord$x, type = "l", xlim = c(-1.0, 1.0), asp = 1,
         ann = F, axes = F, col = "black", lwd = lw)
  }
  }
rm(coord, ef, mt, i, lw, x)
