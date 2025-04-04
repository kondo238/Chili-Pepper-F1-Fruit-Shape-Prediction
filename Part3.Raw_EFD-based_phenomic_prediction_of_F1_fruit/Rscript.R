##########Part 3.Prediction_of_EFDs_for_the_159_F1_accessions_by_PPδ_based_on_raw_EFD_data_of_the_parental_accessions######################################################################
setwd(getwd())
library(stringr)
library(progress)


#################################1.Data_loading###################################################################################################################
Acclist <- as.data.frame(read.csv("./Dataset/Accession_list.csv")) #Load_accession_list_data
All_data <- as.data.frame(read.csv("./Dataset/Raw_EFD_data.csv"))# Load_raw_EFD_data_(80 variables)
Inbred_ID <- Acclist[Acclist$PopulationType == "Inbred",3]
F1parent_ID <- Acclist[which(str_detect(Acclist$Note, pattern = "F1 parent")),3]
combi <- as.data.frame(read.csv("./Dataset/Parental_combinations_of_F1.csv")) #Load_list_for_parental_combinations_of_F1
Rep_R <- as.data.frame(read.csv("./Dataset/Representative_ratio_between_dominance_and_additive_effect.csv")) #Load_representative_ratio_between_dominance_and_additive_effect



################################2.Prediction_of_EFDs_of_F1_by_PPdelta#############################################################################################
#Setting
Direction <- c("a", "b") #Setting direction
Years <- c("2021", "2022", "2023") #Setting cultivation year 
Replicate <- c("P1", "P2", "P3", "P4", "P5") #Setting fruit repricates
F1_ID <- Acclist[Acclist$PopulationType == "F1",3] #Extract F1 ID

#Prepare_matrix_for_saving_predicted_EFDs
df <- as.data.frame(matrix(NA,
                           ncol = 82,
                           nrow = length(Direction)*length(F1_ID)*(length(Years)*length(Replicate))^2))
df[,1] <- c(rep(Direction[1], length(F1_ID)*(length(Years)*length(Replicate))^2),
            rep(Direction[2], length(F1_ID)*(length(Years)*length(Replicate))^2))
colnames(df) <- c("Direction", "ID", colnames(All_data)[-c(1:4)])

g <- vector()
for(i in 1:length(F1_ID)){
  a <- rep(F1_ID[i], (length(Replicate)*length(Years))^2)
  g <- c(g, a)
}
g <- rep(g, length(Direction))
df[,2] <- g


#Perform_PPδ
for(x in 1:length(Direction)){
  pheno <- All_data[All_data$Direction == Direction[x],]
  pb <- progress_bar$new(total = length(F1_ID)) #Set progress bar
  for(i in 1:length(F1_ID)){
    pb$tick()
    p1 <- as.character(combi[i,1]) #Mother_parent_ID
    p2 <- as.character(combi[i,2]) #Father_parent_ID
    f1 <- as.character(combi[i,3]) #F1_ID
    p1h <- pheno[pheno$ID == p1,] #Mother_parent_EFD
    p2h <- pheno[pheno$ID == p2,] #Father_parent_EFD
    
    pre_F1 <- as.data.frame(df[0,])
    for(r in 1:(length(Years)*length(Replicate))){
      p1h_2 <- p1h[r,-c(1:4)]
      for(q in 1:(length(Years)*length(Replicate))){
        p2h_2 <- p2h[q,-c(1:4)]
        midpoint <- (p1h_2 + p2h_2)/2 #Midpoint_calculation
        a <- abs(p1h_2 - p2h_2)/2 #Additive_effect_calculation
        d_est <- a*Rep_R[Rep_R$Direction == Direction[x],-c(1:2)] #Estimate_dominance_effects_in_the_crossing_combination
        pre_pheno <- midpoint + d_est #Calculate_predicted_EFDs_for_the_F1_by_adding_estimated_dominance_effects_on_the_midpoint_EFDs
        pre_F1_2 <- cbind(Direction=(Direction[x]),
                        ID=F1_ID[i],
                        pre_pheno)
        pre_F1 <- rbind(pre_F1, pre_F1_2)
      }
    }
      
      df[df$Direction == Direction[x] &
           df$ID == F1_ID[i],] <- pre_F1
      
      Sys.sleep(1 / i) #Update progress bar
  }
}

df[,(1+2)] <- 1
df[,(21+2)] <- 0
df[,(41+2)] <- 0


PPdelta_predicted_EFD <- df
write.csv(PPdelta_predicted_EFD,  "Predicted_F1_EFDs_with_all_combination_with_dominance_effects.csv") #Save_predicted_EFDs_as_csv_format
rm(df, mean, mean_mt, ph1, ph1_2, ph2, ph2_2, pheno, a, Direction, f1, g, i, no2, p1, p2, q, r, Replicate, x, Years)






#################################3.Draw_fruit_contours_based_on_predicted_and_real_EFDs_of_F1######################################################
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
Direction <- "a" #Determine_the_direction (a=Direction1 & b=Direction2)
par(mfrow = c(16, 15), mar = c(0.0, 0.0, 0.0, 0.0)) #Determine_the_number_of_row(reft)_and_column(right)
F1 <- F1_ID[1] #Put_one_F1_ID(1~159)
Pre_df <- PPdelta_predicted_EFD[PPdelta_predicted_EFD$ID == F1,] #Extract_predicted_EFDs_of_the_choosed_F1
Real_df <- All_data[All_data$ID == F1,] #Extract_real_raw_EFDs_of_the_choosed_F1

#Draw_real_and_predicted_EFD-based fruit contours of the F1
#Draw_Real_EFD-based_fruit_contours(Black)
for(i in 1:(length(Years)*length(Replicate))){
  df <- Real_df[Real_df$Direction == Direction,]
  x <- df[i,-c(1:4)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "black", lwd = lw)
}

#Draw_Real_EFD-based_fruit_contours(Black)
for(i in 1:(nrow(Pre_df)/2)){
  df <- Pre_df[Pre_df$Direction == Direction,]
  x <- df[i,-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "red", lwd = lw)
}

rm(x, ef, coord, lw, i, Direction)
