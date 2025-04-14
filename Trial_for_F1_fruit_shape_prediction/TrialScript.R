########## Trial script: Perform fruit shape prediction using PPδ ######################################################################
###################### NOTE!!!################################################################################################################################
# 1.This trial uses EFD data of the two accessions (Parent1 and Parent2) as parental examples and predicts their F1 fruit shape contours
# 2.The EFD data was obtained by SHAPE program saved as ".nef" format (The number of harmonics was 20)
# 3.In the shape program, EFDs for one direction (Direction 1 or Direction 2) were extracted from five fruits per one accession and saved as same ".nef" file
# 4.For implementation of PPδ, representative ratio between dominance and additive_effect derived from the 156 F1 accession in author's research.
# 5. fruit shape prediction was performed based on parental raw EFD and averaged EFD data

################################# Step.1 Data loading from ".nef" format file and Data formatting ###################################################################################################################
# Parental raw EFD data was saved in "df"
# Parental Averaged EFD data was saved in "Ave df"
setwd(getwd())
ID <- c("Parent1", "Parent2") # Define parental ID
Direction <- c("a", "b") # Define direction (a=Direction1, b=Direction2)
Replicate <- c("P1", "P2", "P3", "P4", "P5") # Define fruit replicates (In this trial, we have five data for one direction per one accession)

# Load EFD data
Parent1_a <- as.data.frame(read.delim("./Dataset/Parent1_a.nef", sep = "")) # Load EFD data (Direction 1) from five fruits of Parent1 using SHAPE program
Parent1_b <- as.data.frame(read.delim("./Dataset/Parent1_b.nef", sep = "")) # Load EFD data (Direction 2) from five_fruits of Parent1 using SHAPE program

Parent2_a <- as.data.frame(read.delim("./Dataset/Parent2_a.nef", sep = "")) # Load_EFD_data(Direction 1) from five fruits of Parent2 using SHAPE program
Parent2_b <- as.data.frame(read.delim("./Dataset/Parent2_b.nef", sep = "")) # Load_EFD_data(Direction 2) from_five_fruits_of Parent2 using SHAPE program

# Data collection formatting for raw EFD data
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


# Data collection of formatting for averaged EFD data
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

################################# Step 2. Draw parental fruit contours based on the raw and averaged EFD #####################################################
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

ID <- c("Parent1", "Parent2") # Define parental ID
Direction <- c("a", "b") # Define direction (a=Direction1, b=Direction2)
Replicate <- c("P1", "P2", "P3", "P4", "P5") # Define fruit replicates (In this trial, we have five data for one direction per one accession)

# Setting for illustration
lw <- 2 # Line thickness
par(mfrow = c(length(ID), 12), mar = c(0.0, 0.0, 0.0, 0.0)) # Determine the number of row (left) and column (right)

# Draw raw and averaged EFD-based fruit contours of the two parents (Parent1 and Parent2)
# MEMO in plotted figure: Averaged EFD-based contours were drawn as red lines & Raw EFD-based contours were drawn as black lines
# MEMO in plotted figure: The contours in the first row showed parent1 and those in second row showed parent2
# MEMO in plotted figure: For raw EFD-based contours, left five contours showed direction1 and right five contours showed direction2.

for(i in 1:length(ID)){
  mt <- df[df$ID == ID[i],]
  ave_mt <- Ave_df[Ave_df$ID == ID[i],]
  
  # For Averaged EFD (Direction 1)
  x <- ave_mt[ave_mt$Direction =="a",-c(1:2)]
  x <- as.numeric(x[1,])
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "red", lwd = lw)
  
  # For raw EFD (Direction 1)
  for(r in 1:length(Replicate)){
  x <- mt[mt$Direction =="a",-c(1:3)]
  x <- as.numeric(x[r,])
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "black", lwd = lw)
  }
  
  # For Averaged EFD (Direction 2)
  x <- ave_mt[ave_mt$Direction =="b",-c(1:2)]
  x <- as.numeric(x[1,])
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "red", lwd = lw)
  
  # For raw EFD (Direction 1)
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


################################ Step 3. Prediction of EFDs of F1 by PPδ using parental raw and averaged EFD data #############################################################################################
# Predicted F1 EFD based on parental raw EFD data was saved in "Pre_df"
# Predicted F1 EFD based on parental averaged EFD data was saved in "Ave_Pre_df"

ID <- c("Parent1", "Parent2") # Define parental ID
Direction <- c("a", "b") # Define direction (a=Direction1, b=Direction2)
Replicate <- c("P1", "P2", "P3", "P4", "P5") # Define fruit replicates (In this trial, we have five data for one direction per one accession)
Rep_R <- as.data.frame(read.csv("./Dataset/Representative_ratio_between_dominance_and_additive_effect.csv")) # Load representative ratio between dominance and additive effect derived from 156 F1 accessions in author's research

# Prepare matrix for saving predicted EFDs
Pre_df <- as.data.frame(matrix(NA,
                           ncol = 82,
                           nrow = 0))
colnames(Pre_df) <- c("Direction", "Combination", colnames(df)[-c(1:3)])
Ave_Pre_df <- Pre_df
Pre_df2 <- Ave_Pre_df2 <- Pre_df[1,]
Pre_df2[1,] <- Ave_Pre_df2[1,] <- NA

# Perform_PPδ for Averaged EFD
for(x in 1:length(Direction)){
  pheno <- Ave_df[Ave_df$Direction == Direction[x],]
      p1 <- ID[1] #Mother_parent_ID
      p2 <- ID[2] #Father_parent_ID
      p1h <- as.numeric(pheno[pheno$ID == p1,-c(1:2)]) # Mother parent EFD
      p2h <- as.numeric(pheno[pheno$ID == p2,-c(1:2)]) # Mother parent EFD
      
      midpoint <- (p1h + p2h)/2 # Midpoint calculation
      a <- abs(p1h - p2h)/2 # Additive effect calculation
      d_est <- a*Rep_R[Rep_R$Direction == Direction[x],-c(1:2)] # Estimate dominance effects in the crossing combination
      d_est[,c(1,21,41)] <- 0
      pre_pheno <- as.numeric(midpoint + d_est) # Calculate predicted EFDs for the F1 by adding estimated dominance effects on the midpoint EFDs
      Ave_Pre_df3 <- Ave_Pre_df2
      Ave_Pre_df3[1,] <- c(Direction[x],
                           "Parent1_X_Parent2",
                       pre_pheno)
      Ave_Pre_df <- rbind(Ave_Pre_df, Ave_Pre_df3)
}

Ave_Pre_df[,(1+2)] <- 1 # Put 1 on_a1 which should be constant values
Ave_Pre_df[,(21+2)] <- 0 # Put 0 on b1 which should be constant values
Ave_Pre_df[,(41+2)] <- 0 # Put 0 on c1 which should be constant values

for(i in 1:80){
  Ave_Pre_df[,c(i+2)] <- as.numeric(Ave_Pre_df[,c(i+2)])
}

rm(x, a, pheno, p1, p2, p1h, p2h, midpoint, d_est, pre_pheno, Ave_Pre_df2, Ave_Pre_df3)

# Perform PPδ for raw EFD
for(x in 1:length(Direction)){
  pheno_P1 <- df[df$Direction == Direction[x] & df$ID == "Parent1",]
  pheno_P2 <- df[df$Direction == Direction[x] & df$ID == "Parent2",]
  for(i in 1:nrow(pheno_P1)){
    for(j in i:nrow(pheno_P2)){
    p1 <- Replicate[i] # Mother parent ID
    p2 <- Replicate[j] # Father parent ID
    p1h <- as.numeric(pheno_P1[pheno_P1$Replicate == p1,-c(1:3)]) # Mother parent EFD
    p2h <- as.numeric(pheno_P2[pheno_P2$Replicate == p2,-c(1:3)]) # Mother parent EFD
    
    midpoint <- (p1h + p2h)/2 # Midpoint calculation
    a <- abs(p1h - p2h)/2 # Additive effect calculation
    d_est <- a*Rep_R[Rep_R$Direction == Direction[x],-c(1:2)] # Estimate dominance effects in the crossing combination
    d_est[,c(1,21,41)] <- 0
    pre_pheno <- as.numeric(midpoint + d_est) # Calculate predicted EFDs for the F1 by adding estimated dominance effects on the midpoint EFDs
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

Pre_df[,(1+2)] <- 1 # Put 1 on_a1 which should be constant values
Pre_df[,(21+2)] <- 0 # Put 0 on_b1 which should be constant values
Pre_df[,(41+2)] <- 0 # Put 0 on_c1 which should be constant values

for(i in 1:80){
  Pre_df[,c(i+2)] <- as.numeric(Pre_df[,c(i+2)])
}

rm(d_est, pheno_P1, pheno_P2, Pre_df2, Pre_df3, a, Direction, i, ID, j, midpoint, p1, p1h, p2, p2h, pre_pheno, Replicate, x)


################################# Step 3. Draw fruit contours based on predicted EFDs of F1 ######################################################
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
ID <- c("Parent1", "Parent2") # Define parental ID
Direction <- c("a", "b") # Define direction (a=Direction1, b=Direction2)
Replicate <- c("P1", "P2", "P3", "P4", "P5") # Define fruit replicates (In this trial, we have five data for one direction per one accession)
lw <- 2 # Line thickness
par(mfrow = c(2, 16), mar = c(0.0, 0.0, 0.0, 0.0)) # Determine the number of row (left) and column (right)

# Draw real and predicted EFD-based fruit contours of the F1
# MEMO in plotted figure: Averaged EFD-based contours were drawn as red lines & Raw EFD-based contours were drawn as black lines
# MEMO in plotted figure: The contours in the first row showed Direction 1 and those in second row showed Direction 2
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
