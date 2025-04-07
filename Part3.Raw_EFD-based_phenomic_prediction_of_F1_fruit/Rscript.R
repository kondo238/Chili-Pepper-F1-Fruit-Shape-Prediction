########## Part 3. Part 3. Prediction of EFDs for the 159 F1 accessions using PPδ based on raw EFD data of parental accessions ######################################################################
# Set working directory
setwd(getwd())
# Load required libraries
library(stringr)
library(progress)


################################# 1 Data loading ###################################################################################################################
# Load list of accessions
Acclist <- as.data.frame(read.csv("./Dataset/Accession_list.csv")) 

# Load raw Elliptical Fourier Descriptor (EFD) data for all accessions (80 variables)
All_data <- as.data.frame(read.csv("./Dataset/Raw_EFD_data.csv"))

# Extract IDs for inbred lines
Inbred_ID <- Acclist[Acclist$PopulationType == "Inbred",3]

# Extract IDs of parental lines used in F1 crosses
F1parent_ID <- Acclist[which(str_detect(Acclist$Note, pattern = "F1 parent")),3]

# Load list of parental combinations used to generate F1s (columns: mother, father, F1)
combi <- as.data.frame(read.csv("./Dataset/Parental_combinations_of_F1.csv")) 

# Load ratio between dominance and additive effects, estimated previously
Rep_R <- as.data.frame(read.csv("./Dataset/Representative_ratio_between_dominance_and_additive_effect.csv")) 

################################ 2. Prediction of F1 EFDs using PPδ #############################################################################################

# Define directions (a = horizontal, b = vertical or vice versa)
Direction <- c("a", "b") 

# Define cultivation years and replicate labels
Years <- c("2021", "2022", "2023")
Replicate <- c("P1", "P2", "P3", "P4", "P5")

# Extract F1 accession IDs
F1_ID <- Acclist[Acclist$PopulationType == "F1",3] #Extract F1 ID

# Prepare a data frame to store predicted EFDs (82 columns: Direction + ID + 80 EFDs)
df <- as.data.frame(matrix(NA,
                           ncol = 82,
                           nrow = length(Direction)*length(F1_ID)*(length(Years)*length(Replicate))^2))
df[,1] <- c(rep(Direction[1], length(F1_ID)*(length(Years)*length(Replicate))^2),
            rep(Direction[2], length(F1_ID)*(length(Years)*length(Replicate))^2))
colnames(df) <- c("Direction", "ID", colnames(All_data)[-c(1:4)])

# Fill ID column in the prediction dataframe
g <- vector()
for(i in 1:length(F1_ID)){
  a <- rep(F1_ID[i], (length(Replicate)*length(Years))^2)
  g <- c(g, a)
}
g <- rep(g, length(Direction))
df[,2] <- g

# Main prediction loop using PPδ
for(x in 1:length(Direction)){
  pheno <- All_data[All_data$Direction == Direction[x],] # Filter by direction
  pb <- progress_bar$new(total = length(F1_ID)) # Set progress bar
  for(i in 1:length(F1_ID)){
    pb$tick()
    
    # Get parent and F1 IDs from combination list
    p1 <- as.character(combi[i,1]) # Mother
    p2 <- as.character(combi[i,2]) # Father
    f1 <- as.character(combi[i,3]) # F1
    # Get EFDs of each parent
    p1h <- pheno[pheno$ID == p1,] # Mother
    p2h <- pheno[pheno$ID == p2,] # Father

    # Initialize temporary data frame for current F1 predictions
    pre_F1 <- as.data.frame(df[0,])

    # Loop through all combinations of p1 and p2 fruit samples (years x replicates)
    for(r in 1:(length(Years)*length(Replicate))){
      p1h_2 <- p1h[r,-c(1:4)] # Remove metadata
      for(q in 1:(length(Years)*length(Replicate))){
        p2h_2 <- p2h[q,-c(1:4)] # Remove metadata

        # Midpoint of the two parents (purely additive expectation)
        midpoint <- (p1h_2 + p2h_2)/2

        # Additive effect: half the absolute difference
        a <- abs(p1h_2 - p2h_2)/2 
        # Dominance effects in the crossing combination
        d_est <- a*Rep_R[Rep_R$Direction == Direction[x],-c(1:2)] 

        # Final predicted F1 phenotype (midpoint + dominance)
        pre_pheno <- midpoint + d_est 

        # Assemble one prediction row
        pre_F1_2 <- cbind(Direction=(Direction[x]),
                        ID=F1_ID[i],
                        pre_pheno)
        pre_F1 <- rbind(pre_F1, pre_F1_2)
      }
    }
      # Insert the predicted EFDs into the global df
      df[df$Direction == Direction[x] &
           df$ID == F1_ID[i],] <- pre_F1
      
      Sys.sleep(1 / i) # Update progress bar
  }
}

# Post-processing: force some coefficients to fixed values (corresponding to translation/scale)
df[,(1+2)] <- 1
df[,(21+2)] <- 0
df[,(41+2)] <- 0

# Save final predicted EFD matrix
PPdelta_predicted_EFD <- df
write.csv(PPdelta_predicted_EFD,  "Predicted_F1_EFDs_with_all_combination_with_dominance_effects.csv") 

rm(df, p1h, p1h_2, p2h, p2h_2, pheno, a, Direction, f1, g, i, p1, p2, q, r, x, d_est, midpoint, pb, pre_F1, pre_F1_2, pre_pheno) # Clean up the workspace






################################# 3. Draw fruit contours based on predicted and real EFDs of F1 ######################################################
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


# Set plot layout and parameters
lw <- 2 # Line thickness
Direction <- "a" # Choose direction to visualize (a=Direction1 & b=Direction2)
par(mfrow = c(16, 15), mar = c(0.0, 0.0, 0.0, 0.0)) # Plot grid layout (rows x columns)

# Choose one F1 accession to visualize
F1 <- F1_ID[1] # Put one F1 ID (1~159)

# Extract real and predicted EFDs for selected F1
Pre_df <- PPdelta_predicted_EFD[PPdelta_predicted_EFD$ID == F1,]
Real_df <- All_data[All_data$ID == F1,] 

# Plot real fruit shapes (in black)
for(i in 1:(length(Years)*length(Replicate))){
  df <- Real_df[Real_df$Direction == Direction,]
  x <- df[i,-c(1:4)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "black", lwd = lw)
}

# Plot predicted fruit shapes (in red)
for(i in 1:(nrow(Pre_df)/2)){
  df <- Pre_df[Pre_df$Direction == Direction,]
  x <- df[i,-c(1:2)]
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "red", lwd = lw)
}

rm(x, ef, coord, lw, i, Direction, F1, Pre_df, Real_df, df) # Clean up the workspace
