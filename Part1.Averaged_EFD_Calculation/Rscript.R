##########################Analysis.1_Calculate_Averaged_EFD_data##############################################################################
# 1.Data loading
setwd(getwd())

Acclist <- as.data.frame(read.csv("./Dataset/Accession_list.csv")) # Load accession list data
All_data <- as.data.frame(read.csv("./Dataset/Raw_EFD_data.csv"))# Load raw EFD data (80 variables)

# 2.Calculate Averaged EFD values
# Define parameters
Direction <- c("a", "b") # a=Direction1 & b=Direction2
Year <- c(2021, 2022, 2023) # Not used in this step
ID <- Acclist$ID # Extract unique IDs

#Create an empty dataframe for average EFD values
aveEFD <- as.data.frame(matrix(NA,
                               ncol = 82, # 2 columns (Direction, ID) + 80 EFD variables
                               nrow = length(ID)*length(Direction)))

#Assign Direction and ID columns
aveEFD[,1] <- c(rep(Direction[1], length(ID)), # First half: Direction "a"
                rep(Direction[2], length(ID))) # Second half: Direction "b"
aveEFD[,2] <- rep(ID, length(Direction)) # Repeat IDs for both directions

# Define column names: Direction, ID, and 80 EFD variables (a1-a20, b1-b20, c1-c20, d1-d20)
colnames(aveEFD) <- c("Direction",
                      "ID",
                      paste("a", as.character(1:20), sep = ""),
                      paste("b", as.character(1:20), sep = ""),
                      paste("c", as.character(1:20), sep = ""),
                      paste("d", as.character(1:20), sep = ""))

# Compute average EFD values per ID and Direction
for(d in 1:length(Direction)){
  df <- All_data[All_data$Direction == Direction[d],]
for(i in 1:length(ID)){
  df2 <- df[df$ID == ID[i],]
  for(n in 1:80){
    aveEFD[aveEFD$Direction == Direction[d] &
                 aveEFD$ID == ID[i],(2+n)] <- mean(df2[,(4+n)])
  }
}
}

# Remove temporary variables 
rm(df, df2, d, i, n)

# Save the averaged EFD data as a CSV file
write.csv(aveEFD, "Averaged_EFD_data.csv", row.names = F) #Save_Averaged_EFD_data


# 3.Draw_fruit_contours_based_on_raw_and_averaged_EFDs
# Function to convert EFD coefficients to coordinates
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
par(mfrow = c(10, 30), mar = c(0.0, 0.0, 0.0, 0.0)) # Determine row number and column

# Draw averaged fruit contours (red) on raw contours (gray)
for(i in 1:length(ID)){

  # Extract raw EFD data for the given IA and Direction 
  raw_df <- All_data[All_data$Direction == Direction &
                       All_data$ID == ID[i],-c(1:4)]
  
  # Draw raw contours (grey)
  x <- as.matrix(raw_df[1,])
  ef <- as.matrix(x)
  coord <- ef2coord(ef)
  plot(coord$y, -coord$x, type = "l", xlim = c(-1.0, 1.0), ylim = c(-1.2, 1.2), asp = 1,
       ann = F, axes = F, col = "gray", lwd = lw)
  for(r in 2:nrow(raw_df)){
    x <- as.matrix(raw_df[r,])
    ef <- as.matrix(x)
    coord <- ef2coord(ef)
    lines(coord$y, -coord$x, col = "gray", lwd = lw)
  }
  
  # Draw averaged contours (red)
  x <- aveEFD[aveEFD$Direction == Direction &
                     aveEFD$ID == ID[i],-c(1:2)]
  ef <- as.matrix(x[1,])
  coord <- ef2coord(ef)
  lines(coord$y, -coord$x, col = "red", lwd = lw)
  }
  
# Remove temporary variables
rm(coord, raw_df, x, i, r, lw, ef)
