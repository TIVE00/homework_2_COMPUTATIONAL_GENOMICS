#raw_count.txt is a tab delimited file that contains raw data 
#and patients are divided in 2 groups: 1) asymptomatic, 2) symptomatic

#gene_annot.txt provides additional info on different genes, including their length

## FIRST FUNCTION##

#input: #exprData = numeric data matrix w\ same format as raw_count (gene IDs = rows, subject IDs = columns)
        #pdffilename = name of .pdf file where to save plots
        #pcolor = color used to draw points
        #lcolor = color used to draw line
#output: .pdf file w\ MvA plots of each sample VS. sample 1 (one plot per page). For each plot:
         #provide meaningful names for x-axis and y-axis and title
         #draw horizontal line corresponding to x-axis (line color must be lcolor)
         #tune the size of the plot and size/type of the points to make plot readable

MvAplot <- function ( exprData, pdffilename, pcolor = 'black', lcolor = 'red'){
  
  # Read exprData as a data frame
  exprData <- read.table(exprData, header = TRUE, row.names = 1, sep = "\t")
  
  # Open a PDFto save plots
  pdf(pdffilename, width = 6, height = 4) 
  
  samples <- ncol(exprData) #number of columns in exprData matrix
  
  # Loop through each sample to create MA plots
  for (sample in 2:samples) { #start from sample 2 since sample 1 is the reference sample 
    
    # Calculate M and A values
    M <- log2(exprData[, sample] +1) - log2(exprData[, 1]) # M = log2 mean #add 1 in both terms in order to exclude cases of log(0)
    A <- (log2(exprData[, sample] +1) + log2(exprData[, 1] +1 )) / 2  #  A = log2 average intensity #add 1 in both terms in order to exclude cases of log(0)
    
    # Print each plot on a new page in the PDF
    plot(A, M, pch = 16, col = pcolor, xlab = 'A = (log2(r*1) + log2(r*2))/2', ylab ='M = log2(r*1) - log2(r*2)', main = paste("MvA Plot - Sample", colnames(exprData)[sample], "vs Sample", colnames(exprData)[1]))
    
    # Horizontal line corresponding to x-axis
    abline(h = 0, col = lcolor, lwd = 1.5, lty = 2)
  }
  
  # Close the PDF device
  dev.off()
  
}

#test
#MvAplot('/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 2/raw_count.txt','plot_MVA',pcolor='black',lcolor='red')

##SECOND FUNCTION##

#input:#exprData = numeric data matrix w\ same format as raw_count (gene IDs = rows, subject IDs = columns)
       #annot = data frame that provides additional info on different genes, including their length
       #Mtrim = number between 0 and 1 indicating fraction of observations to be deleted from each end (positive and negative values) of the sorted vector M before calculating the mean (Ms are the log ratios defined as in the MvA plot)
       #Atrim = vector of 2 elements indicating lower and upper thresholds to trim the most extreme values of A ( A is the average in log2 scale)
# what it does: #scales data by their sequencing depth and multiply by 10^6
                #calculates scaling factors SF (w\ respect to sample 1) by trimming the most extreme values of A and taking the trimmed means of M values ( use R function mean)
                #normalize data by their scaling factor w\ respect to sample 1
                #scales genes by their length and multiply by 10^3
#output: #list of 2 elements: i) normalized matrix (in original scale, not log) , ii) vector of scaling factors (in original scale, not log)

TMMnorm <- function (exprData, annot, Mtrim =0.2, Atrim = c(0,8)){
  
  # Read exprData as data frame
  exprData <- read.table(exprData, header = TRUE, row.names = 1, sep = "\t")
  
  # Read annot as data frame 
  annot <- read.table(annot, header = TRUE, row.names = 1, sep = "\t")
  
  gene_lengths <- annot[ ,ncol(annot)] #column of gene lengths in annot is always the last column
  
  # Scale data by sequencing depth and multiply by 10^6 ( = calculate CPM )
  scaled_data <- t(t(exprData)/colSums(exprData))*(10^6)
  
  samples <- ncol(scaled_data) #number of columns in scaled_data matrix
  
  # Initialize storage for M and A values
  M_values <- list()
  A_values <- list()
  
  # Loop through each sample to calculate M and A values
  for (sample in 2:ncol(scaled_data)) { #start from sample 2 since sample 1 is the reference sample 
    
    M <- log2(scaled_data[, sample] +1) - log2(scaled_data[, 1]) # M = log2 mean #add 1 in both terms in order to exclude cases of log(0)
    A <- (log2(scaled_data[, sample] +1) + log2(scaled_data[, 1] +1 )) / 2  #  A = log2 average intensity #add 1 in both terms in order to exclude cases of log(0)
    
    M_values[[sample]] <- M
    A_values[[sample]] <- A
  }
  
  # Combine all M and A values across samples in order to trim them
  all_M <- unlist(M_values)
  all_A <- unlist(A_values)
  
  # Trim A values : trheshold Atrim
  A_trimmed <- (all_A >= Atrim[1]) & (all_A <= Atrim[2])
  
  #keep only M values within the accepted A values
  M_to_sort <- all_M[A_trimmed]
  
  # Trim M values : threshold Mtrim
  M_sorted <- sort(M_to_sort)
  M_lower <- floor(Mtrim*length(M_sorted)) #lower bound
  M_upper <- ceiling ((1-Mtrim)*length(M_sorted)) #upper bound
  M_trimmed <- M_sorted[(M_lower +1):M_upper]
  
  # Initialize a vector to store scaling factors
  scaling_factors <- numeric(ncol(scaled_data))
  scaling_factors[1] <- 1  # reference sample for SF is sample 1
  
  # Loop trhough each sample to calculate SFs
  for(sample in 2:ncol(scaled_data)){
    
    sample_M <- M_values[[sample]]
    sample_M <- sample_M[A_trimmed]
    
    scaling_factors[sample] <- 2 ^ mean(sample_M, na.rm = TRUE)
  }
  
  # Normalize each sample by its SF with respect to sample 1
  normalized_data <- scaled_data / scaling_factors
  
  # Adjust by gene length and multiply by 10^3
  normalized_data <- sweep(normalized_data, 1, gene_lengths / (10^3), FUN = "/")
  
  # Return normalized data and scaling factors (both in original scale)
  list(normalized_matrix = normalized_data, scaling_factors = scaling_factors)
  
}
#test
#TMMnorm('/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 2/raw_count.txt','/Users/aurorativeron/Desktop/COMPUTER ENGINEERING/SECOND YEAR/FIRST SEMESTER/COMPUTATIONAL GENOMICS/homework 2/gene_annot.txt',0.2,c(2,8))
