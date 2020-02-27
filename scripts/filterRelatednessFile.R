#1) Load in Packages, take in command line arguments
#Packages

#Load in command line arguments
args <- commandArgs(trailingOnly = T)

#produce an error if 0 command line arguments entered
stopifnot(length(args) > 0)

#2) Read variables from args

#take in first argument as input file
input <- args[1]

#default values for threshold
threshold <- 0.1

#default value for output. Uses input string, chops off suffix, and adds on expected suffix for next step of snakefile. 
output <- paste0(strsplit(input,split = ".relatedness2")[[1]],"_relatedness2_filtered.txt")

#If second two arguments are specified, take in threshold and output respectively
if (length(args) > 1) {
  threshold <- as.double(args[2])
}
if (length(args) > 2) {
  output <- args[3]
}

#read in file specified by input. Both input and output include paths to their files from the folder where the snakefile is

relatedness <- read.table(input, sep = "\t", header=TRUE, stringsAsFactors = FALSE)

#3) Find unique combinations of Relatedness that pass the threshold

#create empty vectors that will get the values of significant relatedness chromosomes.
pass1 <- c()
pass2 <- c()

#Filters relatedness file based on threshold, and fills the two vectors from above. Does not consider duplicate sample combinations
for (i in 1:nrow(relatedness)) {
  if (relatedness[i,7][[1]] > threshold) {
    if (relatedness[i,1][[1]] != relatedness[i,2][[1]]) {
      #if this combination was previously stored in pass1 and 2, then the first of the combination will be in pass2, and the second will be in pass1
      #checking for the intersection of index1 and index2 will see if there is a common point where these were stored. If so, this combination was already recorded
      index1 <- grep(paste0("^",relatedness[i,1][[1]],"$"),pass2)
      index2 <- grep(paste0("^",relatedness[i,2][[1]],"$"),pass1)
      
      #if index1 or index2 are empty, this also passes. Pass1 is the first number of each combo. Pass2 is the second. 
      if (length(intersect(index1,index2)) == 0 | length(index1) == 0 | length(index2) == 0) {
        pass1 <- c(pass1,relatedness[i,1][[1]])
        pass2 <- c(pass2,relatedness[i,2][[1]])
      }
    }
  }
}

#4) Remove samples that are related to the most number of other samples. Repeat is a looped search to remove these.
removal <- c()
rand = function(d) sample(1:d,1,replace=T)
repeat{
  
  #build a table of # of occurances of each value
  occurances <- table(c(pass1,pass2))
  
  #if occurances is empty, break
  if (length(occurances) == 0) {
    break
  }
  
  #find the most frequently occuring individual in relatedness and add them to a vector to be removed
  maxNum <- grep(max(occurances),occurances)
  maxVal <- names(occurances)[maxNum[rand(length(maxNum))]]
  
  removal <- c(removal, maxVal)
  
  #remove the individual from pass1 and pass2 and reconsider the occurances
  index1 <- grep(paste0("^",maxVal,"$"),pass1)
  index2 <- grep(paste0("^",maxVal,"$"),pass2)
  indexAll <- c(index1,index2)
  pass1 <- pass1[-indexAll]
  pass2 <- pass2[-indexAll]
}

#5) Write final output
write.table(removal, file = output, quote = F, col.names = F, row.names = F, sep = "\t")
