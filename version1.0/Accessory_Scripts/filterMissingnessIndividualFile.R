#1) Load in Packages, take in command line arguments
#Packages
library(readr)

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
output <- paste0(strsplit(input,split = ".imiss")[[1]],"_imiss_filtered.txt")

#If second two arguments are specified, take in threshold and output respectively
if (length(args) > 1) {
  threshold <- args[2]
}
if (length(args) > 2) {
  output <- args[3]
}

#read in file specified by input. Both input and output include paths to their files from the folder where the snakefile is
missingnessINDI <- read_delim(input, "\t", escape_double = FALSE, trim_ws = TRUE)

#3) filter missingness file

#create empty vector that will get the values of significant missingness chromosomes.
indi <- c()

#Filters missingness file based on threshold, and fills the two vectors from above
for (i in 1:nrow(missingnessINDI)) {
  if (missingnessINDI[i,5][[1]] > threshold) {
    indi <- c(indi,missingnessINDI[i,1][[1]])
  }
}

#4) Write final output
write.table(indi, file = output, quote = F, col.names = F, row.names = F, sep = "\t")
