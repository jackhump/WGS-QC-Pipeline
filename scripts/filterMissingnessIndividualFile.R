#1) Load in Packages, take in command line arguments

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
missingnessINDI <- read.table(input, header=TRUE, sep = "\t", stringsAsFactors=FALSE)

missing_samples <- missingnessINDI[ missingnessINDI$F_MISS > threshold, "INDV"]

writeLines(missing_samples, con = output)
