

library(readr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]

DP <- as.numeric(readLines( paste0(input, ".DP.txt" ) ))

MQ <- as.numeric(readLines( paste0(input, ".MQ.txt" ) ))

VQSLOD <- as.numeric(readLines( paste0(input, ".VQSLOD.txt" ) ))

# calculate median and sd of DP - this depends on sample size rather than quality
DP_med <- median(DP)
DP_sd <- sd(DP)


# plot histogram of DP with range 0 - median(DP) + 2*sd(DP)

pdf( paste0( input, "_vcf_metrics.pdf") )
hist(DP,500)
hist(DP[ DP < (DP_med + 2*DP_sd) ], 500, main = "DP 0 - median + 2 sd" )
hist(MQ,60)
hist(VQSLOD)
hist(VQSLOD[VQSLOD > 0 & VQSLOD < 40],45)

sink( paste0( input, "_vcf_metrics.txt") )
print("* DP ")
print(" % > median - 2 sd ")
length(DP[DP> DP_med - 2*DP_sd ])/length(DP)
print("% > median")
length(DP[DP>DP_med])/length(DP)
print("% > median + 2 sd")
length(DP[DP> DP_med + 2*DP_sd])/length(DP)

#hist(MQ,60)
print("* MQ")
print(" % > 58.75" )
length(MQ[MQ > 58.75])/length(MQ)
print(" % > 50" )
length(MQ[MQ > 50])/length(MQ)

#hist(VQSLOD)
#hist(VQSLOD[VQSLOD > 0],45)
print("* VQSLOD")
print(" % > 7.81")
length(VQSLOD[VQSLOD > 7.81])/length(VQSLOD)
print(" % > 2 ")
length(VQSLOD[VQSLOD > 2])/length(VQSLOD)

