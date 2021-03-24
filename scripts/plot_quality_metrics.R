library(readr)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
cov <- args[2] 

DP <- as.numeric(readLines( paste0(input, ".DP.txt") ))

MQ <- as.numeric(readLines( paste0(input, ".MQ.txt") ))

VQSLOD <- as.numeric(readLines( paste0(input, ".VQSLOD.txt") ))

SAMPLES <- readLines( paste0(input, ".samples.txt") )

# calculate median and sd of DP - this depends on sample size rather than quality
DP_med <- median(DP)
DP_sd <- sd(DP)

# plot histogram of DP with range 0 - median(DP) + 2*sd(DP)
#hist(DP,500)
DP_plot1 <- ggplot(data.frame(DP = DP), aes(x = DP)) +
  geom_histogram(bins = 500) +
  geom_vline(xintercept = c(5000,cov*length(SAMPLES)), color = c("blue","red")) +
  annotate("text", y = Inf, 
           x = c(5000,cov*length(SAMPLES)), 
           label = c("Default","Recomended"), size=3, angle=-90, vjust=-0.4, hjust=0) +
  xlim(0,600000) +
  theme_classic() +
  labs(title = "Read depth distribution", subtitle = paste0("Recomended threshold: ", cov*length(SAMPLES)),
       y = "Count")

#hist(DP[ DP < (DP_med + 2*DP_sd) ], 500, main = "DP 0 - median + 2 sd" )
DP_plot2 <- ggplot(data.frame(DP = DP[ DP < (DP_med + 2*DP_sd) ]), aes(x = DP)) +
  geom_histogram(bins = 500) +
  geom_vline(xintercept = c(5000,cov*length(SAMPLES)), color = c("blue","red")) +
  annotate("text", y = Inf, 
           x = c(5000,cov*length(SAMPLES)), 
           label = c("Default","Recomended"), size=3, angle=-90, vjust=-0.4, hjust=0) +
  theme_classic() +
  labs(title = "Read depth distribution", subtitle = paste0("Recomended threshold: ", cov*length(SAMPLES)),
       y = "Count", x = "DP (0 - median + 2 sd)")

#hist(MQ,60)
MQ_plot <- ggplot(data.frame(MQ = MQ), aes(x = MQ)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(58.75,61.25), color = c("blue","blue")) +
  xlim(c(55,65)) +
  theme_classic() +
  labs(title = "Mapping quality distribution", y = "Count")

#hist(VQSLOD)
#hist(VQSLOD[VQSLOD > 0 & VQSLOD < 40],45)
VQSLOD_plot <- ggplot(data.frame(VQSLOD = na.omit(VQSLOD)), aes(x = VQSLOD)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(7.81), color = c("blue")) +
  xlim(c(0,30)) +
  theme_classic() +
  labs(title = "VQSLOD distribution", y = "Count")

pdf( paste0( input, "_vcf_metrics.pdf"), width = 9, height = 3 )
DP_plot2 + MQ_plot + VQSLOD_plot
dev.off()

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
