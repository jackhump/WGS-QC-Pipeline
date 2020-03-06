# for a WGS QC run
# find the number of variants in each stats file
# get the total number of remaining variants at each step of the pipeline

library(dplyr)
library(stringr)
library(readr)
#statsFolder <- "test/stats/"

args <- commandArgs(trailingOnly = TRUE)

statsFolder <- args[1]
outFile <- args[2]

cmd <- paste0("grep \"number of records:\" ", statsFolder, "*stats.txt" )
all <- system(cmd, intern = TRUE)

# parse output
stat <- str_split_fixed(all, "\t", 4)[,1]
stat <- gsub(".stats.txt:SN", "", stat)
stat <- basename(stat)

# split file name into chr, chunk and step
chr <- str_split_fixed(stat, "_", 3)[,1]
chunk <- str_split_fixed(stat, "_", 3)[,2]
step <- str_split_fixed(stat, "_", 3)[,3]

# get out number of variants
n_var <- as.numeric(str_split_fixed(all, "\t", 4)[,4])

df <- data.frame(file = stat, chr = chr, chunk = chunk, step = step,  n_var = n_var, stringsAsFactors = FALSE)

# combine Filter6 for SNPs and Indels
indels <-  filter(df, step == "Filter6_Indels")

# remove these counts as they are same as Filter6 when summed together
# Filter_7 SNPs is the same as Filter7 + the Filter6_Indels
df <- filter(df,! step %in% c("Filter6_Indels", "Filter6_SNPs", "Filter7_SNPs", "Chunk") )

# for filter 7 - sum Filter 7 and indel counts together
df$step <- ifelse( df$step == "CombineSNPsIndels", yes = "Filter7", no = df$step)

# fix chrAll steps and chunking steps
df$step <- ifelse(grepl("_filtered_stats", df$file), yes = "FilterBlacklist", no = df$step)


df$step <- ifelse(df$step == "", yes = df$chunk, no = df$step)

# remove chromosome concat count - duplicate
df <- filter(df, file != "chrAll_Recombined")

# sum chunked steps together
summary <- df %>% group_by(step) %>%
    summarise( n_var = sum(n_var) ) %>% 
    arrange( desc(n_var) )

# put together dictionary of steps:
filter_dict <- c(
"Initial" = "Initial joint VCF",
"BlacklistFiltered" = "Remove samples and variants with blacklist region",
"separateBiallelic" = "Only biallellic variants kept",
"Filter1" = "GATK PASS",
"Filter2" = "Genotype Depth (GDP)",
"Filter3" = "Genotype Quality (GQ)",
"Filter4" = "Variant Missingness",
"Filter5" = "Total Read Depth (DP)",
"Filter6" = "Mapping Quality (MQ)",
"Filter7" = "Variant Quality (VQSLOD)",
"Filter8" = "Inbreeding Coefficient",
"Filter9" = "Sample Missingness",
"full" = "Sample Relatedness",
"MAF0.01" = "Minor alelle frequency >= 0.01"
)

summary$filter <- filter_dict[ match(summary$step, names(filter_dict) ) ]

readr::write_tsv(summary, outFile)
save.image(file = paste0(statsFolder, "all_stats.RData"))
