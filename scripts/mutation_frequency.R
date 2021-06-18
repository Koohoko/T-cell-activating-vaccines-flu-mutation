#!/usr/bin/Rscript
library(tidyverse)
library('seqinr')

files.name <- list.files("../results/", "lofreq")
sample_names <- sapply(files.name, function(x) {
    tmp <- strsplit(x, "_", fixed = TRUE)[[1]][2]
    tmp <- strsplit(tmp, ".vcf", fixed = TRUE)[[1]][1]
    return(tmp)
})
files.path <- list.files("../results/", "vcf", full.names = T)

data <- lapply(files.path, function(x){
    tmp <- read_delim(x, delim = "\t", comment = "##")[,1:8]
    sample_name <- strsplit(x, "/", fixed = TRUE)[[1]]
    sample_name <- sample_name[length(sample_name)]
    sample_name <- strsplit(sample_name, "_", fixed = TRUE)[[1]][2]
    sample_name <- strsplit(sample_name, ".vcf", fixed = TRUE)[[1]][1]
    tmp$sample <- sample_name
    
    var_caller <- strsplit(x, "/", fixed = TRUE)[[1]]
    var_caller <- var_caller[length(var_caller)]
    var_caller <- strsplit(var_caller, "_", fixed = TRUE)[[1]][1]
    tmp$from <- var_caller
    
    return(tmp)
})

data <- bind_rows(data)

data$read_depth <- vapply(data$INFO, function(x){
    tmp <- strsplit(x, "DP=", fixed = T)[[1]][2]
    strsplit(tmp, ";", fixed = TRUE)[[1]][1]
}, character(1), USE.NAMES = F)

data$allele_frequency <- vapply(data$INFO, function(x){
    tmp <- strsplit(x, ";AF=", fixed = T)[[1]][2]
    strsplit(tmp, ";", fixed = TRUE)[[1]][1]
}, character(1), USE.NAMES = F)

names(data)[1] <- "Segment"
data <- data %>% select(Segment, POS, REF, ALT, allele_frequency, 
                        read_depth, sample, from)

## filter low frequency
data <- data %>% filter(allele_frequency > 0.005)

## remove INDELs
data_snps <- data %>% filter(nchar(REF) == 1) %>% filter(nchar(ALT) == 1)

data_indels <- setdiff(data, data_snps)
data_indels <- data_indels %>% arrange(sample, Segment, POS, from)
write_csv(data_indels, "../results/data_indels.csv")
rm(data_indels)

## add codon
orf <- read.csv("../data/ORF_WT.csv", na.strings = '')
fasta <- read.fasta("../data/reference.fasta")

tmp <- apply(data_snps, 1, function(x){
    tmp <- orf[,1] == x[["Segment"]]
    pos <- as.numeric(x[["POS"]])
    ## check ORF
    start <- orf[tmp,2]
    stop <- orf[tmp,3]
    check <- all(pos >= start, pos <= stop)
    
    if(check){
        orf_position <- pos - start + 1
    } else {
        orf_position <- NA
    }
    
    ## extract codon
    if(!is.na(orf_position)){
        sequence_tmp <- fasta[[x[["Segment"]]]][start:stop]
        if (orf_position %% 3 == 0){
            ori_codon <- paste(sequence_tmp[(orf_position-2):orf_position], 
                           collapse = "")
            alt_codon <- paste0(substr(ori_codon, 1, 2), x[["ALT"]])
        } else if (orf_position %% 3 == 1){
            ori_codon <- paste(sequence_tmp[orf_position:(orf_position+2)], 
                           collapse = "")
            alt_codon <- paste0(x[["ALT"]], substr(ori_codon, 2, 3))
        } else {
            ori_codon <- paste(sequence_tmp[(orf_position-1):(orf_position+1)], 
                           collapse = "")
            alt_codon <- paste0(substr(ori_codon, 1, 1), x[["ALT"]], 
                                substr(ori_codon, 3, 3))
        }
        ori_AA <- seqinr::translate(s2c(ori_codon))
        alt_AA <- seqinr::translate(s2c(alt_codon))
        return(c(orf_position, ori_codon, alt_codon, ori_AA, alt_AA))
    } else {
        return(c(NA, NA, NA, NA, NA))
    }
})

tmp <- as_tibble(t(tmp))
names(tmp) <- c("orf_position", "ori_codon", "alt_codon", "ori_AA", "alt_AA")
data_snps <- bind_cols(data_snps, tmp)
data_snps$Silent_mutation <- data_snps$ori_AA == data_snps$alt_AA

### add Allele frequency
files.path <- list.files("../FastQC/", ".mpileup.txt",full.names = T)
df_readcount <- lapply(files.path, function(x){ 
    tmp <- read_delim(x, delim = "\t", comment = "##")
    sample_name <- strsplit(x, "/", fixed = TRUE)[[1]]
    sample_name <- sample_name[length(sample_name)]
    sample_name <- strsplit(sample_name, "_", fixed = TRUE)[[1]][1]
    tmp$sample <- sample_name
    
    return(tmp)
})

df_readcount <- bind_rows(df_readcount)
names(df_readcount)[1:2] <- c("Segment", "POS")
df_readcount <- df_readcount %>% select(-identifier, -indelcount)
df_readcount <- df_readcount[!duplicated(df_readcount),]

write_csv(df_readcount, "../results/readcount.csv")
write_csv(data_snps, "../results/data_snps.csv")

## by tools
data_snps_by_tools <- data_snps %>% select(-read_depth) %>% 
    mutate(allele_frequency = 1) %>% 
    spread(from, allele_frequency, fill = 0)
data_snps_by_tools$n_found_in <- data_snps_by_tools$freebayes +
    data_snps_by_tools$lofreq + data_snps_by_tools$vardict
data_snps_by_tools_rc <- left_join(data_snps_by_tools,
                                   df_readcount, c("sample", "Segment", "POS"))
AF_alignment <- apply(data_snps_by_tools_rc, 1, function(x){
    alt <- tolower(x[["ALT"]])
    tmp <- paste0(alt, "count")
    return(as.numeric(x[[tmp]]) / as.numeric(x[["depth"]]))
})
data_snps_by_tools$AlleleFrequency_alignment <- AF_alignment

## by sample
data_snps_by_samples_1 <- data_snps_by_tools %>% 
    mutate(sample = paste0(sample, "_foundin")) %>% 
    select(-freebayes, -vardict, -lofreq, -AlleleFrequency_alignment) %>% 
    spread(sample, n_found_in, fill = 0)

sample_names <- unique(df_readcount$sample)
data_snps_by_samples_2 <- data_snps_by_samples_1 %>% select(Segment:Silent_mutation)
tmp <- matrix(nrow = nrow(data_snps_by_samples_2), ncol = length(sample_names))
tmp <- as_tibble(tmp)
names(tmp) <- sample_names
data_snps_readcount <- lapply(sample_names, function(x){
    return(data_snps_by_samples_2)
})
data_snps_readcount <- bind_rows(data_snps_readcount)
data_snps_readcount <- bind_cols(data_snps_readcount,
                                    gather(tmp, key = sample)) %>% select(-value)
data_snps_readcount <- left_join(data_snps_readcount, 
                                 df_readcount,  c("sample", "Segment", "POS"))
AF_alignment <- apply(data_snps_readcount, 1, function(x){
    alt <- tolower(x[["ALT"]])
    tmp <- paste0(alt, "count")
    return(as.numeric(x[[tmp]]) / as.numeric(x[["depth"]]))
})
data_snps_readcount$AlleleFrequency_alignment <- AF_alignment
check <- data_snps_readcount %>% select(Segment:ALT, sample) %>% duplicated()
data_snps_readcount <- data_snps_readcount[!check,]
data_snps_readcount <- data_snps_readcount %>% 
    mutate(sample = paste0(sample, "_AF")) %>% 
    select(Segment:sample, AlleleFrequency_alignment) %>%
    spread(sample, AlleleFrequency_alignment, fill = "")

data_snps_by_samples_2 <- data_snps_by_tools %>% 
    mutate(sample = paste0(sample, "_AF_masked")) %>% 
    select(-freebayes, -vardict, -lofreq, -n_found_in) %>% 
    spread(sample, AlleleFrequency_alignment, fill = "")
data_snps_by_samples <- left_join(data_snps_by_samples_1, data_snps_by_samples_2)
data_snps_by_samples <- left_join(data_snps_by_samples, data_snps_readcount)

data_snps_by_samples <- data_snps_by_samples %>% arrange(Segment, POS)

write_csv(data_snps_by_samples, "../results/data_snps_by_samples.csv")

