# load the library
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rsamtools")
library(Rsamtools)
library(tidyverse)
library(seqinr)
library(reshape2)
library(ggplot2)
library(viridis)
library(future.apply)
plan(multiprocess)
options(stringsAsFactors = F)

# setwd("G:/My Drive/work/2019/2019-12-04_Maireid_epitope_diveristy/script")
readBAM <- function(bamFile) {
    bam <- scanBam(bamFile, param=ScanBamParam(what=c("rname", "pos", "seq", "cigar")))

    # A function for collapsing the list of lists into a single list
    # as per the Rsamtools vignette
    .unlist <- function(x) {
        x1 <- x[[1L]]
        if (is.factor(x1)) {
            structure(unlist(x), class = "factor", levels = levels(x1))
        } else {
            do.call(c, x)
        }
    }

    bam_field <- names(bam[[1]])

    list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

    bam_df <- do.call("DataFrame", list)
    names(bam_df) <- bam_field

    # return a list that can be called as a data frame
    return(bam_df)
}

# correct the peptide file
ref <- read.fasta("../data/reference.fasta")
peptides0 <- read_csv("../data/diversity_peptides_index.csv",  na = "")
peptides0$Name <- as.character(peptides0$Name)
peptides1 <- read_csv("../data/diversity_peptides_index_120.csv",  na = "")
peptides <- bind_rows(peptides0, peptides1)
# 120
nt_start <- sapply(peptides[[2]], function(x){
    as.numeric(strsplit(x, "-")[[1]][1])
})
nt_end <- sapply(peptides[[2]], function(x){
    as.numeric(strsplit(x, "-")[[1]][2])
})
peptides[[2]] <- 
    paste0(round((nt_start + nt_end)/2) - 60, "-", round((nt_start + nt_end)/2) + 60)

rst <- peptides

sample.list <- list.files("../data/", "bam$", recursive = F, full.names = F)
sample.list_f <- list.files("../data/", "bam$", recursive = F, full.names = T)
for (dir.i in 1:length(sample.list)) {
    # Load the bam file
    print(dir.i)
    sample.name <- strsplit(sample.list[dir.i], "_")[[1]][1]
    bam <- readBAM(sample.list_f[dir.i])

    hs_nt <- rep(NA, nrow(peptides))
    hgs_nt <- hs_nt
    hs_aa <- hs_nt
    hgs_aa <- hs_nt

    for (i in 1:nrow(peptides)) {
        # print(i)
        # sapply(peptides$Nucleotide.pos..on.Consensus.seq.,function(x){
        #     tmp = strsplit(x, "-")[[1]]
        #     tmp = as.numeric(tmp)
        #     (tmp[2]-tmp[1]+1)/3
        #     })

        # determine segment
        cur_segment <- peptides$Segment[i]
        ref.name <- names(ref)[grep(cur_segment, names(ref))]
        ref.seq <- ref[[grep(cur_segment, names(ref))]]

        seq <- bam$seq[which(bam$rname == ref.name)]
        cigar <- bam$cigar[which(bam$rname == ref.name)]
        # flag = bam$flag[which(bam$rname==ref.name)]
        pos <- bam$pos[which(bam$rname == ref.name)]
        # width = bam$qwidth[which(bam$rname==ref.name)]

        which.t <- which(is.na(cigar))
        if(length(which.t)>0){
            seq <- seq[-which.t]
            cigar <- cigar[-which.t]
            pos <- pos[-which.t]
        }
        
        # get the seqs by cigar
        seq.clip <- future_mapply(function(x, y) {
            cigar.1 <- unlist(strsplit(x, "\\d+"))
            cigar.2 <- unlist(strsplit(x, "\\D"))
            seq.clip.tmp <- c()
            start <- 1
            print(x)
            for (j in 1:length(cigar.2)) {
                num <- as.numeric(cigar.2[j])
                move <- cigar.1[j + 1]
                print(move)
                if (move %in% c("M", "=", "X")) {
                    seq.clip.tmp <- paste0(seq.clip.tmp, c2s(s2c(as.character(y))[start:(start + num - 1)]))
                    start <- start + num
                } else if (move == "D") {
                    seq.clip.tmp <- paste0(seq.clip.tmp, c2s(rep("-", num)))
                } else {
                    start <- start + num
                }
            }
            return(seq.clip.tmp)
        }, cigar, seq)

        length.seq <- sapply(seq.clip, nchar)

        peptide.name <- as.character(peptides[i, 1])
        ## These positions are in NT sense
        peptide.pos1 <- as.numeric(strsplit(peptides[[i, 2]], "-")[[1]][1])
        peptide.pos2 <- as.numeric(strsplit(peptides[[i, 2]], "-")[[1]][2])

        check1 <- pos <= peptide.pos1
        check2 <- (pos + length.seq - 1) >= peptide.pos2
        check.t <- check1 & check2
        # sum(check1)
        # sum(check2)
        # sum(check.t)
        if(sum(check.t)>0){
            pep.tmp <- seq.clip[check.t]
            pos.tmp <- pos[check.t]

            ## NT
            pep.nt.tmp <- mapply(
                function(x, y) {
                    substr(x, peptide.pos1 - y + 1, peptide.pos2 - y + 1)
                },
                pep.tmp, pos.tmp
            )
            df.tmp <- sapply(pep.nt.tmp, function(x) {
                s2c(x)
            })
            df.tmp <- as.data.frame(t(df.tmp))

            df.tmp.agg <- aggregate(cbind(df.tmp[0], Count = 1), by = df.tmp, length)
            df.tmp.agg <- df.tmp.agg[order(df.tmp.agg[, ncol(df.tmp.agg)], decreasing = T),]
            names(df.tmp.agg)[1:(ncol(df.tmp.agg) - 1)] <- peptide.pos1:peptide.pos2
            write_csv(df.tmp.agg, paste0(
                "../results/",
                sample.name, "_", cur_segment, "_",
                peptide.name, "_nt.csv"
            ))
            ### Simpson's and Shannon's diversity index
            freq_tmp <- df.tmp.agg$Count / sum(df.tmp.agg$Count)
            # freq_tmp <- c(0.463,0.062,0.039,0.027,0.037,0.016,0.033,0.054,0.248,0.02)

            #### HS: Shannon's entropy
            hs_nt[i] <- -sum(freq_tmp * log(freq_tmp))
            #### Gini-Simpson index
            hgs_nt[i] <- 1 - sum(freq_tmp * freq_tmp)

            ## AA
            pep.aa.tmp <- sapply(pep.nt.tmp, function(x) {
                translate(s2c(x))
            })
            df.tmp <- as.data.frame(t(pep.aa.tmp))

            df.tmp.agg <- aggregate(cbind(df.tmp[0], Count = 1), by = df.tmp, length)
            df.tmp.agg <- df.tmp.agg[order(df.tmp.agg[, ncol(df.tmp.agg)], decreasing = T),]
            aa_start_pos <- round(as.numeric(strsplit(peptides[[i, 2]], "-")[[1]][1])/3)
            names(df.tmp.agg)[1:(ncol(df.tmp.agg) - 1)] <- aa_start_pos:(aa_start_pos+ncol(df.tmp.agg) - 2)
            write_csv(df.tmp.agg, paste0(
                "../results/",
                sample.name, "_", cur_segment, "_",
                peptide.name, "_aa.csv"
            ))
            ### Simpson's and Shannon's diversity index
            freq_tmp <- df.tmp.agg$Count / sum(df.tmp.agg$Count)
            #### HS: Shannon's entropy
            hs_aa[i] <- -sum(freq_tmp * log(freq_tmp))
            #### Gini-Simpson index
            hgs_aa[i] <- 1 - sum(freq_tmp * freq_tmp)
        }         
    }

    ##ToRun
    rst <- cbind(rst, hs_nt, hgs_nt, hs_aa, hgs_aa)
    names(rst)[(ncol(rst) - 3):ncol(rst)] <- paste0(
        sample.name, "_",
        names(rst)[(ncol(rst) - 3):ncol(rst)]
    )
}

write.csv(rst, "../results/diveristy_result.csv", row.names = F)

