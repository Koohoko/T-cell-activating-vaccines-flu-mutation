#!/usr/bin/Rscript

if (!requireNamespace("ggplot2")){
    install.packages("ggplot2")
}
library('ggplot2')
library(ggrepel)

# bam_read_count ----------------------------------------------------------

# ridges_plot -------------------------------------------------------------
if (!requireNamespace("ggridges")){
    install.packages("ggridges")
}
if (!requireNamespace("viridis")){
    install.packages("viridis")
}
library('ggridges')
library('viridis')
library("tidyverse")

files = list.files("../FastQC", "mpileup", full.names = T)

#combine readcounts of different samples
df = data.frame()
for(i in 1:length(files)){
    sample.name = strsplit(files[i], "_")[[1]][1]  
    sample.name = strsplit(sample.name, "/")[[1]][3] 
    tmp = read.table(files[i], stringsAsFactors = F, header = T)
    tmp = cbind(rep(sample.name,nrow(tmp)),tmp)
    df = rbind(df,tmp)
}

#data clean and calculate alle_freq
names(df)[1:2] = c('Sample', "Segment")
df <- as_tibble(df) %>% filter(grepl("^d", Sample))
df$Sample <- factor(df$Sample)

Alle_freq = df$altcount / df$depth
Alle_freq[is.na(Alle_freq)] <- 0

Nt = unlist(lapply(table(df$Sample), function(x) {
    1:x
}))
sample_seq <- seq_along(table(df$Sample))
Y = rep(sample_seq, table(df$Sample))
df = cbind(df, Alle_freq, Nt, Y)
df$log10.depth <- log10(df$depth)
head(df)
y.height <- quantile(df$depth, seq(0, 1, 0.01))[96]
y.log.height <- quantile(df$log10.depth, seq(0, 1, 0.01))[96]
#combined Plot
df$Sample <- factor(df$Sample)
table(df$Sample)

p <- ggplot(data=df) + 
    geom_ridgeline_gradient(aes(x=Nt,y=Y*y.height,
        height=depth,group=Sample,fill=Segment),alpha=0.8)+
    scale_fill_viridis(discrete = TRUE, direction = -1,alpha=0.8)+
    scale_y_continuous(breaks = (sample_seq)*y.height, 
        label = levels(df$Sample)) + 
    scale_x_continuous(breaks=seq(0,13562,1000)) +
    ylab('Sample') #+theme_ridges()
ggsave('../FastQC/Coverage_d357.pdf', plot = p,height = 17,width = 9)

p <- ggplot(data=df) + 
    geom_ridgeline_gradient(aes(x=Nt,y=Y*y.log.height,
        height=log10.depth,group=Sample,fill=Segment),alpha=0.8)+
    scale_fill_viridis(discrete = TRUE, direction = -1,alpha=0.8)+
    scale_y_continuous(breaks = (sample_seq)*y.log.height, 
        label = levels(df$Sample)) + 
    scale_x_continuous(breaks=seq(0,13562,1000)) +
    ylab('Sample') #+theme_ridges()
ggsave('../FastQC/Coverage_log10_d357.pdf', plot = p, height = 17,width = 9)

ggplot(data=df) +
    geom_ridgeline_gradient(aes(x=Nt,y=Y,
        height=Alle_freq,group=Sample,fill=Segment))+
    scale_fill_viridis(discrete = TRUE, direction = -1)+
    scale_y_continuous(breaks = sample_seq,label = levels(df$Sample)) +
    ylab('Sample')
ggsave('../FastQC/Alle_freq.jpeg',device = 'jpeg',height = 17,width = 9)

times = round(df$Alle_freq*1000)+1
df.density = df[rep(1:nrow(df),times),]
rownames(df.density)=c()

ggplot(data=df.density,aes(x=Nt,y=Sample,fill=Segment)) + 
    geom_density_ridges(alpha=0.8)+
    scale_fill_viridis(discrete = TRUE, direction = -1)+
    scale_x_continuous(breaks=seq(0,13562,1000))
ggsave('../FastQC/Alle_freq_smooth.jpeg',
    device = 'jpeg',height = 17,width = 9,dpi = 300)


library("tidyverse") 
data <- read_csv("../results/readcount.csv") 
tmp <- data %>% group_by(sample, Segment) %>% summarise(d=mean(depth)) %>% arrange(d) 
tmp_sample <- data %>% group_by(sample) %>% summarise(d=mean(depth)) %>% arrange(d) 
write_csv(tmp, "../results/mean_depth_sample.csv")
write_csv(tmp_sample, "../results/mean_depth_sample_segment.csv")

# Q33 checking
fastqc_zip_files <- list.files("../FastQC/", "fastqc.zip", full.names = T)
fastqc_zip_files_name <- list.files("../FastQC/", "fastqc.zip")
fastqc_zip_sample_name <- sapply(fastqc_zip_files_name, function(x){
    strsplit(x, "_fastqc")[[1]][1]
})
df <- lapply(seq_along(fastqc_zip_files), function(i){
    zip_file_tmp <- fastqc_zip_files[i]
    zip_file_name_tmp <- fastqc_zip_files_name[i]
    zip_sample_name_tmp <- fastqc_zip_sample_name[i]
    system(paste0("unzip ", zip_file_tmp))
    target_line_num <- grep(">>Per sequence quality scores", readLines(paste0("./", zip_sample_name_tmp, "_fastqc/fastqc_data.txt")))
    target_line_num_end <- grep(">>Per base sequence content", readLines(paste0("./", zip_sample_name_tmp, "_fastqc/fastqc_data.txt")))
    df_tmp <- read_delim(paste0("./", zip_sample_name_tmp, "_fastqc/fastqc_data.txt"), "\t", skip = target_line_num, n_max = target_line_num_end-target_line_num-3)
    names(df_tmp)[1] <- "Quality"
    df_tmp$propotion <- df_tmp$Count/sum(df_tmp$Count)
    df_tmp$sample <- zip_sample_name_tmp
    system(paste0("rm -rf ./", zip_sample_name_tmp, "_fastqc/"))
    return(df_tmp)
})
df <- bind_rows(df)
df_q33 <- df %>% group_by(sample) %>% summarise(proportion_over_or_equal_Q33 = sum(propotion[Quality>=33]), Quality = 33, propotion = propotion[which(Quality==33)]) %>% arrange(proportion_over_or_equal_Q33)
df_q33_for_plot <- df_q33 %>% filter(proportion_over_or_equal_Q33 < 0.9)

df %>% ggplot(aes(x = Quality, y = propotion, color = sample))+
    geom_line()+
    theme(legend.position = "none") +
    geom_text_repel(aes(label = sample), data = df_q33_for_plot)

write_csv(df_q33[,1:2], "../results/df_q33.csv")
write_csv(df, "../results/df_q_all.csv")
