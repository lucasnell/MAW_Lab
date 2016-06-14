
library("dplyr")
library("tidyr")
library("stringr")
library("readr")
library("ggplot2")


tagDF <- read_tsv('/Volumes/MW_18TB/Lucas_Nell/musGBS/allTags.tsv.gz', 
                  progress = TRUE) %>% 
    arrange(Chromosome, desc(Strand), Start)

tagDF %>% group_by(Chromosome, Strand)

tagDF %>%
    filter(Chromosome %in% paste0('chr', c(seq(19), 'X'))) %>%
    group_by(Sample_ID) %>%
    summarize(depMean = mean(Reads), 
              depMedian = median(Reads),
              depMin = min(Reads),
              depMax = max(Reads))


tagDF %>% filter(Reads == max(tagDF$Reads))



# If you want to plot densities
tagDF %>%
    filter(Chromosome != 'chrM') %>%
    ggplot(aes(x = Reads, color = factor(Sample_ID))) +
    geom_density(fill = NA) + 
    scale_x_log10('Number of reads per "stack"', breaks = c(4, 10, 100, 1000)) +
    ylab('Density') +
    scale_color_manual('sample', 
                       values = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
                                  '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'))



# If you want to plot actual count numbers
tagDF %>%
    filter(Chromosome != 'chrM') %>%
    ggplot(aes(x = Reads, ..count.., color = factor(Sample_ID))) +
    geom_density(fill = NA) + 
    scale_x_log10('Number of reads per "stack"', breaks = c(4, 10, 100, 1000)) +
    scale_y_continuous('Millions of "stacks"', 
                       labels = function(x){x/1e6}) +
    scale_color_manual('sample', 
                       values = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
                                  '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'))



# Percent and total # of stacks with # reads >= 4
tagDF %>% 
    filter(Chromosome != 'chrM') %>%
    group_by(Sample_ID) %>%
    summarize(percent = (length(Reads[Reads >= 4]) * 100) / n(),
              total = format(length(Reads[Reads >= 4]), big.mark = ','))

# Sample_ID    percent       total
#     (int)      (dbl)       (dbl)
#       697   3.190598      41,808
#       858   7.828415     106,861
#      1137   4.780112      41,581
#      1501   5.801726      59,801
#      1559   6.161988      59,743
#      1810   4.346358      34,104
#      1866   2.432416      13,448
#      1907   1.911653      10,407
#      1941   4.675852      45,206


