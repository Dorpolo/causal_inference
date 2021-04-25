library(dplyr)
library(ggplot2)
library(reshape2)
library(stats)
library(reticulate)

data <- read.csv('~/Desktop/private/TAU/tau_stats/data/rhc.csv', header = T, sep = ';')

install.packages('MatchIt')