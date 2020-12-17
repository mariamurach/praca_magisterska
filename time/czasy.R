library(tidyverse)
time <- read_csv("czasy.csv")
ggplot(time, aes(x = program, y = time$`czas [s]`)) + geom_bar(stat = "identity")
