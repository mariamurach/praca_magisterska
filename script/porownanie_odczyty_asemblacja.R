library(tidyverse)
ath <- read_csv("ath_cnv.csv", skip = 4) %>% filter(`Included in AthCNV dataset?` == "YES") %>% select(CNV_ID, Chr, Start, Stop, size, `Included in AthCNV dataset?`)
ath$Chr <- factor(gsub("Chr", "", ath$Chr))
asemblacja <- read_csv("asemblacja/pokrycie_AthCNV_asemblacja.csv") 
odczyty <- read_csv("pokrycie_AthCNV_odczyty.csv")

all <- bind_rows(asemblacja, odczyty)


all %>% ggplot(aes())
