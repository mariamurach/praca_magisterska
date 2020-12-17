library(tidyverse)
sniffles_minimap <- read_tsv("short_stats/sniffles_minimap.tsv", comment = "#") %>% mutate(program = "Sniffles+minimap")
lumpy <- read_tsv("short_stats/lumpy.tsv", comment = "#") %>% mutate(program = "Lumpy")
nanosv <-  read_tsv("short_stats/nanosv.tsv", comment = "#") %>% mutate(program = "nanoSV")
sniffles <-  read_tsv("short_stats/sniffles.tsv", comment = "#") %>% mutate(program = "Sniffles+NGMLR")

results <- bind_rows(sniffles_minimap, lumpy, nanosv, sniffles)
results$`[3]key` <- sub(':', '', results$`[3]key`)
results$`[3]key` <- sub('number of ', '', results$`[3]key`)
results <- results %>% rename(typ_wariantu  = `[3]key`, n = `[4]value`)
results <- results %>% select(typ_wariantu , n, program)
results$typ_wariantu  <- as_factor(results$typ_wariantu )
warianty <- results %>% filter(!(n == 0)) %>% filter(!(typ_wariantu  == "samples")) %>% filter(!(typ_wariantu  == "records"))

rekordy <- results %>% filter(typ_wariantu == "records")
ggplot(rekordy, aes(program, n, fill = program)) +
  geom_bar(stat = "identity") + 
  labs(y = "ilość rekordów w pliku VCF")

ggplot(warianty, aes(typ_wariantu , n, fill = typ_wariantu )) +
  geom_bar(stat = 'identity') +
  labs(x = "typ wariantu", y = "lczba znalezionych wariantów") +
  facet_wrap(~program)

