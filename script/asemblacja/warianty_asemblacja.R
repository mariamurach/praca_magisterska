library(tidyverse)
ath <- read_csv("../ath_cnv.csv", skip = 4) %>% filter(`Included in AthCNV dataset?` == "YES") %>% select(CNV_ID, Chr, Start, Stop, size, `Included in AthCNV dataset?`)
ath$Chr <- factor(gsub("Chr", "", ath$Chr))

svmu <- read_tsv("../../cnv/svmu_sv.txt")
assemblytics <- read_tsv("../../cnv/Arabidopsis.Assemblytics_structural_variants.bed")
showdiff <- read_tsv("../../cnv/showdiff.diff", skip = 3)

svmu <- svmu %>% mutate(CHROM = REF_CHROM, START = REF_START, STOP = REF_END, type = SV_TYPE, SVLEN = abs(LEN)) %>% select(CHROM, START, STOP, SVLEN, type)
svmu$program <- "svmu"
assemblytics <- assemblytics %>% mutate(CHROM = `#reference`, START = ref_start, STOP = ref_stop, type = type, SVLEN = abs(size)) %>% select(CHROM, START, STOP, SVLEN, type)
assemblytics$program <- "assemblytics" 
showdiff <- showdiff %>% mutate(CHROM = `[SEQ]`, START = `[S1]`, STOP = `[E1]`, type = `[TYPE]`, SVLEN = abs(`[LEN 1]`)) %>% select(CHROM, START, STOP, SVLEN, type)
showdiff$program <- "showdiff"

assemblytics$wariant <- ifelse(assemblytics$type == "Insertion", "INS", 
                                           ifelse(assemblytics$type == "Deletion", "DEL", "INNY"))
assemblytics$type <- factor(assemblytics$type)
assemblytics_insertions <- assemblytics %>% filter(type == "Insertion") %>% nrow
assemblytics_deletions <- assemblytics %>% filter(type == "Deletion") %>% nrow
assemblytics_others <- assemblytics %>% filter(!(type == "Deletion")) %>% filter(!(type == "Insertion")) %>% nrow


showdiff$type <- factor(showdiff$type)
showdiff_insertions <- showdiff %>% filter(type == "BRK") %>% nrow
showdiff_deletions <- showdiff %>% filter(type == "GAP") %>% nrow
showdiff_others <- showdiff %>% filter(!(type == "BRK")) %>% filter(!(type == "GAP")) %>% nrow

showdiff$wariant <- ifelse(showdiff$type == "BRK", "INS", 
                       ifelse(showdiff$type == "GAP", "DEL", "INNY"))


svmu$type <- factor(svmu$type)

svmu$type <- factor(svmu$type)
svmu_insertions <- svmu %>% filter(type == "INS") %>% nrow
svmu_deletions <- svmu %>% filter(type == "DEL") %>% nrow
svmu_others <- svmu %>% filter(!(type == "DEL")) %>% filter(!(type == "INS")) %>% nrow

svmu$wariant <- ifelse(svmu$type == "INS", "INS", 
                       ifelse(svmu$type == "DEL", "DEL", "INNY"))

stats_wariants <- tibble(program = c("assemblytics", "showdiff", "svmu"), 
                         insercje = c(assemblytics_insertions, showdiff_insertions, svmu_insertions),
                         delecje = c(assemblytics_deletions, showdiff_deletions, svmu_deletions),
                         inne_cnv = c(assemblytics_others, showdiff_others, svmu_others))

stats_wariants <- stats_wariants %>% pivot_longer(., cols = c("insercje", "delecje", "inne_cnv"), names_to = "typ_wariantu", values_to = "ilość")
stats_wariants$typ_wariantu <- factor(stats_wariants$typ_wariantu, levels =  c("insercje", "delecje", "inne_cnv") )
stats_wariants$program <- factor(stats_wariants$program)

warianty <- stats_wariants %>%group_by(program) %>% summarize(sum = sum(ilość))
warianty %>% ggplot(aes(x = program, y = sum, fill = program)) +
  geom_bar(stat = "identity") + 
  labs( y = "n") + 
  geom_text(aes(label=sum), position=position_dodge(width=0.9), vjust=-0.25)

stats_wariants %>% ggplot(aes(typ_wariantu, ilość, fill= typ_wariantu)) + 
  geom_bar(stat = "identity", color = "white", width = 1) +
  geom_text(aes(label=ilość), position=position_dodge(width=0.9), vjust=-0.25)+ 
  facet_wrap(stats_wariants$program)

warianty <- bind_rows(svmu, assemblytics, showdiff)
warianty$length <- ""
warianty$length <- ifelse(abs(warianty$SVLEN) >= 500, ">=500bp", 
                          ifelse( abs(warianty$SVLEN) < 500 & abs(warianty$SVLEN) >= 50, "50-499bp", "< 50bp" ))
warianty$length <- factor(warianty$length, levels = c("< 50bp", "50-499bp", ">=500bp" ))
write_csv(warianty, "wszystkie_warianty_asemblacja.csv")
warianty %>% ggplot(aes(x = program, fill = length)) + 
  geom_bar(position = "fill") + 
  labs(x = "program", y = "Ułamek wariantów") + 
  scale_fill_discrete(name = "Długość")

warianty$wariant <- factor(warianty$wariant, levels  = c("INS", "DEL", "INNY"))

warianty %>% ggplot(aes(x = program, fill = wariant)) + 
  geom_bar(position = "fill") + 
  labs(x = "program", y = "Ułamek wariantów") + 
  scale_fill_discrete(name = "Typ wariantu")


warianty_500 <- warianty %>% filter(SVLEN >= 500) %>% filter(!is.na(CHROM)) %>% arrange(CHROM, START)
warianty_500$program <- factor(warianty_500$program)
war <- warianty_500 %>% count(program, sort = TRUE)

warianty_500 %>% ggplot(aes(x = program, fill = wariant)) + 
  geom_bar(position = "fill") + 
  labs(x = "program", y = "Ułamek wariantów") + 
  scale_fill_discrete(name = "Typ wariantu")

war_500 <- warianty_500 %>%group_by(program) %>% count
war_500 %>% ggplot(aes(x = program, y = n, fill = program)) +
  geom_bar(stat = "identity") + 
  labs( y = "n") + 
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25)



coverage <- 0.2
overlap_02 <- ath
overlap_02$assemblytics <- 0
overlap_02$showdiff <- 0
overlap_02$svmu <- 0

for(index in 1:nrow(warianty_500))
{
  cnv <- warianty_500 %>% slice(index)
  chrom <- ath %>% 
    filter( Chr == cnv$CHROM, 
            cnv$START < Start & cnv$STOP< Stop & (cnv$STOP - Start >= coverage*size & cnv$STOP - Start >= coverage*cnv$SVLEN) | 
              cnv$START < Start & cnv$STOP > Stop & (Stop - Start >= coverage*size & Stop - Start >= coverage*cnv$SVLEN) |
              cnv$START > Start & cnv$STOP < Stop & (cnv$STOP - cnv$START >= coverage*size & cnv$STOP - cnv$START >= coverage*cnv$SVLEN) | 
              cnv$START > Start & cnv$STOP > Stop & (Stop - cnv$START >= coverage*size & Stop - cnv$START >= coverage*cnv$SVLEN))
  program <- cnv$program
  if(program == "svmu") {overlap_02$svmu[overlap_02$CNV_ID %in% chrom$CNV_ID] <- 1}
  else if(program == "assemblytics") {overlap_02$assemblytics[overlap_02$CNV_ID %in% chrom$CNV_ID] <- 1 }
  else if(program == "showdiff") {overlap_02$showdiff[overlap_02$CNV_ID %in% chrom$CNV_ID] <- 1 }
}  


coverage <- 0.5
overlap_05 <- ath
overlap_05$assemblytics <- 0
overlap_05$showdiff <- 0
overlap_05$svmu <- 0

for(index in 1:nrow(warianty_500))
{
  cnv <- warianty_500 %>% slice(index)
  chrom <- ath %>% 
    filter( Chr == cnv$CHROM, 
            cnv$START < Start & cnv$STOP< Stop & (cnv$STOP - Start >= coverage*size & cnv$STOP - Start >= coverage*cnv$SVLEN) | 
              cnv$START < Start & cnv$STOP > Stop & (Stop - Start >= coverage*size & Stop - Start >= coverage*cnv$SVLEN) |
              cnv$START > Start & cnv$STOP < Stop & (cnv$STOP - cnv$START >= coverage*size & cnv$STOP - cnv$START >= coverage*cnv$SVLEN) | 
              cnv$START > Start & cnv$STOP > Stop & (Stop - cnv$START >= coverage*size & Stop - cnv$START >= coverage*cnv$SVLEN))
  program <- cnv$program
  if(program == "svmu") {overlap_05$svmu[overlap_05$CNV_ID %in% chrom$CNV_ID] <- 1}
  else if(program == "assemblytics") {overlap_05$assemblytics[overlap_05$CNV_ID %in% chrom$CNV_ID] <- 1 }
  else if(program == "showdiff") {overlap_05$showdiff[overlap_05$CNV_ID %in% chrom$CNV_ID] <- 1 }
}  

coverage <- 0.8
overlap_08 <- ath
overlap_08$assemblytics <- 0
overlap_08$showdiff <- 0
overlap_08$svmu <- 0

for(index in 1:nrow(warianty_500))
{
  cnv <- warianty_500 %>% slice(index)
  chrom <- ath %>% 
    filter( Chr == cnv$CHROM, 
            cnv$START < Start & cnv$STOP< Stop & (cnv$STOP - Start >= coverage*size & cnv$STOP - Start >= coverage*cnv$SVLEN) | 
              cnv$START < Start & cnv$STOP > Stop & (Stop - Start >= coverage*size & Stop - Start >= coverage*cnv$SVLEN) |
              cnv$START > Start & cnv$STOP < Stop & (cnv$STOP - cnv$START >= coverage*size & cnv$STOP - cnv$START >= coverage*cnv$SVLEN) | 
              cnv$START > Start & cnv$STOP > Stop & (Stop - cnv$START >= coverage*size & Stop - cnv$START >= coverage*cnv$SVLEN))
  program <- cnv$program
  if(program == "svmu") {overlap_08$svmu[overlap_08$CNV_ID %in% chrom$CNV_ID] <- 1}
  else if(program == "assemblytics") {overlap_08$assemblytics[overlap_08$CNV_ID %in% chrom$CNV_ID] <- 1 }
  else if(program == "showdiff") {overlap_08$showdiff[overlap_08$CNV_ID %in% chrom$CNV_ID] <- 1 }
}  


write_csv(overlap_02, "overlap_20.csv")
write_csv(overlap_05, "overlap_50.csv")
write_csv(overlap_08, "overlap_80.csv")

overlap_02$pokrycie = "20%"
overlap_05$pokrycie = "50%"
overlap_08$pokrycie = "80%"

overlap <- bind_rows(overlap_02, overlap_05, overlap_08)

write_csv(overlap, "asemblacja_AthCNV_detekcja.csv")

overlaps <- overlap %>% 
  group_by(pokrycie) %>% 
  mutate(assemblytics = sum(assemblytics)*100/sum(nrow(ath)), showdiff = sum(showdiff)*100/sum(nrow(ath)), svmu = sum(svmu)*100/sum(nrow(ath))) %>%
  select(assemblytics, showdiff, svmu) %>% 
  unique %>% 
  pivot_longer(., c("assemblytics", "showdiff", "svmu"), names_to = "program", values_to = "%AthCNV")

write_csv(overlaps, "pokrycie_AthCNV_asemblacja.csv")

pd <- position_dodge(0.1)
ggplot(overlaps, aes(pokrycie, `%AthCNV`, color = program, group = program))  +
  geom_point(position = pd, size = 3) + 
  geom_line(alpha = 0.3, position = pd) + 
  labs(x = "% nalożenia CNV'ów", y = "%wariantów AthCNV znalezionych przez program") +
  geom_text(aes(label=round(`%AthCNV`, 1)), position=position_dodge(width=0.7), vjust=-0.4)



