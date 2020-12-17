library(tidyverse)

lumpy<- read_csv("lumpy.csv", col_types = cols(.default = "d", CHROM = "c", SVTYPE = "c", ALT = "c"))
sniffles<- read_csv("sniffles.csv", col_types = cols(.default = "d", CHROM = "c", SVTYPE = "c", ALT = "c"))
sniffles_minimap<- read_csv("sniffles_minimap.csv", col_types = cols(.default = "d", CHROM = "c", SVTYPE = "c", ALT = "c"))
nanosv<- read_csv("nanosv.csv", col_types = cols(.default = "d", CHROM = "c", SVTYPE = "c", ALT = "c"))

lumpy$name = "lumpy"
sniffles$name = "sniffles"
sniffles_minimap$name = "sniffles_minimap"
nanosv$name = "nanosv"

warianty <- bind_rows(lumpy, sniffles, sniffles_minimap, nanosv) %>% filter(!is.na(SVLEN))
warianty$length <- ""
warianty$length <- ifelse(abs(warianty$SVLEN) >= 500, ">=500bp", 
                    ifelse( abs(warianty$SVLEN) < 500 & abs(warianty$SVLEN) >= 50, "50-499bp", "< 50bp" ))
warianty$length <- factor(warianty$length, levels = c("< 50bp", "50-499bp", ">=500bp" ))

war <- warianty %>% count(name, sort = TRUE) %>% rename(program = name)
war %>% ggplot(aes(x = program, y = n, fill = program)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25)
  
 
warianty %>% ggplot(aes(x = name, fill = length)) + 
  geom_bar(position = "fill") + 
  labs(x = "program", y = "Ułamek wariantów") + 
  scale_fill_discrete(name = "Długość")

warianty %>% ggplot(aes(x =name, fill = SVTYPE)) + 
  geom_bar(position = "fill") + 
  labs(x = "program", y = "Ułamek wariantów") + 
  scale_fill_discrete(name = "Typ wariantu")

warianty_500 <- warianty %>% filter(length == ">=500bp")

war_500 <- warianty_500 %>% count(name, sort = TRUE) %>% rename(program = name)
war_500 %>% ggplot(aes(x = program, y = n, fill = program)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25)

warianty_500 %>% ggplot(aes(x =name, fill = SVTYPE)) + 
  geom_bar(position = "fill") + 
  labs(x = "program", y = "Ułamek wariantów") + 
  scale_fill_discrete(name = "Typ wariantu")

ath <- read_csv("ath_cnv.csv", skip = 4)
ath <- ath %>% filter(`Included in AthCNV dataset?` == "YES") %>% select(CNV_ID, Chr, Start, Stop, size, `Included in AthCNV dataset?`)
ath$Chr <- factor(gsub("Chr", "", ath$Chr))


warianty_500$END <- ifelse(is.na(warianty_500$END), warianty_500$START + warianty_500$SVLEN, warianty_500$END) 
#warianty_500$SVLEN <- abs(warianty_500$SVLEN)

coverage <- 0.2
overlap_02 <- ath
overlap_02$lumpy <- 0
overlap_02$sniffles_minimap <- 0
overlap_02$sniffles <- 0
overlap_02$nanosv <- 0


for(index in 1:nrow(warianty_500))
{
  cnv <- warianty_500 %>% slice(index)
  chrom <- ath %>% 
    filter( Chr == cnv$CHROM, 
      cnv$START < Start & cnv$END< Stop & (cnv$END - Start >= coverage*size & cnv$END - Start >= abs(coverage*cnv$SVLEN)) | 
            cnv$START < Start & cnv$END > Stop & (Stop - Start >= coverage*size & Stop - Start >= abs(coverage*cnv$SVLEN)) |
            cnv$START > Start & cnv$END < Stop & (cnv$END - cnv$START >= coverage*size & cnv$END - cnv$START >= abs(coverage*cnv$SVLEN)) | 
            cnv$START > Start & cnv$END > Stop & (Stop - cnv$START >= coverage*size & Stop - cnv$START >= abs(coverage*cnv$SVLEN)))
  program <- cnv$name
  if(program == "lumpy") {overlap_02$lumpy[overlap_02$CNV_ID %in% chrom$CNV_ID] <- 1}
  else if(program == "sniffles") {overlap_02$sniffles[overlap_02$CNV_ID %in% chrom$CNV_ID] <- 1 }
  else if(program == "sniffles_minimap") {overlap_02$sniffles_minimap[overlap_02$CNV_ID %in% chrom$CNV_ID] <- 1 }
  else if(program == "nanosv") { overlap_02$nanosv[overlap_02$CNV_ID %in% chrom$CNV_ID] <- 1 }
}  

coverage <- 0.5
overlap_05 <- ath
overlap_05$lumpy <- 0
overlap_05$sniffles_minimap <- 0
overlap_05$sniffles <- 0
overlap_05$nanosv <- 0


for(index in 1:nrow(warianty_500))
{
  cnv <- warianty_500 %>% slice(index)
  chrom <- ath %>% 
    filter( Chr == cnv$CHROM, 
            cnv$START < Start & cnv$END< Stop & (cnv$END - Start >= coverage*size & cnv$END - Start >= abs(coverage*cnv$SVLEN)) | 
              cnv$START < Start & cnv$END > Stop & (Stop - Start >= coverage*size & Stop - Start >= abs(coverage*cnv$SVLEN)) |
              cnv$START > Start & cnv$END < Stop & (cnv$END - cnv$START >= coverage*size & cnv$END - cnv$START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START > Start & cnv$END > Stop & (Stop - cnv$START >= coverage*size & Stop - cnv$START >= abs(coverage*cnv$SVLEN)))
  program <- cnv$name
  if(program == "lumpy") {overlap_05$lumpy[overlap_05$CNV_ID %in% chrom$CNV_ID] <- 1}
  else if(program == "sniffles") {overlap_05$sniffles[overlap_05$CNV_ID %in% chrom$CNV_ID] <- 1 }
  else if(program == "sniffles_minimap") {overlap_05$sniffles_minimap[overlap_05$CNV_ID %in% chrom$CNV_ID] <- 1 }
  else if(program == "nanosv") { overlap_05$nanosv[overlap_05$CNV_ID %in% chrom$CNV_ID] <- 1 }
}  


coverage <- 0.8
overlap_08 <- ath
overlap_08$lumpy <- 0
overlap_08$sniffles_minimap <- 0
overlap_08$sniffles <- 0
overlap_08$nanosv <- 0


for(index in 1:nrow(warianty_500))
{
  cnv <- warianty_500 %>% slice(index)
  chrom <- ath %>% 
    filter( Chr == cnv$CHROM, 
            cnv$START < Start & cnv$END< Stop & (cnv$END - Start >= coverage*size & cnv$END - Start >= abs(coverage*cnv$SVLEN)) | 
              cnv$START < Start & cnv$END > Stop & (Stop - Start >= coverage*size & Stop - Start >= abs(coverage*cnv$SVLEN)) |
              cnv$START > Start & cnv$END < Stop & (cnv$END - cnv$START >= coverage*size & cnv$END - cnv$START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START > Start & cnv$END > Stop & (Stop - cnv$START >= coverage*size & Stop - cnv$START >= abs(coverage*cnv$SVLEN)))
  program <- cnv$name
  if(program == "lumpy") {overlap_08$lumpy[overlap_08$CNV_ID %in% chrom$CNV_ID] <- 1}
  else if(program == "sniffles") {overlap_08$sniffles[overlap_08$CNV_ID %in% chrom$CNV_ID] <- 1 }
  else if(program == "sniffles_minimap") {overlap_08$sniffles_minimap[overlap_08$CNV_ID %in% chrom$CNV_ID] <- 1 }
  else if(program == "nanosv") { overlap_08$nanosv[overlap_08$CNV_ID %in% chrom$CNV_ID] <- 1 }
}  

write_csv(overlap_02, "overlap_20.csv")
write_csv(overlap_05, "overlap_50.csv")
write_csv(overlap_08, "overlap_80.csv")

overlap_02$pokrycie = "20%"
overlap_05$pokrycie = "50%"
overlap_08$pokrycie = "80%"


overlap <- bind_rows(overlap_02, overlap_05, overlap_08)

write_csv(overlap, "odczyty_AthCNV_detekcja.csv")

overlaps <- overlap %>% 
  group_by(pokrycie) %>% 
  mutate(lumpy = sum(lumpy)*100/sum(nrow(ath)), sniffles = sum(sniffles)*100/sum(nrow(ath)), sniffles_minimap = sum(sniffles_minimap)*100/sum(nrow(ath)), nanosv = sum(nanosv)*100/sum(nrow(ath))) %>%
  select(lumpy,sniffles_minimap, sniffles, nanosv) %>% 
  unique %>% 
  pivot_longer(., c("lumpy", "sniffles_minimap", "sniffles", "nanosv"), names_to = "program", values_to = "%AthCNV")

write_csv(overlaps, "pokrycie_AthCNV_odczyty.csv")

pd <- position_dodge(0.1)
ggplot(overlaps, aes(pokrycie, `%AthCNV`, color = program, group = program))  +
  geom_point(position = pd, size = 3) + 
  geom_line(alpha = 0.3, position = pd) + 
  labs(x = "% nalożenia CNV'ów", y = "%wariantów AthCNV znalezionych przez program") +
  geom_text(aes(label=round(`%AthCNV`, 1)), position=position_dodge(width=0.7), vjust=-0.4)

# 
# 
# ggplot(overlaps, aes(program, `%AthCNV`, fill = program))  +
#   geom_bar(stat = "identity") + 
#   labs(x = "% nalożenia CNV'ów") + 
#   facet_wrap(~pokrycie)

odleglosc <- overlap %>% 
  pivot_longer(., c("lumpy", "sniffles_minimap", "sniffles", "nanosv"), names_to = "program", values_to = "detekcja") %>% 
  filter(detekcja == 1)

odleglosc %>% ggplot(aes(program, size, fill = program)) +
  geom_violin(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(odleglosc$size, c(0.1, 0.9)))


