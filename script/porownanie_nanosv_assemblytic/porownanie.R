library(tidyverse)

nanosv <- read_csv("../wszystkie_warianty_odczyty.csv") %>% filter(name == "nanosv") %>% filter(!(CHROM == "mitochondria")) %>% filter(!(CHROM == "chloroplast")) %>%
  mutate(CHROM = as.numeric(CHROM), type = SVTYPE, program = name)
assemblytics <- read_csv("../asemblacja/wszystkie_warianty_asemblacja.csv") %>% filter(program == "assemblytics") %>% rename(END = STOP) %>% arrange(CHROM)

nanosv$END <- ifelse(is.na(nanosv$END), nanosv$START + nanosv$SVLEN, nanosv$END) 
##################################################
coverage = 0.001 
assemblytics_overlap_01 <- tibble()
nanosv_overlap_01 <- tibble()
for(index in 1:nrow(assemblytics))
{
  cnv <- assemblytics %>% slice(index)
  chrom <- nanosv %>% 
    filter( CHROM == cnv$CHROM, 
            cnv$START < START & cnv$END< END & (cnv$END - START >= 1 & cnv$END - START >= 1) | 
              cnv$START < START & cnv$END > END & (END - START >= 1 & END - START >= 1) |
              cnv$START > START & cnv$END < END & (cnv$END - cnv$START >= 1 & cnv$END - cnv$START >= 1) | 
              cnv$START > START & cnv$END > END & (END - cnv$START >= 1 & END - cnv$START >= 1))
  if(nrow(chrom) >= 1)
  {
    nanosv_overlap_01 <- bind_rows(chrom, nanosv_overlap_01)
    #assemblytics_overlap_01 <- bind_rows(cnv, assemblytics_overlap_01)
  }
}  

for(index in 1:nrow(nanosv))
{
  cnv <- nanosv %>% slice(index)
  chrom <- assemblytics %>% 
    filter( CHROM == cnv$CHROM, 
            cnv$START < START & cnv$END< END & (cnv$END - START >= 1 & cnv$END - START >= 1) | 
              cnv$START < START & cnv$END > END & (END - START >= 1 & END - START >= 1) |
              cnv$START > START & cnv$END < END & (cnv$END - cnv$START >= 1 & cnv$END - cnv$START >= 1) | 
              cnv$START > START & cnv$END > END & (END - cnv$START >= 1 & END - cnv$START >= 1))
  if(nrow(chrom) >= 1)
  {
    #nanosv_overlap_01 <- bind_rows(chrom, nanosv_overlap_01)
    assemblytics_overlap_01 <- bind_rows(chrom, assemblytics_overlap_01)
  }
} 
assemblytics_overlap_01 <- assemblytics_overlap_01 %>% unique
nanosv_overlap_01 <- nanosv_overlap_01 %>% unique

warianty_any <- bind_rows(nanosv_overlap_01, assemblytics_overlap_01)
warianty_any$pokrycie <- ">=1bp"


##############################################################################################
coverage = 0.01 
assemblytics_overlap_1 <- tibble()
nanosv_overlap_1 <- tibble()
for(index in 1:nrow(assemblytics))
{
  cnv <- assemblytics %>% slice(index)
  chrom <- nanosv %>% 
    filter( CHROM == cnv$CHROM, 
            cnv$START < START & cnv$END< END & (cnv$END - START >= coverage*SVLEN & cnv$END - START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START < START & cnv$END > END & (END - START >= coverage*SVLEN & END - START >= abs(coverage*cnv$SVLEN)) |
              cnv$START > START & cnv$END < END & (cnv$END - cnv$START >= coverage*SVLEN & cnv$END - cnv$START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START > START & cnv$END > END & (END - cnv$START >= coverage*SVLEN & END - cnv$START >= abs(coverage*cnv$SVLEN)))
  if(nrow(chrom) >= 1)
  {
    nanosv_overlap_1 <- bind_rows(chrom, nanosv_overlap_1)
  }
}  

for(index in 1:nrow(nanosv))
{
  cnv <- nanosv %>% slice(index)
  chrom <- assemblytics %>% 
    filter( CHROM == cnv$CHROM, 
            cnv$START < START & cnv$END< END & (cnv$END - START >= coverage*SVLEN & cnv$END - START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START < START & cnv$END > END & (END - START >= coverage*SVLEN & END - START >= abs(coverage*cnv$SVLEN)) |
              cnv$START > START & cnv$END < END & (cnv$END - cnv$START >= coverage*SVLEN & cnv$END - cnv$START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START > START & cnv$END > END & (END - cnv$START >= coverage*SVLEN & END - cnv$START >= abs(coverage*cnv$SVLEN)))
  if(nrow(chrom) >= 1)
  {
    #nanosv_overlap_01 <- bind_rows(chrom, nanosv_overlap_01)
    assemblytics_overlap_1 <- bind_rows(chrom, assemblytics_overlap_1)
  }
} 
assemblytics_overlap_1 <- assemblytics_overlap_1 %>% unique
nanosv_overlap_1 <- nanosv_overlap_1 %>% unique
warianty_1 <- bind_rows(assemblytics_overlap_1, nanosv_overlap_1)
warianty_1$pokrycie <- ">=1%"


##############################################################################################

coverage = 0.1 
assemblytics_overlap_10 <- tibble()
nanosv_overlap_10 <- tibble()
for(index in 1:nrow(assemblytics))
{
  cnv <- assemblytics %>% slice(index)
  chrom <- nanosv %>% 
    filter( CHROM == cnv$CHROM, 
            cnv$START < START & cnv$END< END & (cnv$END - START >= coverage*SVLEN & cnv$END - START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START < START & cnv$END > END & (END - START >= coverage*SVLEN & END - START >= abs(coverage*cnv$SVLEN)) |
              cnv$START > START & cnv$END < END & (cnv$END - cnv$START >= coverage*SVLEN & cnv$END - cnv$START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START > START & cnv$END > END & (END - cnv$START >= coverage*SVLEN & END - cnv$START >= abs(coverage*cnv$SVLEN)))
  if(nrow(chrom) >= 1)
  {
    nanosv_overlap_10 <- bind_rows(chrom, nanosv_overlap_10)
  }
}

for(index in 1:nrow(nanosv))
{
  cnv <- nanosv %>% slice(index)
  chrom <- assemblytics %>% 
    filter( CHROM == cnv$CHROM, 
            cnv$START < START & cnv$END< END & (cnv$END - START >= coverage*SVLEN & cnv$END - START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START < START & cnv$END > END & (END - START >= coverage*SVLEN & END - START >= abs(coverage*cnv$SVLEN)) |
              cnv$START > START & cnv$END < END & (cnv$END - cnv$START >= coverage*SVLEN & cnv$END - cnv$START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START > START & cnv$END > END & (END - cnv$START >= coverage*SVLEN & END - cnv$START >= abs(coverage*cnv$SVLEN)))
  if(nrow(chrom) >= 1)
  {
    assemblytics_overlap_10 <- bind_rows(chrom, assemblytics_overlap_10)
  }
}

assemblytics_overlap_10 <- assemblytics_overlap_10 %>% unique
nanosv_overlap_10 <- nanosv_overlap_10%>% unique
warianty_10 <- bind_rows(assemblytics_overlap_10, nanosv_overlap_10)
warianty_10$pokrycie <- ">=10%"


##############################################################################################

coverage = 0.9 
assemblytics_overlap_90 <- tibble()
nanosv_overlap_90 <- tibble()
for(index in 1:nrow(assemblytics))
{
  cnv <- assemblytics %>% slice(index)
  chrom <- nanosv %>% 
    filter( CHROM == cnv$CHROM, 
            cnv$START < START & cnv$END< END & (cnv$END - START >= coverage*SVLEN & cnv$END - START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START < START & cnv$END > END & (END - START >= coverage*SVLEN & END - START >= abs(coverage*cnv$SVLEN)) |
              cnv$START > START & cnv$END < END & (cnv$END - cnv$START >= coverage*SVLEN & cnv$END - cnv$START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START > START & cnv$END > END & (END - cnv$START >= coverage*SVLEN & END - cnv$START >= abs(coverage*cnv$SVLEN)))
  if(nrow(chrom) >= 1)
  {
    nanosv_overlap_90 <- bind_rows(chrom, nanosv_overlap_90)
  }
}

for(index in 1:nrow(nanosv))
{
  cnv <- nanosv %>% slice(index)
  chrom <- assemblytics %>% 
    filter( CHROM == cnv$CHROM, 
            cnv$START < START & cnv$END< END & (cnv$END - START >= coverage*SVLEN & cnv$END - START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START < START & cnv$END > END & (END - START >= coverage*SVLEN & END - START >= abs(coverage*cnv$SVLEN)) |
              cnv$START > START & cnv$END < END & (cnv$END - cnv$START >= coverage*SVLEN & cnv$END - cnv$START >= abs(coverage*cnv$SVLEN)) | 
              cnv$START > START & cnv$END > END & (END - cnv$START >= coverage*SVLEN & END - cnv$START >= abs(coverage*cnv$SVLEN)))
  if(nrow(chrom) >= 1)
  {
    assemblytics_overlap_90 <- bind_rows(chrom, assemblytics_overlap_90)
  }
}

assemblytics_overlap_90 <- assemblytics_overlap_90 %>% unique
nanosv_overlap_90 <- nanosv_overlap_90 %>% unique
warianty_90 <- bind_rows(assemblytics_overlap_90, nanosv_overlap_90)
warianty_90$pokrycie <- ">=90%"


###########################################
warianty <- bind_rows(warianty_any, warianty_1, warianty_90, warianty_10)
warianty$pokrycie <- factor
warianty %>% ggplot(aes(x = program, fill = length)) + 
  geom_bar(position = "fill") + 
  labs(x = "program", y = "Ułamek wariantów") + 
  scale_fill_discrete(name = "Długość") + 
  facet_wrap(~pokrycie)


warianty %>% ggplot(aes(x =program, fill = type)) + 
  geom_bar(position = "fill") + 
  labs(x = "program", y = "Ułamek wariantów") + 
  scale_fill_discrete(name = "Typ wariantu") + 
  facet_wrap(~pokrycie)


warianty %>% ggplot(aes(x = type, fill = length)) + 
  geom_bar(position = "fill")  + 
  facet_wrap(~program)
  

########################################################################

ath <- read_csv("ath_cnv.csv", skip = 4)
ath <- ath %>% filter(`Included in AthCNV dataset?` == "YES") %>% select(CNV_ID, Chr, Start, Stop, size, `Included in AthCNV dataset?`)
ath$Chr <- factor(gsub("Chr", "", ath$Chr))

ath <- ath %>%  mutate(CHROM = Chr, START = )
