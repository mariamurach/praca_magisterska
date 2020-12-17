library(tidyverse)

SV_colnames <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE1')

ssplit <- function(s, split = '='){
  unlist(strsplit(s, split = split))
}

# note, capital letters just respect original naming conventions of the VCF file
getSVTYPE <- function(info){
  ssplit(grep("SVTYPE", info, value = T))[2]
}

getSVLEN <- function(info){
  SVLEN <- ssplit(grep("SVLEN", info, value = T))
  ifelse(length(SVLEN) == 0, NA, as.numeric(SVLEN[2]))
}

getEND <- function(info){
  END <- ssplit(grep("END", info, value = T))
  ifelse(length(END) == 0, NA, as.numeric(END[2]))
}

load_sv <- function(file){
  vcf_sv_table <- read.table(file, stringsAsFactors = F)
  colnames(vcf_sv_table) <- SV_colnames
  # possible filtering
  # vcf_sv_table <- vcf_sv_table[vcf_sv_table$FILTER == 'PASS',]
  vcf_sv_table_info <- strsplit(vcf_sv_table$INFO, ';')
  vcf_sv_table$SVTYPE <- unlist(lapply(vcf_sv_table_info, getSVTYPE))
  vcf_sv_table$SVLEN <- unlist(lapply(vcf_sv_table_info, getSVLEN))
  vcf_sv_table$END <- unlist(lapply(vcf_sv_table_info, getEND))
  return(vcf_sv_table)
}

lumpy <- load_sv('../cnv/lumpy.vcf') %>% as_tibble() %>% select(CHROM, ID, START = POS, END, SVLEN,  ALT, SVTYPE)
sniffles <- load_sv('../cnv/sniffles.vcf') %>% as_tibble() %>% select(CHROM, ID, START = POS, END, SVLEN,  ALT, SVTYPE)
sniffles_minimap <- load_sv('../cnv/sniffles_minimap.vcf') %>% as_tibble() %>% select(CHROM, ID, START = POS, END, SVLEN,  ALT, SVTYPE)
nanosv <- load_sv('../cnv/nanosv.vcf') %>% as_tibble() %>% select(CHROM, ID, START = POS, END, SVLEN,  ALT, SVTYPE)
write_csv(lumpy, "lumpy.csv")
write_csv(sniffles, "sniffles.csv")
write_csv(sniffles_minimap, "sniffles_minimap.csv")
write_csv(nanosv, "nanosv.csv")
lumpy_variants <- as.data.frame(table(lumpy$SVTYPE)) %>%
  as_tibble %>% 
  rename(lumpy = Freq)
lumpy_variants <- lumpy_variants %>%
  gather(key = program, value = value, 2:ncol(lumpy_variants)) %>% 
  spread(key = names(lumpy_variants)[1], value = 'value')

sniffles_variants <- as.data.frame(table(sniffles$SVTYPE)) %>%
  as_tibble %>% 
  rename(sniffles = Freq)
sniffles_variants <- sniffles_variants %>%
  gather(key = program, value = value, 2:ncol(sniffles_variants)) %>% 
  spread(key = names(sniffles_variants)[1], value = 'value')

sniffles_minimap_variants <- as.data.frame(table(sniffles_minimap$SVTYPE)) %>%
  as_tibble %>% 
  rename(sniffles_minimap = Freq)
sniffles_minimap_variants <- sniffles_minimap_variants %>%
  gather(key = program, value = value, 2:ncol(sniffles_minimap_variants)) %>% 
  spread(key = names(sniffles_minimap_variants)[1], value = 'value')

nanosv_variants <- as.data.frame(table(nanosv$SVTYPE)) %>%
  as_tibble %>% 
  rename(nanosv = Freq)
nanosv_variants <- nanosv_variants %>%
  gather(key = program, value = value, 2:ncol(nanosv_variants)) %>% 
  spread(key = names(nanosv_variants)[1], value = 'value')

variants <- bind_rows(nanosv_variants, sniffles_variants, sniffles_minimap_variants, lumpy_variants)
variants_longer <- variants %>% pivot_longer(., cols = c("BND", "DUP", "INS", "DEL", "INV", "INVDUP"), names_to = "wariant", values_to = "n")

write_csv(variants, "variants.csv")

variants_longer %>% 
  ggplot(aes(program, n, fill = program)) +
  geom_bar(stat="identity", width=1) + 
  facet_wrap(vars(wariant))



  