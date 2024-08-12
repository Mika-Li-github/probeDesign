library(tidyverse)
library(biomaRt)
library(BiocManager)
setwd("C:/Users/mikali/Desktop/ProbeDesignOutput/ProbeDesignOutput/data")

mart = useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
ensembl <- useEnsembl(biomart = "snps")
ensembl <- useDataset(dataset = "hsapiens_snp", mart = ensembl)
snps386 <- read_csv("unique_snp.csv")

#獲得序列
snpmart<- useMart("ENSEMBL_MART_SNP",dataset = "hsapiens_snp")
snp_sequence <- getBM(attributes = c("refsnp_id", 
                                     "snp"), 
                      filters = c("snp_filter", "chr_name", "start", "upstream_flank", "downstream_flank"), 
                      values = list(snps386$`dbSNP ID`, snps386$Chromosome, snps386$Position, 800, 800), 
                      mart = snpmart,
                      checkFilters = FALSE 
)
#Chromosome和Position是必要的
#哪些是沒有的，挑出來
# 使用 dplyr 套件篩選 var1 為 NA 的觀測值
merged_df <- merge(snps386, snp_sequence, by.x = "dbSNP ID", by.y = "refsnp_id", all.x = TRUE)


na_snp_sequence <- merged_df %>% filter(is.na(snp))
na_snp_sequence


na_snp_sequence <- tibble(
  `dbSNP ID` = c("rs5837881", "rs72867732", "rs28533662", "rs11363316", "rs2066847"),
  Chromosome = c("2", "6", "6", "16", "16"),
  Position = c("203843041", "31392512", "31358156", "30159695", "50729868")
)
print(na_snp_sequence)
na_sequence <- getBM(attributes = c("refsnp_id", 
                                     "snp"), 
                      filters = c("snp_filter", "chr_name", "start", "upstream_flank", "downstream_flank"), 
                      values = list(na_snp_sequence$`dbSNP ID`, na_snp_sequence$Chromosome, na_snp_sequence$Position, 800, 800), 
                      mart = snpmart,
                      checkFilters = FALSE 
)
filtered_snps386 <- snps386 %>% filter(!`dbSNP ID` %in% c("rs11306716", "rs115879499", "rs117486637", "rs34670647", "rs5743293"))
filtered_snps386 <- rbind(filtered_snps386, na_snp_sequence)
combined_seq <- rbind(snp_sequence, na_sequence)

merged_df <- merge(filtered_snps386, combined_seq, by.x = "dbSNP ID", by.y = "refsnp_id", all.x = TRUE)
#補完序列了，好耶

#檢查一下序列形式
# 使用正则表达式提取"%%"之间的字符
merged_df$Extracted <- str_extract(merged_df$`snp`, "%(.*?)%")

# 去掉提取结果中的"%"
merged_df$Extracted <- gsub("%", "", merged_df$Extracted)

write_csv(merged_df, "seq386.csv")
#然後就去設計探針吧