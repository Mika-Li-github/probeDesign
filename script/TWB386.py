# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 12:32:01 2024

@author: mikali
"""

import os

os.chdir('C:/Users/mikali/Desktop/gitfamily/TWB/TWBv2.0_SNPs位點相關資訊')

import pandas as pd

TWBSNP位點 = pd.read_csv('TWBv2.0 SNPs位點相關資訊.csv')

TWBSNP位點.rename(columns={'dbSNP RS ID': 'dbSNP ID'}, inplace=True)

# 選取需要的變量
print(TWBSNP位點.columns)
TWBSNP = TWBSNP位點[['dbSNP ID', 'Chromosome', 'Physical Position', 'Ref Allele', 'Alt Allele', "Minor Allele", 'Minor Allele Frequency']]

os.chdir('C:/Users/mikali/Desktop/ProbeDesignOutput/ProbeDesignOutput/data')
#先前獲得序列的386SNP
seq386 = pd.read_csv('seq386.csv')
print(seq386.columns)

filtered_TWBSNP = TWBSNP[TWBSNP['dbSNP ID'].isin(seq386['dbSNP ID'])]
#238，比之前多兩筆

filtered_TWBSNP = filtered_TWBSNP.drop_duplicates()
#事實上只有205個在TWBSNP

seq386['Position'] = seq386['Position'].astype(str)
filtered_TWBSNP['Physical Position'] = filtered_TWBSNP['Physical Position'].astype(str)
filtered_TWBSNP['Chromosome'] = filtered_TWBSNP['Chromosome'].astype(str)

# 選取需要比較的欄位,並重新命名
seq386_cols = ['dbSNP ID', 'Chromosome', 'Position']
filtered_TWBSNP_cols = ['dbSNP ID', 'Chromosome', 'Physical Position']

filtered_TWBSNP_1 = filtered_TWBSNP[filtered_TWBSNP_cols].rename(columns={'Physical Position':'Position'})
seq386_1 = seq386[seq386_cols]

# 以主變數為依據合併兩個 DataFrame
merged_df = pd.merge(filtered_TWBSNP_1, seq386_1, on='dbSNP ID', how='inner')

# 比較小表格的 var1 和 var2 與大表格的 other_var1 和 other_var2
comparison_result = (merged_df['Chromosome_x'] == merged_df['Chromosome_y']) & (merged_df['Position_x'] == merged_df['Position_y'])

# 輸出結果
if comparison_result.all():
    print("所有變數均相等")
else:
    print("存在不相等的變數")
#竟然有不相等的變數

# 比較小表格的 var1 和 var2 與大表格的 other_var1 和 other_var2
not_equal_Chromosome = merged_df[merged_df['Chromosome_x'] != merged_df['Chromosome_y']]
not_equal_Position = merged_df[merged_df['Position_x'] != merged_df['Position_y']]

# 提取不相等的行
not_equal_rows = pd.concat([not_equal_Chromosome, not_equal_Position]).drop_duplicates()

# 輸出不相等的行
print("不相等的項目:")
print(not_equal_rows)

'''
現在要來看差補的資訊
'''
# 讀取TSV文件
tsv_file_path = 'C:/Users/mikali/Desktop/gitfamily/TWB/TWBv2.0_差補位點資訊.v2/TWB2.hg38.infoscore.v2.txt'  # 替換為你的TSV文件路徑
TWB2_hg38_infoscore_v2 = pd.read_csv(tsv_file_path, sep='\t')
#這次一樣先全386都找看看
TWB2_hg38_infoscore_v2.columns

filtered_TWB2_hg38_infoscore_v2 = TWB2_hg38_infoscore_v2[TWB2_hg38_infoscore_v2['[3]ID'].isin(seq386['dbSNP ID'])]
filtered_TWB2_hg38_infoscore_v2 = filtered_TWB2_hg38_infoscore_v2.drop_duplicates()
#這次315個

#那有沒有不在那205裡面的?
filtered_TWB2_hg38_infoscore_v2_2 = filtered_TWB2_hg38_infoscore_v2[~filtered_TWB2_hg38_infoscore_v2['[3]ID'].isin(filtered_TWBSNP['dbSNP ID'])]
#136個

diff_from_TWB2_136 = filtered_TWB2_hg38_infoscore_v2_2[['[3]ID', '# [1]CHROM', '[2]POS', '[4]REF', '[5]ALT']]

#保險起見去掉重複值
diff_from_TWB2_136 = diff_from_TWB2_136.drop_duplicates()
#呵呵，沒有

new_dir = "C:/Users/mikali/Desktop/ProbeDesignOutput/ProbeDesignOutput/data"
os.chdir(new_dir)

diff_from_TWB2_136.to_csv('diff_from_TWB2_136.csv', index=False)

#386中還有哪些不在TWB
filtered_seq386 = seq386[~seq386['dbSNP ID'].isin(filtered_TWBSNP['dbSNP ID'])]
#四大天王中，不在TWBv2.0 SNPs位點相關資訊的，181個

filtered_seq386 = filtered_seq386[~filtered_seq386['dbSNP ID'].isin(diff_from_TWB2_136['[3]ID'])]
#45個
filtered_seq386.to_csv('TWB2missing_SNPS45.csv', index=False)

#來合併阿
#首先變數要一致
diff_from_TWB2_136.rename(columns={'[3]ID': 'dbSNP ID', '# [1]CHROM': 'Chromosome', '[2]POS': 'Position', '[4]REF': 'ref allele', '[5]ALT': 'minor allele (Alternative)'}, inplace=True)

diff_from_TWB2_136['Chromosome'] = diff_from_TWB2_136['Chromosome'].str.replace(r'[chr]', '', regex=True)
diff_from_TWB2_136['Position'] = diff_from_TWB2_136['Position'].astype(str)
#檢查位置是否相等
# 選取需要比較的欄位,並重新命名
seq386_cols = ['dbSNP ID', 'Chromosome', 'Position']
diff_from_TWB2_136_cols = ['dbSNP ID', 'Chromosome', 'Position']

diff_from_TWB2_136_1 = diff_from_TWB2_136[diff_from_TWB2_136_cols]
seq386_1 = seq386[seq386_cols]


# 以主變數為依據合併兩個 DataFrame
merged_df = pd.merge(diff_from_TWB2_136_1, seq386_1, on='dbSNP ID', how='inner')

# 比較小表格的 var1 和 var2 與大表格的 other_var1 和 other_var2
comparison_result = (merged_df['Chromosome_x'] == merged_df['Chromosome_y']) & (merged_df['Position_x'] == merged_df['Position_y'])

# 輸出結果
if comparison_result.all():
    print("所有變數均相等")
else:
    print("存在不相等的變數")
#竟然有不相等的變數

# 比較小表格的 var1 和 var2 與大表格的 other_var1 和 other_var2
not_equal_Chromosome = merged_df[merged_df['Chromosome_x'] != merged_df['Chromosome_y']]
not_equal_Position = merged_df[merged_df['Position_x'] != merged_df['Position_y']]

# 提取不相等的行
not_equal_rows = pd.concat([not_equal_Chromosome, not_equal_Position]).drop_duplicates()

# 輸出不相等的行
print("不相等的項目:")
print(not_equal_rows)
# 創建 filtered_TWBSNP 的副本
filtered_TWBSNP = filtered_TWBSNP.copy()

# 現在可以安全地重命名列而不會出現警告
filtered_TWBSNP.rename(columns={
    'rsId': 'dbSNP ID',
    'Physical Position': 'Position',
    'Ref Allele': 'ref allele',
    'Alt Allele': 'minor allele (Alternative)',
    'Minor Allele': 'minor allele',
    'Minor Allele Frequency': 'MAF'
}, inplace=True)
#和有序列的一起合併，但對不上位置的那兩個要先去掉
filtered_TWBSNP = filtered_TWBSNP[filtered_TWBSNP['dbSNP ID'] != "rs2066847"]
diff_from_TWB2_136 = diff_from_TWB2_136[diff_from_TWB2_136['dbSNP ID'] != "rs5837881"]


# 先合併大表格和中表格
# 正確的方式是將列名放在列表中
filtered_TWBSNP_1 = filtered_TWBSNP[['dbSNP ID', 'ref allele', 'minor allele (Alternative)', 'minor allele', 'MAF']]

merged_df = pd.merge(seq386, filtered_TWBSNP_1, on='dbSNP ID', how='outer')

# 再合併結果與小表格
diff_from_TWB2_136_1 = diff_from_TWB2_136[['dbSNP ID', 'ref allele', 'minor allele (Alternative)']]

final_df = pd.merge(merged_df, diff_from_TWB2_136_1, on='dbSNP ID', how='outer')
final_df['ref allele_x'] = final_df['ref allele_x'].fillna(final_df['ref allele_y'])
final_df['minor allele (Alternative)_x'] = final_df['minor allele (Alternative)_x'].fillna(final_df['minor allele (Alternative)_y'])

final_df = final_df.rename(columns={'ref allele_x': 'ref allele', 'minor allele (Alternative)_x': 'minor allele (Alternative)'})

final_df = final_df.drop(['ref allele_y', 'minor allele (Alternative)_y'], axis=1)

final_df.to_csv('seq386allele.csv', index=False)

























