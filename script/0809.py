# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 14:51:06 2024

@author: mikali
"""

'''
從pandasgwas套件獲得associations資料
'''

import pandas as pd

from pandasgwas.get_associations import get_associations_by_efo_id
associations = get_associations_by_efo_id('EFO_0003898')

#有時候會The request for https://www.ebi.ac.uk/gwas/rest/api/efoTraits/EFO_0003898/associations failed: response code was 500

# 按照A欄位合併表格(strongest_risk_alleles有426個，associations有其他資訊)
panda_A = pd.merge(associations.strongest_risk_alleles, associations.associations, on='associationId', how='outer')

# 使用 str.rsplit() 以最右邊的 '-' 分割 A 欄位
panda_A[['rsId', 'riskAllele']] = panda_A['riskAlleleName'].str.rsplit('-', n=1, expand=True)

#移除不想要的變數(重複)
column_names = panda_A.columns
print(column_names)

panda_A = panda_A.drop(columns=['riskFrequency_y'])
panda_A = panda_A.rename(columns={'riskFrequency_x': 'riskFrequency'})


panda_A_dropped = panda_A.drop(index=[361, 362, 363, 364])

'''
#	rsId
361	rs30187
362	rs116488202
364	rs116488202
363	rs10045403

'''

#loci是否有大於1的列
panda_A['list_length'] = panda_A['loci'].apply(len)
#沒錯，有四個
#是否要將它集中處理?
#好啊
# 提取特定行
panda_A_EX = panda_A.iloc[[361, 362, 363, 364]]
print(panda_A_EX)
panda_A_EX.to_csv('panda_A_EX.csv', index=False)

#既然去掉特定行了，我來把loci裡的欄位提取一下
panda_A_dropped.columns
panda_A_dropped = panda_A_dropped.rename(columns={'description': 'description_1'})


# 使用 apply() 函数，提取所有键
extracted_data = panda_A_dropped['loci'].apply(lambda row: pd.Series(row[0]) if isinstance(row, list) and len(row) > 0 else pd.Series())

# 将提取的值合并到原 DataFrame
panda_A_dropped = pd.concat([panda_A_dropped, extracted_data], axis=1)

import os

new_dir = "C:/Users/mikali/Desktop/ProbeDesignOutput/ProbeDesignOutput/data"
os.chdir(new_dir)

#index=False可以避免匯出行索引。
panda_A_dropped.to_csv('panda_A.csv', index=False)


'''
從pandasgwas套件獲得variants資料
'''

from pandasgwas.get_variants import get_variants_by_efo_id
snps = get_variants_by_efo_id('EFO_0003898')

panda_S = snps.variants

#locations 是否有大於1的列
panda_S['list_length'] = panda_S['locations'].apply(len)
#沒有，但倒是有O的，先刪除
'''
#	rsId
28	HLA-B*2705
29	HLA-B*2705
295	HLA-B*2702
332	chr18:14723700
370	rs67025039
371	imm_16_28525386
372	imm_16_28525386
373	HLA-B*2707
'''
panda_S_dropped = panda_S.drop(index=[28, 29, 295, 332, 370, 371, 372, 373])

# 使用 apply() 函数，提取所有键
extracted_data = panda_S_dropped['locations'].apply(lambda row: pd.Series(row[0]) if isinstance(row, list) and len(row) > 0 else pd.Series())

# 将提取的值合并到原 DataFrame
panda_S_dropped = pd.concat([panda_S_dropped, extracted_data], axis=1)


# 使用 apply() 提取字典中的键，转换为新的列
extracted_columns = panda_S_dropped['region'].apply(pd.Series)

# 将提取的列合并到原 DataFrame
panda_S_dropped = pd.concat([panda_S_dropped, extracted_columns], axis=1)



#接下來要處理panda_S其他欄
# 定義一個函數來過濾符合條件的genomicContexts的 dict
# 定義一個函數來篩選字典
def filter_dicts(dict_list):
    return [gene for gene in dict_list if gene['source'] == 'Ensembl' and gene['distance'] == 0]

# 使用 apply 方法來應用篩選函數
panda_S_dropped['filtered_genomicContexts'] = panda_S_dropped['genomicContexts'].apply(filter_dicts)

# 使用 apply 方法來篩選列表長度大於 1 的行
more_than_1_gene_panda_S = panda_S_dropped[panda_S_dropped['filtered_genomicContexts'].apply(len) > 1]
#到這裡會斷掉
#既然才16個，要不提取出來另外處理?
more_than_1_gene_panda_S.to_csv('more_than_1_gene_panda_S.csv', index=False)

#去除不需要的行
panda_S_dropped = panda_S_dropped.drop(columns=['locations', 'region'])
#region 重命名
panda_S_dropped = panda_S_dropped.rename(columns={'name': 'regionName'})

#提取其他基因資訊
panda_S_dropped.columns
# 使用 apply() 函数，提取所有键
extracted_data = panda_S_dropped['filtered_genomicContexts'].apply(lambda row: pd.Series(row[0]) if isinstance(row, list) and len(row) > 0 else pd.Series())

# 将提取的值合并到原 DataFrame
panda_S_dropped = pd.concat([panda_S_dropped, extracted_data], axis=1)

# 使用 apply() 函数，提取所有键
extracted_data = panda_S_dropped['gene'].apply(lambda row: pd.Series(row[0]) if isinstance(row, list) and len(row) > 0 else pd.Series())

# 将提取的值合并到原 DataFrame
panda_S_dropped = pd.concat([panda_S_dropped, extracted_data], axis=1)

# 使用 apply() 提取字典中的键，转换为新的列
extracted_columns = panda_S_dropped['gene'].apply(pd.Series)

# 将提取的列合并到原 DataFrame
panda_S_dropped = pd.concat([panda_S_dropped, extracted_columns], axis=1)
panda_S_dropped = panda_S_dropped.drop(columns=[0])

# 使用 apply() 函数，提取所有键
extracted_data = panda_S_dropped['entrezGeneIds'].apply(lambda row: pd.Series(row[0]) if isinstance(row, list) and len(row) > 0 else pd.Series())

# 将提取的值合并到原 DataFrame
panda_S_dropped = pd.concat([panda_S_dropped, extracted_data], axis=1)

# 使用 apply() 函数，提取所有键
extracted_data = panda_S_dropped['ensemblGeneIds'].apply(lambda row: pd.Series(row[0]) if isinstance(row, list) and len(row) > 0 else pd.Series())

# 将提取的值合并到原 DataFrame
panda_S_dropped = pd.concat([panda_S_dropped, extracted_data], axis=1)

#來去掉不必要的欄位吧

#去除不需要的行
#提取其他基因資訊
panda_S_dropped.columns

panda_S_dropped = panda_S_dropped.drop(columns=['list_length', 'filtered_genomicContexts', 'gene', 'entrezGeneIds', 'ensemblGeneIds'])

#搞定
panda_S_dropped.to_csv('panda_S.csv', index=False)


'''
結合來自ensembl網頁，GWAS網頁，pandasGWAS的資料
'''
new_dir = "C:/Users/mikali/Desktop/ProbeDesignOutput/ProbeDesignOutput/data"
os.chdir(new_dir)

#先從非特例的開始
import pandas as pd

ensembl = pd.read_csv('ensembl-export.csv')
panda_A = pd.read_csv('panda_A.csv')
panda_S = pd.read_csv('panda_S.csv')
gwas = pd.read_csv("gwas-association-downloaded_2024-08-09-EFO_0003898.tsv", sep="\t")

#先來刪除沒有在人類genome裡的
ensembl = ensembl[ensembl['Genomic location (strand)'].str.match(r'^[0-9]|^X|^Y')]
    #...pandas_A沒有 chromosome資訊。
    #...pandas_S之前就處理好了。
    # gwas 要先把 CHR_ID 5 x 6 抓出來當特例
gwas_EX = gwas[gwas['CHR_ID'] == '5 x 6']
gwas_EX.to_csv('gwas_EX.csv', index=False)
gwas_dropped = gwas[gwas['CHR_ID'] != '5 x 6']
gwas_dropped.to_csv('gwas.csv', index=False)

gwas_dropped = gwas_dropped.dropna(subset=['CHR_ID', 'CHR_POS'])
#沒有在人類genome裡的都刪除了
#處理變數以便對照
# 定義一個函數來處理觀測值
def process_observation(obs):
    parts = obs.split(':')
    chromosome = str(parts[0])
    
    start_end = parts[1].split('-')
    start = str(start_end[0])
    end = str(start_end[1][:-1])  # 去除最後的+或-
    
    sign = str(start_end[1][-1])  # 獲取最後的+或-
    
    # 明確指定列名
    return pd.Series([chromosome, start, end, sign], index=['ChromosomeName', 'Start', 'End', 'Strand'])
result = ensembl['Genomic location (strand)'].apply(process_observation)

# 使用 .loc 進行賦值
ensembl.loc[:, ['ChromosomeName', 'Start', 'End', 'Strand']] = result

#現在要從 panda_S, gwas, ensembl 找出unique rsId
#先提取不同表格所需項，再看有沒有重複
ensembl.columns
panda_S.columns
gwas_dropped.columns
# 提取特定變數並重新命名
extracted_ensembl = ensembl[['Name(s)', 'ChromosomeName', 'Start']].rename(columns={'Name(s)': 'dbSNP ID', 'ChromosomeName': 'Chromosome', 'Start': 'Position'})
extracted_panda_S = panda_S[['rsId', 'chromosomeName', 'chromosomePosition']].rename(columns={'rsId': 'dbSNP ID', 'chromosomeName': 'Chromosome', 'chromosomePosition': 'Position'})
extracted_gwas_dropped = gwas_dropped[['SNPS', 'CHR_ID', 'CHR_POS']].rename(columns={'SNPS': 'dbSNP ID', 'CHR_ID': 'Chromosome', 'CHR_POS': 'Position'})

# 合併這些變數
combined_df = pd.concat([extracted_ensembl, extracted_panda_S, extracted_gwas_dropped], ignore_index=True)

# 將變數A轉換為str型
combined_df['Position'] = combined_df['Position'].astype('str')

# 去掉重複的觀測值
unique_df = combined_df.drop_duplicates()

#輸出一下表格
new_dir = "C:/Users/mikali/Desktop/ProbeDesignOutput/ProbeDesignOutput/data"
os.chdir(new_dir)
unique_df.to_csv('unique_snp.csv', index=False)
#去R抓其他欄位
























