# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:31:03 2024

@author: mikali
"""

#from Bio.SeqUtils import MeltingTemp as mt

import os
import pandas as pd

os.chdir('C:/Users/mikali/Desktop/ProbeDesignOutput/ProbeDesignOutput/data')
seq386allele = pd.read_csv('seq386allele.csv')


# 使用最近邻法计算Tm
#tm_nn = mt.Tm_NN(seq386allele['Variant sequences'].values[0])

#seq386allele['tm_nn'] = seq386allele['Variant sequences'].apply(mt.Tm_NN)

#我現在要來拆兩對了
#我覺得要補ref，minor, 還有把沒序列的去掉
#biopython似乎做得到?
#可去來源看看
#rs116046827

from Bio import Entrez
import xml.etree.ElementTree as ET
import re

#b5eca0fdc6438dd7367e13025cf005045d08	
Entrez.email = "akispen01@gmail.com"
Entrez.api_key = "b5eca0fdc6438dd7367e13025cf005045d08"


# fetch info
#rsid = 'rs12769205'  # multiple alt bases
rsid = "rs116046827"  # single alt base
handle = Entrez.efetch(db="SNP", id=rsid, retmode="text")
# result is in xml format
xml_str = handle.readline().strip()

# parse the xml string
myroot = ET.fromstring(xml_str)

# parsing from DOCSUM
# SEQ=[G/A/...]  The first character is the reference base,
# followed by alternative bases
docsum_txt = myroot.find('DOCSUM').text
print(docsum_txt)
# regex to punch out A/G/... in memory parentheses
ptn = re.compile(r'SEQ=\[(.+)\]')
bases = ptn.search(docsum_txt).groups()[0].split('/')  # assume match

# print in "ref>alt" format
print(', '.join([bases[0] + '>' + base for base in bases[1:]]))
#好感動，是ok的


#先把缺alt allele的取出來
seq386allele.columns
nan_seq386allele = seq386allele[seq386allele['minor allele (Alternative)'].isnull()]


'''
#還有一位長這樣，不算。
	dbSNP ID	Chromosome	Position	snp	Extracted	ref allele	minor allele (Alternative)
41	rs11363316	16	30159695	GTGTTGGCCGGGCTGGTCTCCAGCTCCTGACCTTGAGTGATCTGCCCGTCTCAGCCTCCTGAGGTGCTGGGATTGCAGACGGAGTCTCGCTCACTCAATGCTCAATGTTGCCCAGGCTGGAGTGCAGTGGCGTGATCTCGGCTCGCTAAAACCTCCACCTCCCAGCCGCCTGCCTTGGCCTCCTAAAGTGCTAAGATTACAGCCTCTGCCCGGCCGCCACCCCGTCTGGGAAGTGAGGAGCGCTTCTGCCCGGCCGCCACCCCGTCTGTGCTGGGTGTGGTGGTGGGCGCCTGTAGTCCCAGCTTTTTGGAAGGCTGAGGCAGGAGAATCGCTTGAACCCAGGAGGCGGAGGTTGCTGTGAGCTGAGATTACACCACTGCACTCCAGCCTGGGCGACAGAGCAAGACTCCATCTTAAAAAAAAAAAAAGGCCGGGCACAGTGGCTCACGCCTCTAATCCCAGCACTTTGGGAGGACAAGGCAGGTGGATCACGAGGTCAGGAGATTGAGATCATCCTGGCTAACACGATGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCCTGGTGGCAGGCACCTGTAGTACCAGCTACTTGGGAGGCTGAGGCAGGAGAATGACGTGAACCTGGGAGGCGGAGCTGGAAGTGAGCCAAGATCGTGCCACTGCACTCCAGCCTGGGCGACAAAGCGAAACTCCATCTCAGAAAAAAAAAAAAAAAAGTAAGAAAAGAAAAAAAGAAATGTAAAATAGTAAACAGAAATATCTCCCAGCTCTTAAAATTAGGCAATTCTG%TT/T%GGAAAGGACATAACTCTGAATAGAGGTTATAAGGAAGAGTTTGGATTTAACAGGAGAGTGCCAAGAAGCTCTTCCAGTTTACAAAGCAGGAGAGTATCATTAAGCCTGTGTTTTAGAAGATGACTGCATTCTGTGATTCTGAAATAGCTATATATAGGTTGCTCGTGAAATCTGTGCTGCTAGTCTTCTTTTTTTTTTTTTTGAGACAGAATCTTGCTCTGTCCCCCAAGCTGGAGTGCAGTGGCACGATCTCAGCTCACTGCAACCTCCGCTTCCTGGGTTCAAGCAATTCTCCTGTCTTGGTCTCCTGAGTAGCTGGGATTACAGGAGTGCACCACCACGCCCAGCTAATTTTTGTATTTTTAATAGAGAAAGGGTTTCCCTGTGTTGACCAGGCTGGTCTTAAACTCCTGACCTCAGGTGATCTGCCCGCCTCGGCCTTCCAAAGTGCTGGGATTACAGGCGTGAGCCACTGTGCCCGGCGTCGCTAGTCTTCTCTTACTGAAGTGAAACTCACCTCTATTCTTTTTTTTTTTCTTTTTCTGTTTTTGAGACAGAGTCTCAGTCTGTCACCCAGGCTGGAGTGCAGTAGCTCGATCTCGGCTCACTGCAGCCTTCACCTCCTGGGTTCAAGTGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGATTATAGGAGTGCATCACCACGCCCAGATAATTTTTTTTTTGTATTTTTAGTACAGGTGGGGTTTCACCATGTTGGCCAGGCTGGTTTTGAACTCCTGACCTCAAATGATCCACCGGCCTCGGCCTCCCAA	TT/T	T	-
他算完整
#試著找找看
# 假設我們知道要提取的觀測值的索引
indices_to_extract = [41]  # 例如提取索引 1 和 3 的觀測值

# 使用索引選擇觀測值
extracted_rows = seq386allele.loc[indices_to_extract]

# 創建一個新的 DataFrame，並將提取的觀測值新增到其中
nan_seq386allele = pd.concat([nan_seq386allele, extracted_rows], ignore_index=True)
'''

# 定義一個函數來獲取 SNP 信息
def fetch_snp_info(rsid):
    handle = Entrez.efetch(db="SNP", id=rsid, retmode="text")
    xml_str = handle.readline().strip()

    # 解析 XML 字符串
    myroot = ET.fromstring(xml_str)

    # 從 DOCSUM 中解析
    docsum_elem = myroot.find('DOCSUM')
    if docsum_elem is not None:
        docsum_txt = docsum_elem.text
        ptn = re.compile(r'SEQ=\[(.+)\]')
        match = ptn.search(docsum_txt)

        if match:
            bases = match.groups()[0].split('/')  # 假設匹配
            # 返回 "ref>alt" 格式的字符串
            return ', '.join([bases[0] + '>' + base for base in bases[1:]])
    return None  # 如果沒有匹配或 DOCSUM 不存在，返回 None

# 將 SNP 信息添加到 DataFrame 中的新列
#nan_seq386allele['snp_info'] = nan_seq386allele['dbSNP ID'].apply(fetch_snp_info)
#Bad Request
nan_seq386allele = nan_seq386allele.copy()
nan_seq386allele['snp_info'] = nan_seq386allele['dbSNP ID'].apply(fetch_snp_info)
# 計算變數 B 中逗號的數量並加 1，然後存入新變數 'comma_count_B'
nan_seq386allele['count_snp_info'] = nan_seq386allele['snp_info'].str.count(',') + 1

#這玩意兒挺有參考價值
nan_seq386allele.to_csv('nan_seq386allele47.csv')


#20個有2 allele，4個有3 allele
#先取出只有一個的
nan_seq386allele_23 = nan_seq386allele[~nan_seq386allele['snp_info'].str.contains('[,]', na=False)]

nan_seq386allele_23 = nan_seq386allele_23.copy()
nan_seq386allele_23[['ref allele', 'minor allele (Alternative)']] = nan_seq386allele_23['snp_info'].str.split('>', expand=True)

nan_seq386allele_23 = nan_seq386allele_23.drop(['snp_info', 'count_snp_info'], axis=1)

# 獲取小表格中變數A的唯一值
small_a_values = nan_seq386allele_23['dbSNP ID'].unique()

# 過濾大表格，去掉變數A在小表格中出現的行
filtered_large_table = seq386allele[~seq386allele['dbSNP ID'].isin(small_a_values)]

# 將小表格的行追加到過濾後的大表格
final_table = pd.concat([filtered_large_table, nan_seq386allele_23], ignore_index=True)

#可以先dropna(24個)了
final_table_dropped_subset = final_table.dropna(subset=['ref allele', 'minor allele (Alternative)'])

final_table_dropped_subset.to_csv('allele362.csv')
#去設計probe和primer吧































