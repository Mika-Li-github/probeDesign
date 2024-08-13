# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:35:55 2024

@author: mikali
"""

import os

# 打印当前工作目录
print("当前工作目录:", os.getcwd())

# 更改工作目录
new_directory = 'C:/Users/mikali/Desktop/ProbeDesignOutput/ProbeDesignOutput/data'  # 替换为你想要的目录路径
os.chdir(new_directory)

import pandas as pd

seq362 = pd.read_csv('allele362.csv', index_col=0)
#這次就算沒有MAF也先做

#第一步，切「一半」
seq362[['left_part', 'right_part']] = seq362['snp'].str.extract(r'(.+?)%[^%]+%(.+)')

#第二步，把allele接至left_part的最後
seq362.columns
# 將A欄位字符接在B欄位字串後面，形成C欄位
#小心-

seq362['minor allele (Alternative)'] = seq362['minor allele (Alternative)'].replace('-', '')

seq362['left_ref'] = seq362['left_part'] + seq362['ref allele']
seq362['left_min'] = seq362['left_part'] + seq362['minor allele (Alternative)']
seq362['ref_right'] = seq362['ref allele'] + seq362['right_part']
seq362['min_right'] = seq362['minor allele (Alternative)'] + seq362['right_part']

#要把ref_right和min_right順序反過來，且要AT互換，CG互換
# 定义一个函数来反转字符串，要反過來才可以算Tm
def reverse_string(s):
    return s[::-1]

# 应用函数到变量A
seq362['right_ref'] = seq362['ref_right'].apply(reverse_string)
seq362['right_min'] = seq362['min_right'].apply(reverse_string)


#接下來是ATCG互換
def swap_nucleotides(dna_sequence):
    # 创建一个映射字典，用于互换
    swap_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '-':''}
    
    # 使用列表推导式进行互换
    swapped_sequence = ''.join(swap_dict[nucleotide] for nucleotide in dna_sequence)
    
    return swapped_sequence


# 进行互换
seq362['right_ref_conv'] = seq362['right_ref'].apply(swap_nucleotides)
seq362['right_min_conv'] = seq362['right_min'].apply(swap_nucleotides)

#可以用同一套邏輯去計算Tm值了
#left_ref, left_min, right_ref_conv, right_min_conv

#先把left_ref叫出來看看
test_seq = seq362['left_ref'].iloc[0]
#801字元，如預期的。

#想要的輸出:Fwd_ref, Tm_Fr
from Bio.SeqUtils import MeltingTemp as mt


def Tm_Seq(sequences):
    minimun_diff = float('inf')
    for i in range(len(sequences)-16, -1, -1):#第一次取就16字元，取到0(0-1=-1)，也就是整列，每次迭代都減1
            sub_seq = sequences[i:] #子序列
            current_tm_Wallace = mt.Tm_Wallace(sub_seq) #子序列Tm值
            current_diff = abs(current_tm_Wallace - 62.5) #子序列Tm值與62.5的距離
            if 61 <current_tm_Wallace < 63.9: #如果在範圍內
                if current_diff < minimun_diff: #而且當前距離小於最小距離
                    minimun_diff = current_diff #則記錄當前距離
                    best_sub = sub_seq #現在的序列就是最好的序列
                    best_tm = current_tm_Wallace #現在的Tm值就是最好的Tm值                        
                    results = [best_sub, best_tm]
    return results

Tm_Seq(seq362['left_ref'].iloc[0])
mt.Tm_Wallace(seq362['left_ref'].iloc[0][798:])
seq362 = seq362.copy()
seq362[['Fwd_ref', 'Tm_Fr']] = seq362['left_ref'].apply(Tm_Seq).apply(pd.Series)

seq362[['Fwd_min', 'Tm_Fm']] = seq362['left_min'].apply(Tm_Seq).apply(pd.Series)
seq362[['Rev_ref', 'Tm_Rr']] = seq362['right_ref_conv'].apply(Tm_Seq).apply(pd.Series)
seq362[['Rev_min', 'Tm_Rm']] = seq362['right_min_conv'].apply(Tm_Seq).apply(pd.Series)

#要來驗證
validating = seq362["Rev_ref"].iloc[0]
source = seq362["right_ref_conv"].iloc[0]
more1 = mt.Tm_Wallace(source[len(source)-len(validating)-1:])
Tm = mt.Tm_Wallace(source[len(source)-len(validating):])
less1 = mt.Tm_Wallace(source[len(source)-len(validating)+1:])

result = abs(Tm - 62.5) < abs(less1 - 62.5) and abs(Tm - 62.5) < abs(more1 - 62.5)
if result:
    print("A和B都为True")
else:
    print("A和B有一个或全部为False")

#那麼，現在要來建立primer表了
import pandas as pd

seq362.columns
# 创建一个字典primer_Fwd_ref
primer_Fwd_ref = {
    'rs_id': seq362['dbSNP ID'],
    'primer_id': seq362['dbSNP ID'].astype(str) + '_Fwd' + seq362['ref allele'],
    'Tm': seq362['Tm_Fr'],
    'GC content(%)': (seq362['Fwd_ref'].apply(lambda x: x.count('G') + x.count('C')) / seq362['Fwd_ref'].apply(len)) * 100,
    'genotype_label':seq362['ref allele'],
    'primer_len': seq362['Fwd_ref'].apply(len),
    'SNP_chr': seq362['Chromosome'],
    'SNP_position': seq362['Position'],
    'primer_sequence':seq362['Fwd_ref']
}

# 从字典创建DataFrame
primertable_Fwd_ref = pd.DataFrame(primer_Fwd_ref)

#驗證GC content(%)
(seq362['Fwd_ref'].apply(lambda x: x.count('G') + x.count('C')) / seq362['Fwd_ref'].apply(len)) * 100
(17/44)*100


mt.Tm_NN('TATGATTTTATTTAAACAAATAATGAAATATTACTTTTTGGATAATACTATTTTTATGTTTACATTATTTTAGAGACTTAAAAAAATCACTGAAATATTTACCATGATTAAA')
# 62.49590180554668
mt.Tm_GC('TATGATTTTATTTAAACAAATAATGAAATATTACTTTTTGGATAATACTATTTTTATGTTTACATTATTTTAGAGACTTAAAAAAATCACTGAAATATTTACCATGATTAAA')
# 61.13504492912075
mt.Tm_Wallace('TATGATTTTATTTAAACAAATAATGAAATATTACTTTTTGGATAATACTATTTTTATGTTTACATTATTTTAGAGACTTAAAAAAATCACTGAAATATTTACCATGATTAAA')
#如果用Wallace，此行260

'''
####################
'''
seq362.columns
# 创建一个字典primer_Fwd_min
primer_Fwd_min = {
    'rs_id': seq362['dbSNP ID'],
    'primer_id': seq362['dbSNP ID'].astype(str) + '_Fwd' + seq362['minor allele (Alternative)'],
    'Tm': seq362['Tm_Fm'],
    'GC content(%)': (seq362['Fwd_min'].apply(lambda x: x.count('G') + x.count('C')) / seq362['Fwd_min'].apply(len)) * 100,
    'genotype_label':seq362['minor allele (Alternative)'],
    'primer_len': seq362['Fwd_min'].apply(len),
    'SNP_chr': seq362['Chromosome'],
    'SNP_position': seq362['Position'],
    'primer_sequence':seq362['Fwd_min']
}

# 从字典创建DataFrame
primertable_Fwd_min = pd.DataFrame(primer_Fwd_min)
'''
####################
.apply(swap_nucleotides)
'''
seq362.columns
# 创建一个字典 primer_Rev_ref
primer_Rev_ref = {
    'rs_id': seq362['dbSNP ID'],
    'primer_id': seq362['dbSNP ID'].astype(str) + '_Rev' + seq362['ref allele'].apply(swap_nucleotides),
    'Tm': seq362['Tm_Rr'],
    'GC content(%)': (seq362['Rev_ref'].apply(lambda x: x.count('G') + x.count('C')) / seq362['Rev_ref'].apply(len)) * 100,
    'genotype_label':seq362['ref allele'].apply(swap_nucleotides),
    'primer_len': seq362['Rev_ref'].apply(len),
    'SNP_chr': seq362['Chromosome'],
    'SNP_position': seq362['Position'],
    'primer_sequence':seq362['Rev_ref']
}

# 从字典创建DataFrame
primertable_Rev_ref = pd.DataFrame(primer_Rev_ref)
'''
####################
.apply(swap_nucleotides)
'''
seq362.columns
# 创建一个字典 primer_Rev_min
primer_Rev_min = {
    'rs_id': seq362['dbSNP ID'],
    'primer_id': seq362['dbSNP ID'].astype(str) + '_Rev' + seq362['minor allele (Alternative)'].apply(swap_nucleotides),
    'Tm': seq362['Tm_Rm'],
    'GC content(%)': (seq362['Rev_min'].apply(lambda x: x.count('G') + x.count('C')) / seq362['Rev_min'].apply(len)) * 100,
    'genotype_label':seq362['minor allele (Alternative)'].apply(swap_nucleotides),
    'primer_len': seq362['Rev_min'].apply(len),
    'SNP_chr': seq362['Chromosome'],
    'SNP_position': seq362['Position'],
    'primer_sequence':seq362['Rev_min']
}

# 从字典创建DataFrame
primertable_Rev_min = pd.DataFrame(primer_Rev_min)
'''
####################
'''

#合併表格
primer = pd.concat([primertable_Fwd_ref, primertable_Fwd_min, primertable_Rev_ref, primertable_Rev_min], ignore_index=True)
primer.to_csv('primer.csv', index=False, encoding='utf-8')

"""
##########################################################
以下做probe
##########################################################

"""

#要先把 Ref 和 Min 分別對應
seq362['ref_len'] = seq362['Fwd_ref'].apply(len) + seq362['Rev_ref'].apply(len) - seq362['ref allele'].apply(len)
seq362['min_len'] = seq362['Fwd_min'].apply(len) + seq362['Rev_min'].apply(len) - seq362['minor allele (Alternative)'].apply(len)


import numpy as np

seq362['probe_which'] = np.where(seq362['ref_len'] >= seq362['min_len'], 'ref_len', 'min_len')

seq362['probe_len'] = np.where(seq362['ref_len'] >= seq362['min_len'], seq362['ref_len'], seq362['min_len'])

#probe前半已經在Fwd_ref或Fwd_min了
#後半是要..right_part 取Rev_ref或Rev_min的長度-min_allele長度
seq362['Fwd_ref_Rev'] = seq362["Fwd_ref"] + seq362.apply(
    lambda row: row['right_part'][:len(row['Rev_ref']) - len(row['ref allele'])] if len(row['Rev_ref']) > 1 else '',
    axis=1
)
seq362['Fwd_min_Rev'] = seq362["Fwd_min"] + seq362.apply(
    lambda row: row['right_part'][:len(row['Rev_min']) - len(row['minor allele (Alternative)'])] if len(row['Rev_min']) > 1 else '',
    axis=1
)

seq362['probe_sequence'] = np.where(seq362['probe_which'] == 'ref_len', seq362['Fwd_ref_Rev'], seq362['Fwd_min_Rev'])

# 创建一个字典 primer_Rev_min
probe = {
    'rs_id': seq362['dbSNP ID'],
    'Tm': seq362['probe_sequence'].apply(lambda x: mt.Tm_Wallace(x)), 
    'GC content(%)': (seq362['probe_sequence'].apply(lambda x: x.count('G') + x.count('C')) / seq362['probe_sequence'].apply(len)) * 100,
    'probe_len': seq362['probe_sequence'].apply(len),
    'SNP_chr': seq362['Chromosome'],
    'SNP_position': seq362['Position'],
    'probe_sequence':seq362['probe_sequence']
}

# 从字典创建DataFrame
probe = pd.DataFrame(probe)






#沒錯，要改wallace




		
























