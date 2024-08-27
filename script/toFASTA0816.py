# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 14:55:33 2024

@author: mikali
"""

#先把序列檔案讀入
import os
new_directory = 'C:/Users/mikali/Desktop/ProbeDesignOutput/ProbeDesignOutput/data'  # 替換為你想要的目錄路徑
os.chdir(new_directory)
print("目前工作目錄:", os.getcwd()) #列印目前工作目錄
import pandas as pd
primer_GC = pd.read_csv('primer_GC.csv')
probe_GC = pd.read_csv('probe_GC.csv')

# 指定輸出 FASTA 文件的名稱
output_file = 'primer_GC.fasta'

# 將 DataFrame 轉換為 FASTA 格式並寫入文件
with open(output_file, 'w') as f:
    for index, row in primer_GC.iterrows():
        f.write(f">{row['primer_id']}\n")  # 寫入標題行
        f.write(f"{row['primer_sequence']}\n")  # 寫入序列行
        
print(f"FASTA 文件已保存為 {output_file}")

# 指定輸出 FASTA 文件的名稱
output_file = 'probe_GC.fasta'

# 將 DataFrame 轉換為 FASTA 格式並寫入文件
with open(output_file, 'w') as f:
    for index, row in probe_GC.iterrows():
        f.write(f">{row['rs_id']}\n")  # 寫入標題行
        f.write(f"{row['probe_sequence']}\n")  # 寫入序列行
        
print(f"FASTA 文件已保存為 {output_file}")