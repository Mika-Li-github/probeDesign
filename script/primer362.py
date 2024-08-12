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

