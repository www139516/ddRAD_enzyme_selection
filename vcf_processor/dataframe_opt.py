#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/22 18:12
# @Author : Yuancong Wang
# @Site : 
# @File : dataframe_opt.py
# @Software: PyCharm
# @E-mail: wangyuancong@163.com
'''
Description:
Output: 
Input:
Other notes:
'''


import pandas as pd
import os
import re


class DataframeOpt:
    def __init__(self):
        self._df_lst = []

    def fit(self, in_dpath, end_pattern):
        all_fnames = os.listdir(in_dpath)
        df_fnames = []
        df_fpaths = []
        for file in all_fnames:
            if file.endswith(end_pattern):
                df_fnames.append(file)
        df_fpaths = [os.path.join(in_dpath, fname) for fname in df_fnames]
        self._df_lst = [pd.read_csv(df_fpath, sep='\t', header=True) for df_fpath in df_fpaths]

    def combine_vertical(self):
        pass
