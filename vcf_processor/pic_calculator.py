#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/22 9:28
# @Author : Yuancong Wang
# @Site : 
# @File : pic_calculator.py
# @Software: PyCharm
# @E-mail: wangyuancong@163.com
"""
Description: Calculate pic values based on the input vcf file
Output: dataframe containing pic values
Input: path of vcf file
Other notes:
"""

import pandas as pd
import os
import re


class VcfCalculator:

    @staticmethod
    def extract_val_in_list(freq_list):
        val_list = []
        for words in freq_list:
            if re.findall(r'0.\d+', words):
                num = re.findall(r'0.\d+', words)[0]
                val_list.append(float(num))
            else:
                val_list.append(0)
        return val_list

    @staticmethod
    def cal_pic(pi_list):
        pic = 1 - sum([x * x for x in pi_list])
        return pic

    def __init__(self):
        self._in_vcf_fpath = None
        self._in_vcf_fname = None
        self._prefix_fname = None
        self._out_dpath = None
        self._out_fname = None
        self._out_fpath = None
        self._df_vcf = pd.DataFrame()

    def fit(self, in_vcf_pi, out_dpath):
        print('Initializing...')
        self._in_vcf_fpath = in_vcf_pi
        self._in_vcf_fname = os.path.basename(self._in_vcf_fpath)
        self._prefix_fname = self._in_vcf_fname.split('.')[0]
        if not out_dpath:
            self._out_dpath = os.path.dirname(self._in_vcf_fpath)
        else:
            self._out_dpath = out_dpath
        self._out_fname = self._prefix_fname + '_pic.tsv'
        self._out_fpath = os.path.join(self._out_dpath, self._out_fname)
        return self

    def sel_col_from_dataframe(self):
        print('Reading the vcf files.')
        df_vcf_pi = pd.read_csv(self._in_vcf_fpath, skiprows=1, names=['merged'])
        print('Selecting info from vcf files...')
        df_vcf_pi['chr'] = df_vcf_pi.merged.str.split('\t').apply(lambda x: x[0])
        df_vcf_pi['pos'] = df_vcf_pi.merged.str.split('\t').apply(lambda x: x[1])
        df_vcf_pi['n_allele'] = df_vcf_pi.merged.str.split('\t').apply(lambda x: x[2])
        df_vcf_pi['freq'] = df_vcf_pi.merged.str.split('\t').apply(lambda x: x[4:])
        df_vcf_pi = df_vcf_pi.loc[:, ['chr', 'pos', 'n_allele', 'freq']]
        self._df_vcf = df_vcf_pi
        return self



    def add_pic(self):
        print('Calculating pic values based on vcf...')
        self._df_vcf['vals'] = self._df_vcf.freq.apply(lambda x: VcfCalculator.extract_val_in_list(x))
        self._df_vcf['pic'] = self._df_vcf.vals.apply(lambda x: VcfCalculator.cal_pic(x))
        return self

    def write_df_pic(self):
        print('Writing dataframe with pic values...')
        self._df_vcf = self._df_vcf.loc[:, ['chr', 'pos', 'pic']]
        self._df_vcf.to_csv(self._out_fpath, sep='\t', index=None)

