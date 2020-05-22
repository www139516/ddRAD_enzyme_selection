#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/22 8:02
# @Author : Yuancong Wang
# @Site : 
# @File : cal_pic_from_vcf.py
# @Software: PyCharm
# @E-mail: wangyuancong@163.com
"""
Description: Calculate the pic value of every polymorphism site based on the vcf file
Output: The pic value of every site (dataframe)
Input: The vcf file
Other notes: The pic value is used as ref for further enzyme selection
"""


import argparse
from vcf_processor.pic_calculator import VcfCalculator


def main():
    parser = argparse.ArgumentParser(description="Calculate the pic values for the polymorphism sites.")
    parser.add_argument('-v', '--vcf', help='The path of the vcf file')
    parser.add_argument('-d', '--directory', help='The directory for the output file containing pic values', default='')
    args = parser.parse_args()
    pic_cal = VcfCalculator()
    pic_cal = pic_cal.fit(args.vcf, args.directory)
    pic_cal = pic_cal.sel_col_from_dataframe()
    pic_cal = pic_cal.add_pic()
    pic_cal.write_df_pic()


if __name__ == '__main__':
    main()



