#!/usr/bin/env python
# -*- coding: utf-8 -*
'''
2019/11/25
The test file for the e_enzyme project
version: 1.0
Auther: Yuancong Wang
'''

# import classes
#
# my_dir = classes.EDigestionDir('I:/实验相关/E_enzyme/demo_results/demo_01_v3_b73/results_keep')
#
# sam_files = my_dir.get_the_sam_files()
#
# res_dir = my_dir.mk_res_dir()
#
# print(sam_files)
#
# for file_path in sam_files:
#
#     my_unique_sam = classes.SamFile(file_path)
#     my_unique_sam_path = my_unique_sam.get_unique_reads_range(res_dir)
#
# print(my_unique_sam_path)

import re

ss = 'CM018594.1 Zea mays cultivar K0326Y ecotype Quality Protein Maize chromosome 12, whole genome shotgun sequence'
aa = re.findall(r'chromosome (\d+)', ss)
print(aa)


