#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time : 2020/5/22 8:02
# @Author : Yuancong Wang
# @Site : 
# @File : cal_pic_from_vcf.py
# @Software: PyCharm
# @E-mail: wangyuancong@163.com
'''
Description: Calculate the pic value of every polymorphism site based on the vcf file
Output: The pic value of every site (dataframe)
Input: The vcf file
Other notes: The pic value is used as ref for further enzyme selection
'''

import pandas as pd
import os



