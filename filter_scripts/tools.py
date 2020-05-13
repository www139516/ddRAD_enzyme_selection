# coding:utf-8
'''
This file contains all the functions that are used in the 'main_pipe.py'

'''
import os 
import pandas as pd 
import re


# 获取指定目录下所有的 .sam 文件
def get_sam_files (dir_path):
    '''
    This function is used for obtaining all the '.sam' files under the 'dir_path'
    :return: sam file paths
    '''
    # dir_path: 传入存放SAM文件的文件夹路径
    # 获得sam文件路径
    # return: str, sam文件路径列表
    
    all_files = os.listdir(dir_path)
    
    list_sam_files = list()

    for file in all_files:
                       
        file_path = os.path.join(dir_path, file)

        if file.endswith('.sam'): 
       	
            list_sam_files.append(file_path)
            
    return list_sam_files



# 截取仅比对到一个位置的reads
def unique_sam (sam_file_path):
    # sam_file_path: 目标sam文件的路径
    # 过滤所有比对到超过一个位点的reads, 仅保留比对到唯一位点的reads
    # return: SAM
    

    # 获得SAM文件的目录路径和文件名
    sam_file_dir = os.path.split(sam_file_path)[0]

    sam_hap_dir = os.path.join(sam_file_dir, 'sam_hap/')

    if not os.path.exists(sam_hap_dir):

    	os.mkdir(sam_hap_dir)

    sam_file_name = os.path.split(sam_file_path)[1]
    
    # 根据SAM文件名得到unique_sam文件名
    sam_file_unique_name = sam_file_name[0:-4] + '.unique.' + 'sam'
    
    # 获得unique_sam文件名的路径
    sam_file_unique_path = os.path.join (sam_hap_dir,sam_file_unique_name)
    
    with open(sam_file_path, 'r') as f:

        
        # 判断nuique_sam文件是否存在,如果存在则删除
        if os.path.exists(sam_file_unique_path):

            print('The target file exists.')

            print('The target file is removed.')

            os.remove (sam_file_unique_path)

        
        # 将唯一比对到的reads信息写入新文件
        with open (sam_file_unique_path, 'a') as f_unique:
    
            for line in f.readlines():

                if (line[0].isdigit() or line[0:3] == 'chr' or line[0:3] == 'Chr'):

                    if ('XS:i' not in line and '@' not in line):

                        f_unique.writelines(line) 
                    
    return sam_file_unique_path




# 判断两个相邻位置差是否大于设定LD, 如果小于LD则算一个位点
# 设定 LD 大小, 默认1000 bp
def ld_sam(unique_sam_path, ld = 1000):
    # 计算相邻位点差距, 如大于LD则保留, 否则剔除
    # unique_sam_path 唯一比对上reads文件的存放路径
    # ld 默认1000bp
    # return: 剔除后reads的存放路径, str
    
    interval = ld 

    
    # 打开unique sam 的文件
    with open(unique_sam_path, 'r') as f:

        contents = f.readlines()

        num_lines = len(contents)
 
        i = 0
        
        ld_sam_dir = os.path.split(unique_sam_path)[0]

        unique_sam_name = os.path.split(unique_sam_path)[1] 
        
        unique_sam_ld_name = unique_sam_name[0:-4] + '.unique.ld_{}'.format(interval) + '.sam'

        unique_sam_ld_path = os.path.join(ld_sam_dir, unique_sam_ld_name) 


        if os.path.exists(unique_sam_ld_path):

                os.remove(unique_sam_ld_path)


        with open(unique_sam_ld_path,'a') as f_ld:

            if len(contents) > 0:
                f_ld.writelines(contents[0])


            for i in range(num_lines-1):

                j = i + 1

                read = contents[i].split('\t')

                chrom_num = read[2] 

                read_pos = read[3] 

                next_read = contents[j].split('\t')

                next_chrom_num = next_read[2]

                next_read_pos = next_read[3] 

                    
                # 计算染色体, 如果染色体位置不同, 则直接写入
                if (chrom_num != next_chrom_num): 

                    f_ld.writelines(contents[j])

                    continue

                # 计算间隔, 如果间隔大于指定间隔, 则写入新文件
                elif ((int(next_read_pos) - int(read_pos)) >= interval):

                    f_ld.writelines(contents[j]) 
                    
            return unique_sam_ld_path



# 计算每个reads比对的起始位置和终止位置
def cal_read_range(sam_ld_path):
    # sam_ld_path
    # 用于计算每条read在基因组上比对的范围
    # return sam_range路径


    with open(sam_ld_path, 'r') as f:

        contents = f.readlines()

        num_lines = len(contents)
        
        # 分隔ld.sam文件的路径
        sam_ld_dir = os.path.split(sam_ld_path)[0]

        sam_ld_name = os.path.split(sam_ld_path)[1]
        
        # 生成range文件的文件名
        sam_range_name = sam_ld_name[0:-4] + '_range.' + 'sam' 

        sam_range_path = os.path.join(sam_ld_dir,sam_range_name) 
        
        

    if os.path.exists(sam_range_path):

        os.remove(sam_range_path)


    # 打开文件range.sam句柄
    with open(sam_range_path,'a') as f_range:
        
        # 计算每一个read覆盖的范围(range)
        for i in range(num_lines):

            read = contents[i].split('\t')

            # 获得比对的起始位点
            read_pos = read[3]

            # 获得read长度
            read_len = len(read[9])

            # 比对read 终止位点为 起始位点 + 长度
            read_end_pos = int(read_pos) + int(read_len) 

            content_range = (contents[i].strip('\n') + '\t' + str(read_end_pos))

            # 写入终止位点
            f_range.writelines(content_range + '\n') 
                
    return sam_range_path



# 截取hapmap文件多态性的位点, 其它的信息删除
def get_the_poly_pos(hapmap_path, pic_val):
    '''
    Example:
    CHROM POS N_ALLELES N_CHR 1
    chr10 269 2 492 0.377025
    chr10 283 2 504 0.0871049
    chr10 284 2 506 0.0649365

    '''

    # hapmap_path: 传入一个hapmap的路径
    # 仅保留多态性在基因组中的位置信息
    # 将位置信息写入新文件
    # return: 新文件路径

    hmp_dir = os.path.dirname(hapmap_path)
    hmp_file_name = os.path.basename(hapmap_path)

    # 根据传入文件获得新文件名
    hmp_file_name_reduced = hmp_file_name[0:-4] + '_reduced_pic{}'.format(pic_val) + '.tab'
    hmp_file_path_reduced = os.path.join(hmp_dir, hmp_file_name_reduced)

    # 如果新文件存在, 先删除
    if os.path.exists(hmp_file_path_reduced):
        print('The file already exists')
    else:
        print('Generating the file containing polymorphism sites....')
        # col_names = ['chrom', 'pos', 'n_alleles', 'loc', 'pic_val']
        # all_pic_file = pd.read_csv(hapmap_path, sep='\t', names=col_names)
        # sel_pic_file = all_pic_file[all_pic_file.pic_val >= pic_val and all_pic_file.chrom.str.contains('chr')]
        # sel_pic_file.to_csv(hmp_file_path_reduced, sep='\t')
        with open(hapmap_path, 'r') as hmp_file:

            contents = hmp_file.readlines()

        num_lines = len(contents)

        hmp_dir = os.path.split(hapmap_path)[0]

        hmp_file_name = os.path.split(hapmap_path)[1]

        # 根据传入文件获得新文件名
        hmp_file_name_reduced = hmp_file_name[0:-4] + \
                                '_reduced_pic{}'.format(pic_val) + '.tab'

        hmp_file_path_reduced = os.path.join(hmp_dir, hmp_file_name_reduced)

        with open(hmp_file_path_reduced, 'a') as hmp_reduc:

            for line in range(num_lines):

                # list_line = contents[line].split('\t')[:5]
                list_line = re.split(r'\s+', contents[line])[:5]

                if (list_line[0].startswith('c') and float(list_line[-1]) >= pic_val):
                    content = '\t'.join(list_line) + '\n'

                    hmp_reduc.writelines(content)

    return hmp_file_path_reduced


# 比对hapmap, 仅保留有多态性的reads
def get_reads_fit_hmp(hapmap_reduced_path, read_range_path):
    # 读取过滤完的位点文件
    pd_hmp_reduced = pd.read_csv(hapmap_reduced_path, sep='\t', names=['chrom', 'pos', 'n_alle', 'n_chr', 'pic_val'])

    # 获取目标目录
    sam_with_hmp_dir = os.path.split(read_range_path)[0] 

    # 获取sam_range的原文件名
    sam_range_file = os.path.split(read_range_path)[1] 

 	# 根据原文件名获得新文件名和路径
    sam_ham_file = sam_range_file[0:-4] + '_hmp.' + 'sam'  

    sam_ham_path = os.path.join(sam_with_hmp_dir, sam_ham_file)

    with open(read_range_path, 'r') as sam_range:

        contents = sam_range.readlines()

        num_lines = len(contents)

    # 如果有同名文件存在, 先删除
    if os.path.exists(sam_ham_path):
        print('The file exists, deleting...')
        os.remove(sam_ham_path)

    with open(sam_ham_path, 'a') as f_hmp_sam:

        for i in range(num_lines):

            cur_read_chr = contents[i].split('\t')[2]

            cur_read_pos_start = contents[i].split('\t')[3]

            cur_read_pos_end = contents[i].split('\t')[-1]

            print('Searching reads ===> {}/'.format(i) + '{} reads'.format(num_lines))

            # 存储染色体数目信息
            cur_read_chr_digit = ''
        # 如果reads在基因组上则进行下一步操作
            if (cur_read_chr.isdigit()):

                cur_read_chr_digit = cur_read_chr

            elif (cur_read_chr.startswith('chr') and cur_read_chr[3:].isdigit()):

                cur_read_chr_digit = cur_read_chr[3:]

            elif (cur_read_chr.startswith('Chr') and cur_read_chr[3:].isdigit()):

                cur_read_chr_digit = cur_read_chr[3:]

            else:
                print('No match for this read.')
                continue

            if cur_read_chr_digit.isdigit():

            # 截取hammap_reduced中与目的read染色体相同的部分

                # pd_hmp_reduced_fit_chr = pd_hmp_reduced.loc[pd_hmp_reduced.chrom == 'chr{}'.format(cur_read_chr_digit),:]

                # 筛选位于reads覆盖位置范围内的位点
                start = int(cur_read_pos_start)
                end = int(cur_read_pos_end)

                pd_hmp_reduced_fit_chr_fil = pd_hmp_reduced[(pd_hmp_reduced.chrom == 'chr{}'.format(cur_read_chr_digit)) & (pd_hmp_reduced.pos >= start) & (pd_hmp_reduced.pos <= end)]

                num_rows = pd_hmp_reduced_fit_chr_fil.shape[0]

                if num_rows:
                    poly_loc = pd_hmp_reduced_fit_chr_fil.iloc[0, 1]
                    # poly_loc = ', '.join(str(poly_loc))
                    print(str(cur_read_pos_start) + ' ===== ' + str(cur_read_pos_end) + ' =====> ' + str(poly_loc))
                    line = contents[i].strip('\n') + '\t' + str(poly_loc) + '\n'
                    f_hmp_sam.writelines(line)
                else:
                    print('No match for this read')

            else:
                print('No match for this read.')
                continue

    return sam_ham_path
