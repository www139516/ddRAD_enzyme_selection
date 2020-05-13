import os


class Files:
    '''
    File class includes all the file related operations
    1. obtain the list of paths of all sam files in the target directory.
    2. screen out all the unique reads from the sam files
    '''
    @staticmethod
    def remove_exist_file(file_path):
        if os.path.exists(file_path):
            print('The target file exists.')
            print('The target file is removed.')
            os.remove(file_path)

    @staticmethod
    def give_file_name(file_path, file_type):
        name_insert = str(file_type)
        origin_name = os.path.basename(file_path)
        list_name = origin_name.split('.')
        pre_list_name = list_name[:-1]
        new_file_name = '.'.join(pre_list_name) + '.' + name_insert + '.' + list_name[-1]
        return new_file_name

    def __init__(self, dir_path, output_dir=os.getcwd()):
        self.dir_path = dir_path
        self.out_dir = output_dir

    def get_sam_file_paths(self):
        '''
        This function is used for obtaining all the '.sam' files under the 'dir_path'
        :return: sam file paths
        '''

        all_files = os.listdir(self.dir_path)
        list_sam_files_paths = list()

        for file in all_files:
            file_path = os.path.join(self.dir_path, file)
            if file.endswith('.sam'):
                list_sam_files_paths.append(file_path)
        return list_sam_files_paths

    def get_unique_reads_from_sam_files(self, sam_file_paths):
        '''
        screen out the reads that only map to unique positions
        :return: path of file containing unique reads
        '''
        uniq_sam_file_paths = list()
        for sam_file_path in sam_file_paths:
            print(sam_file_path)
            # sam_file_dir = os.path.split(sam_file_path)[0]
            sam_file_name = os.path.basename(sam_file_path)
            uniq_sam_file_name = Files.give_file_name(sam_file_name, 'uniq')
            uniq_sam_file_path = os.path.join(self.out_dir, uniq_sam_file_name)
            Files.remove_exist_file(uniq_sam_file_path)
            f_unique = open(uniq_sam_file_path, 'a')

            with open(sam_file_path, 'r') as my_sam_file:
                for line in my_sam_file.readlines():
                    is_genome = False
                    if line.strip()[0].isdigit() or 'chr' in line[0:9].lower():
                        is_genome = True
                    if is_genome and ('XS:i' not in line and '@' not in line):
                        f_unique.writelines(line)
            f_unique.close()
            uniq_sam_file_paths.append(uniq_sam_file_path)
        return uniq_sam_file_paths

