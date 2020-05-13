# coding: UTF-8
import re
import click


@click.command()
@click.option("--f_ori", help="The path of the original fasta file")
@click.option("--f_tar", help="The path of the target fasta file")
def main(f_ori, f_tar):
    file_path = f_ori
    target_file_path = f_tar
    with open(file_path, 'r') as f_ori:

        with open(target_file_path, 'a') as f_tar:
            for line in f_ori.readlines():
                if 'chromosome' in line:
                    chr_num = re.findall(r'chromosome (\d+)', line)[0]
                    chrom = '>chr{}\n'.format(chr_num)
                    f_tar.writelines(chrom)
                else:
                    f_tar.writelines(line)


if __name__ == '__main__':
    main()



