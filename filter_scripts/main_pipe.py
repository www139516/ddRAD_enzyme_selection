# coding: UTF-8
from files import Files
import __future__
import tools
import os
import click
from ploter import Ploter

# tar_dir = '/data2/wyc/e_enzyme/results/results_keep/'

# hmp_file_path = '/data2/wyc/e_enzyme/genotype_97samples_2.hmp.txt'
@click.command()
@click.option("--directory", help="The path of the directory that contains the output files using digest module")
@click.option("--hmp", default='', help="The path of the hapmap file")
@click.option("--pic", default=0.3, help="The threshold of pic value")
@click.option('--hmpred', default="", help="The file path of prorcessed hmpfile with designated pic value")
@click.option("--ld", default=1000, help="The average length  of ld")
@click.option("--out", default='', help="The directory where you put the output files")
def pipe(directory, hmp, pic, hmpred, ld, out):

	my_files = Files(directory, out)
	my_sam_file_paths = my_files.get_sam_file_paths()
	my_uniq_file_paths = my_files.get_unique_reads_from_sam_files(my_sam_file_paths, out)


	for file in sam_files:

		print(file)

		sam_file_unique_path = tools.unique_sam(file)

		sam_file_range_path = tools.cal_read_range(sam_file_unique_path)

		if hmpred:

			hmp_reduced_path = hmpred

		else:

			hmp_reduced_path = tools.get_the_poly_pos(hmp, pic)

		sam_hmp_file_path = tools.get_reads_fit_hmp(hmp_reduced_path, sam_file_range_path)

		sam_file_ld_path = tools.ld_sam(sam_hmp_file_path, ld)

		my_rep_filter = Ploter(file, sam_file_unique_path, sam_hmp_file_path, sam_file_ld_path)

		my_rep_filter.plot_bar()

	my_rep_filter.merge_pdf()


if __name__ == '__main__':
	pipe()



