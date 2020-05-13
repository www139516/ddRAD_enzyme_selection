# coding: UTF-8
import click
from my_classes import Drawer


@click.command()
@click.option("--dir", help="The path of the directory that contains the output files using digest module")
def main(dir):
    bar_test = Drawer(dir)
    bar_test.get_the_summary_paths()
    bar_test.iter_plot_loci()
    bar_test.iter_frag_num()
    bar_test.iter_frag_cov()
    bar_test.merge_pdf()


if __name__ == '__main__':
    main()

