use RestrictionDigest;
# use List::MoreUtils qw(firstidx);
use List::Util qw(first);

############## 调用bowtie2 的相关变量信息
# my $ref_idx = "/home/wangyc/maizeGenome/v4/v4"; # bowtie2 index 文件位置
my $ref_idx = "/home/wangyc/maizeGenome/v3/v3";
my $output_directory = "/data2/wyc/e_enzyme/results_01/"; # 输出文件夹路径
my $file_res_digest_head = "seq_FB_BF_frags_in_range_v3.fa_by_"; # 输入fasta的文件名前缀

############## 创建目录用于存放保留文件
my $dir_res = $output_directory."results_keep/";

$dir_exist = -d "$output_directory";

if ($dir_exist) {

	print "Directory already exists\n";

} else {

	mkdir( $dir_res );
}



############## 调用 RestrictionDigest 的相关变量信息
my $ref_path = "/home/wangyc/maizeGenome/v3/v3.fa"; # The path of the reference maizeGenome
# my $output_directory = "/data2/wyc/e_enzyme/results/test_6"; # The directory where stores the result files
# my $gff_path = "/home/wangyc/maizeGenome/v4/v4.gff3"; # The path of the GFF file corresponds to the reference genome 

my $gff_path = "/home/wangyc/maizeGenome/v3/v3.gff3";

# bug: SbfI === PstI Illegal division by zero at /home/wangyc/anaconda3/envs/RestrictionDigest/lib/site_perl/5.26.2/RestrictionDigest.pm line 2966, <GEN19> line 34465177.
# bug: EcoRI === MlucI 

my $enzyme_list = "PstI EcoRI HindIII PvuII NlaIII MspI SbfI EcoRI AvaII MspI";

my @enzyme_list = split (/ +/, $enzyme_list);

my @enzyme_list_1 = @enzyme_list;

pop @enzyme_list_1;


print "Testing enzymes: @enzyme_list\n";

print "\n";

# for loop

foreach (@enzyme_list_1) {

	my $e_01 = $_;

	
	shift @enzyme_list;
	my @enzyme_list_2 = @enzyme_list;

	foreach(@enzyme_list_2) {



		my $e_02 = $_;

		print "============================================================$e_01,    $e_02================================================================\n";

		

		my $double_digest=RestrictionDigest::SingleItem::Double->new();

		$double_digest->add_ref(-reference => $ref_path); 

		# $double_digest->new_enzyme(-enzyme_name => "KpnI", -recognition_site => "GGTAC|C");

		# $double_digest->new_enzyme(-enzyme_name => "EcoRV", -recognition_site => "GAT|ATC");

		$double_digest->new_enzyme(-enzyme_name => "SbfI", -recognition_site => "CCTGCA|GG");

		$double_digest->new_enzyme(-enzyme_name => "MlucI", -recognition_site => "|AATT");

		$double_digest->new_enzyme(-enzyme_name => "HindIII", -recognition_site => "A|AGCTT");

		$double_digest->new_enzyme(-enzyme_name => "PstI", -recognition_site => "CTGCA|G");

		# $double_digest->new_enzyme(-enzyme_name => "XbaI", -recognition_site => "T|CTAGA");

		$double_digest->new_enzyme(-enzyme_name => "DraI", -recognition_site => "TTT|AAA");

		$double_digest->new_enzyme(-enzyme_name => "PvuII", -recognition_site => "CAG|CTG");

		$double_digest->new_enzyme(-enzyme_name => "NlaIII", -recognition_site => "CATG|");
		
		$double_digest->new_enzyme(-enzyme_name => "MspI", -recognition_site => "C|CGG");

		# $double_digest->new_enzyme(-enzyme_name => "XspI", -recognition_site => "C|TAG");

		$double_digest->new_enzyme(-enzyme_name => "AfaI", -recognition_site => "GT|AC");

		$double_digest->new_enzyme(-enzyme_name => "AvaII", -recognition_site => "GGWC|C");

		$double_digest->add_enzyme_pair(-front_enzyme =>$e_01, -behind_enzyme =>$e_02);

		$double_digest->change_range(-start=>250,-end=>350);

		$double_digest->change_lengths_distribution_parameters(-front=>200,-behind =>800,-step =>50);

		$double_digest->add_output_dir(-output_dir=> $output_directory);

		$double_digest->double_digest();

		# $double_digest->add_SNPs(-SNPs =>"path to the SNPs file");

		# $double_digest->_SNPs_at_fragments(-sequence_type =>"125SE", -sequence_end =>"front_enzyme");

		$double_digest->add_gff(-gff=>$gff_path);

		$double_digest->frags_in_range_coverage();


		# Call bowtie2 to obtain the SAM file

		my $file_res_digest_path = $output_directory.$file_res_digest_head.$e_01."_and_".$e_02; # 合并得到 输入的fasta路径和文件名
		my $out_sam_name = $output_directory.$e_01."_and_".$e_02.".sort.sam"; # 输出排序后的 SAM 文件
		print "Processing: ======= $file_res_digest_path\n";
		system ("bowtie2 -p 10 -x $ref_idx -fU $file_res_digest_path | samtools sort -O sam -@ 10 -o -> $out_sam_name");



		# 删除不需要的文件
		print "========Deleting the unnecessary files:\n";



		$dir = $output_directory."*"; 
		my @files = glob( $dir ); # 获得output目录下的所有文件名

		foreach (@files ) {

			#遍历每个文件, 如果是双酶切并且在长度范围内的文件就保留, 否则删除
			if ($_ =~ /digestion_summary/ || $_ =~ /coverage_of_FB_BF_frags_in_range/ || $_ =~ /position_FB_BF_frags_in_range/ || $_ =~ /seq_FB_BF_frags_in_range/ || $_ =~ /sam/) {

				system ("mv $_ $dir_res");
			} else {

				unlink ($_);
			}

		}

	};

}
