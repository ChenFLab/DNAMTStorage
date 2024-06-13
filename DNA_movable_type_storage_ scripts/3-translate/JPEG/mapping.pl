#!/usr/bin/perl
use strict;
use warnings;
my $seq_dir = $ARGV[0];
my $mw_dir=$ARGV[1];
my $mapping_res_dir=$ARGV[2];
open(SEQ, "<$seq_dir") or die "file.txt 文件无法打开, $!";
open(mapping_res,">$mapping_res_dir");
my $i=0;
my $seq_name;
my $mapping;
while(<SEQ>){
    if($i % 2==0){
        $seq_name=$_;
    }
    elsif($i % 2==1){
        open temp_file,">$seq_dir.$i.file";
        print  temp_file ">$i\n$_";
        close temp_file;
   ##建立一个序列名的文件，然后建立索引，比对完之后删掉所有索引文件
        system("bwa index $seq_dir.$i.file>nohup.txt");
   ###进行比对
	    system("bwa aln $seq_dir.$i.file $mw_dir > $seq_dir.$i.sai");
	    system(" bwa samse $seq_dir.$i.file $seq_dir.$i.sai $mw_dir>>$mapping_res_dir");
	    system("rm $seq_dir.$i.*");
	    print $i;
    }
	    $i=$i+1;
	###对比对结果进行分析

}
close mapping_res;



