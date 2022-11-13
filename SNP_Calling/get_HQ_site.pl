use strict;
use warnings;

my ($vcf, $out) = @ARGV;

die "perl $0 <in.vcf> <out.HQ.vcf>\n indivdual range 3-10 (at 0.1 speed)\n" unless (@ARGV == 2);

my $cutoff = 0.1;
my $max = 10;
my $min = 3;
my $select;

if ($vcf =~ /\.gz$/) {
    open FL, "gzip -dc $vcf|";
} else {
    open FL, $vcf;
};
open FLS, ">$out";

while (<FL>) {
    chomp;
    my @tmp = split;
    if (/^#/) {
        if (/^#CHROM/) {
            my $number = $#tmp - 8;
            $select = int($number*$cutoff);
            $select = $min if ($select < $min);
            $select = $max if ($select > $max);
            $select = $number if ($number < 3); # total individual number
            print FLS "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        } else {
            print FLS "$_\n";
        };
        next;
    };
    
    next if ($tmp[4] =~ /,/); # multi-allele
    my $indv_num = 0;
    foreach my $p (9 .. $#tmp) {
        my $gt = (split /:/, $tmp[$p])[0];
        $indv_num++ if ($gt eq '0/1' or $gt eq '1/1');
    };
    
    # output
    if ($indv_num >= $select) {
        print FLS "$tmp[0]\t$tmp[1]\t.\t$tmp[3]\t$tmp[4]\t.\tPASS\t.\n";
    };
};
close FLS;
close FL;
