#! /usr/bin/perl -w
use strict;

die "perl $0 <in.fas> <...>\n" unless(@ARGV>=1);

$| = 1;
while( my $fa_file=shift ){
    my $id = "";
    my $len = my $seq_len = 0;
    if($fa_file =~ /\.gz$/){
        open IN,"gzip -dc $fa_file | " or die "$!\n";
    }else{
        open IN,$fa_file or die $!;
    }
    while (<IN>){
        chomp;
        next if(/^\#/);
        s/\r$//;
        if(/^>(\S+)/){
            if($id eq ""){
                $id = $1;
                next;
            }
            print "$id\t$len\t$seq_len\n";
            $id = $1;
            $len = 0;
            $seq_len = 0;
            next;
        }
        $len += length($_);
        s/[^ATCG]//ig;
        $seq_len += length($_);
    }
    print "$id\t$len\t$seq_len\n";
    close IN;
}


#yangxianwei


