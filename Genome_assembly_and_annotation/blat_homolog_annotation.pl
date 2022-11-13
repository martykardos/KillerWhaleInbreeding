#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin);
#use Data::Dumper;
#use File::Basename;
sub usage{
    print STDERR <<USAGE;
    ##############################################
	    Updated by Zyllei

      Usage: $0 -p pep.fa -r ref.fa [option]

        Homolog annotation pipeline use blat instead of blast and solar.
      Advanced:
            Step1 blat:
       -c <int>: split number of pep file. Default=50;
            Step2 filter:
       -t <int>: best hits will be conserved in filter PSL. Default=5
       -d <flo>: max divergence rate permitted to best hit. Default=0.2
       -a <flo>: min align coverage permitted in filter PSL. Default=0.3
       -o <flo>: overlap cutoff in nonredundant PSL. Default=0.2
       -g <str>: give a GFF file to reserve the same exon alignment
                in nonredundant step and output statistics file in
                step4. Default=none.
       -k <int>: Skip the gene with intron length N times than ref (need -g), 0 for not filter. Default=0
       -i <int>: Max total intron length of the gene, 0 for not limited. Default=0
            Step3 genewise:
       -l <int>: lines per shell when running genewise. Default=600
       -e <int>: extend length for both sides of region. Default=1000
            Others:
       -m <int>: max job to qsub at the same time. Default=100
       -vf <int>: blat_memory
       -q <str>: other parameter of qsub.Default=" -q st.q -P st_breeding"
       -s <int>: step to run. Default=12345
       -n    : run in local instead of qsub. Default=none

        eg: nohup perl $0 -p homolog.pep -r reference.fa -t 5 -d 0.3 -l 600 -s 12345 &
         or nohup perl $0 -p homolog.pep -r reference.fa -t 5 -d 0.3 -l 600 -s 12345 -q "-q bc.q -P gagtest" &
    ##############################################
USAGE
    exit;
}
&usage if(@ARGV ==0);

############# global variable
my ($pep_file,$ref_file);# global
my ($cut_number);# step1
my ($top_hit,$diff_rate,$align_rate,$overlap_cutoff,$origin_gff_file,$ref_intron_len_restrict,$max_total_intron_len);# step2
my ($genewise_lines,$extend_len);# step3
my ($max_job,$blat_mem,$other_qsub_para,$step,$run_local);# other

GetOptions(
  "p=s"=>\$pep_file,
  "r=s"=>\$ref_file,

  "c=i"=>\$cut_number,

  "t=i"=>\$top_hit,
  "d=f"=>\$diff_rate,
  "a=f"=>\$align_rate,
  "o=f"=>\$overlap_cutoff,
  "g=s"=>\$origin_gff_file,
  "k=i"=>\$ref_intron_len_restrict,
  "i=i"=>\$max_total_intron_len,

  "l=i"=>\$genewise_lines,
  "e=i"=>\$extend_len,

  "m=i"=>\$max_job,
  "vf=i"=>\$blat_mem,
  "q=s"=>\$other_qsub_para,
  "s=i"=>\$step,
  "n"=>\$run_local,
#  ""=>\,
);
############# Check variable
&usage if(!defined $pep_file or !defined $ref_file);
die("unexist $pep_file or $ref_file\n") unless(-e $pep_file and -e $ref_file);

############ set default values
$cut_number||=50;

$top_hit||=5;
$diff_rate=0.2 if(!defined $diff_rate);
$align_rate||=0.3;
$overlap_cutoff||=0.2;
$ref_intron_len_restrict = 0 if(!defined $ref_intron_len_restrict);
$max_total_intron_len = 0 if(!defined $max_total_intron_len);

$genewise_lines||=600;
$extend_len||=1000;## means 5' 1000bp and 3' 1000bp

$max_job||=100;
$other_qsub_para||=' -q st.q -P st_breeding';
$step||=12345;
$run_local||=0;
$blat_mem||=7;
##### Other global variable
my $genewise_query_per_dir=200;
# my $blat_mem=7;

############ file basename
my $pep_basename=`basename $pep_file `;
chomp $pep_basename;
my $ref_file_basename=`basename $ref_file `;
chomp $ref_file_basename;

############ other program
my $blat=-e "$Bin/blat"?"$Bin/blat":'/share/app/blat-319/blat/blat';
my $genewise=-e "$Bin/genewise"?"$Bin/genewise":'/opt/blc/genome/biosoft/wise2.2.0/bin/genewise';
my $sge_qsub=-e "$Bin/qsub_jobs.pl"?"perl $Bin/qsub_jobs.pl":'perl /ifs1/ST_PLANT/USER/liyuxiang/common_bin/blat_homolog_annotation_pipeline/qsub_jobs.pl';
my $gzip=-e "$Bin/gzip"?"$Bin/gzip -9":'/bin/gzip -9';

############ sub directory and shell files
my $pwd=`pwd`;
chomp $pwd;
my $pep_cut_dir="$pwd/$pep_basename.cut";
my $genewise_dir="$pwd/$pep_basename.genewise";

my $blat_result_psl="$pwd/$pep_basename.psl";
my $filter_blat_psl="$pwd/$pep_basename.psl.filter";
my $nonredundant_psl="$pwd/$pep_basename.psl.filter.nonredundant";

my $blat_shell="$pwd/blat_$pep_basename.sh";
my $genewise_shell="$pwd/genewise_$pep_basename.sh";
my $genewise_result_list="$pwd/genewise_result_$pep_basename.lst";
my $genewise_result="$pwd/genewise_result_$pep_basename.genewise";
my $gff_file="$pwd/$pep_basename.genewise.gff";
my $stat_file="$pwd/$pep_basename.genewise.stat.log";

my $clean_dir="$pwd/clean_while_finished";
mkdir($clean_dir) if(!-e $clean_dir);

############################  BEGIN MAIN  ############################

############### change file dir to abs path and link files
if($pep_file!~m{^/}){
    `ln -s $pep_file $pep_basename ` unless(-e $pep_basename);
    $pep_file="$pwd/$pep_file";
}
if($ref_file!~m{^/}){
    `ln -s $ref_file $ref_file_basename ` unless(-e $ref_file_basename);
    $ref_file="$pwd/$ref_file";
}

######## Step 1, cut pep file and run blat
if($step=~/1/){
    print "Start Step 1, Cut pep file and run blat.  ";&GetTime();
    ############### creat directory
    `rm -r $pep_cut_dir ` if(-e $pep_cut_dir);
    mkdir ($pep_cut_dir);

    ############### pep total number statistics
    my $pep_total_number=`grep -c "^>" $pep_file `;
    chomp $pep_total_number;
    my $number_per_file=$pep_total_number % $cut_number==0 ? $pep_total_number/$cut_number : int($pep_total_number/$cut_number)+1;

    ################################# read and cut pep_file
    my @pep_cut_file_prefix;
    my $pep_cut_output_prefix="$pep_cut_dir/pep_cut";
    ############## Check file and parameters and print
    &GetHelp(1);
    print "Step1 para:\nPep_file: $pep_file\nRef_file: $ref_file\nPep_total_number: $pep_total_number\nNumber_in_cut_file: $number_per_file\nRun_local: $run_local\n";
    ############## Begin to cut pep_file
    &CutFasta($pep_file,$number_per_file,$pep_total_number,$pep_cut_output_prefix,\@pep_cut_file_prefix);

    ############################# set blat
    open OUT,">$blat_shell" or die("$!:$blat_shell\n");
    open CLEAN,">$clean_dir/clean_step1.sh" or die("$!:$clean_dir/clean_step1.sh\n");
    ############## Print out shell
    foreach my $e (@pep_cut_file_prefix){
        print OUT "$blat $ref_file $e.fa -q=prot -t=dnax -noHead $e.psl;\n";
        print CLEAN "rm $e.fa $e.psl\n";
    }
    close OUT;
    print CLEAN "rm -r $pep_cut_dir\nrm -r $blat_shell.*\n";
    close CLEAN;

    print "Finished split pep file,begin to run blat.  ";&GetTime();
    ############################# run blat
    if($run_local==0){
        #`$sge_qsub --maxjob $max_job --resource "vf=$step1_mem" --reqsub $blat_shell `;
        `$sge_qsub $blat_shell -m $max_job -v $blat_mem -j 'blat' -o \" $other_qsub_para\" `;
    }else{
        `sh $blat_shell`;
    }
    `for i in $pep_cut_dir/*.psl;do sort -k 10,10 -k 1,1nr \$i;done >$blat_result_psl `;

    print "Finish Step1 at   ";&GetTime();
}

######## Step 2, filter psl, nonredundance
if($step=~/2/){
    print "\nStart Step2, filter PSL.  ";&GetTime();
    ############### Check files and parameter, print
    &GetHelp(2);
    print "Step 2 para:\nBlat_result: $blat_result_psl\nPsl_after_filter: $filter_blat_psl\nPsl_after_nonredundant: $nonredundant_psl\nBest hit: $top_hit\nMin_align_coverage: $align_rate\nMax_divergence_rate: $diff_rate\nOverlap_cutoff_in_nonredundant: $overlap_cutoff\nStat_tmp: $stat_file.tmp\n";
    print "Set ref_gff_file: $origin_gff_file\n" if($origin_gff_file);
    ############################# get top $top_hit result in psl file
    &FilterPSL($blat_result_psl,$filter_blat_psl,$top_hit,$align_rate,$diff_rate);

    ############################# filter redundant
    &NonRedundantPSL($filter_blat_psl,$nonredundant_psl,$overlap_cutoff,$stat_file,$origin_gff_file);
    print "Finish Step2 at   ";&GetTime();
}

######## Step 3, Set and run genewise
if($step=~/3/){
    print "\nStart Step3,Set and run genewise.   ";&GetTime();
    ############# formating and creat directory
    if(-e $genewise_dir){
        print "Deleting directory $genewise_dir\n";
        `rm -r "$genewise_dir" `;
    }
    mkdir($genewise_dir);

    ############# Check file and para, print
    &GetHelp(3);
    print "Step 3 para:\nPep_file: $pep_file\nRef_file: $ref_file\nNonredundant_psl_file: $nonredundant_psl\nGenewise_dir: $genewise_dir\nGenewise_query_per_dir: $genewise_query_per_dir\nExtend_length: $extend_len bp\nRun_local: $run_local\n";
    ############# Read PSL file
    my (@psl_information,%pep_seq,%scaffold_info);
    &ReadPSL($nonredundant_psl,$extend_len,\@psl_information,\%pep_seq,\%scaffold_info);

    ############# Creat genewise config and directory
    &GenewiseConfig($pep_file,$ref_file,$genewise_dir,$genewise_shell,$genewise_query_per_dir,\@psl_information,\%pep_seq,\%scaffold_info);

    print "Finished genewise config, begin to run genewise.   ";&GetTime();
    ################# Run genewise
    if($run_local==0){
        #`$sge_qsub --maxjob $max_job --lines $genewise_lines --reqsub $genewise_shell `;
        `$sge_qsub $genewise_shell -m $max_job -l $genewise_lines -j 'genewise' -o \" $other_qsub_para\" `;
    }else{
        `sh $genewise_shell`;
    }
    print "Finished Step3 at   ";&GetTime();
}

######## Step 4, Convert genewise result
if($step=~/4/){
    print "\nStart Step4, Convert genewise result to GFF3 format.   ";&GetTime();
    ############ Check file and parameter, print
    &GetHelp(4);
    print "Step 4 para:\nGenewise_result_list: $genewise_result_list\nGenewise_result: $genewise_result\nOut_put_gff: $gff_file\nOut_put_stat: $stat_file\n";
    ############ Merge genewise result
    `for i in \` cat $genewise_result_list \`;do cat \$i ;done >$genewise_result `;

    ############ Convert genewise to GFF3 format
    &ConvertFormat($genewise_result,$gff_file,$stat_file);
    print "Finished Step4 at   ";&GetTime();
}
######## Step 5, Check result and clean temp data
if($step=~/5/){
    print "\nStart Step5, Check result and clean temp data.   ";&GetTime();
    ############ Check Step3 results
    my $step3_qsub_dir=(glob("$genewise_shell.*.qsub"))[0];
    my $condition=&CheckWorkFinished($step3_qsub_dir);
    ############ Clean data
    if($condition==0){
        print "Checked, all works were finished, start clean directory.   ";&GetTime();
        `sh "$clean_dir/clean_step1.sh" 2>/dev/null `;
        `sh "$clean_dir/clean_step3_genewise.sh" 2>/dev/null `;
        `rm $genewise_result_list $blat_shell $genewise_shell $stat_file.tmp 2>/dev/null `;
        `$gzip $blat_result_psl ` if(-e $blat_result_psl);
        `$gzip $genewise_result ` if(-e $genewise_result);
        `rm -r $clean_dir `;
        print "Finished Step5 at   ";&GetTime();
    }else{
        print ("\n".('#' x 80)."\n Discover some mistake in $step3_qsub_dir .\n Exit clean data   ");&GetTime();
        print (('#' x 80)."\n");
    }
}
##################### All finished
print "\nRun $0 step: $step\nPipelie finished at   ";&GetTime();
############################   END  MAIN  ############################

###############  sub CutFasta  ###############
####### &CutFasta($fasta_file,$number_per_file,\@empty_array_of_output_file_prefix);
sub CutFasta{
    my ($fasta_file,$num_per_cut,$total_num,$prefix,$cut_file_array)=@_;
    my $number_counting=0;
    my $cut_file_number=1;

    open OUT,">$prefix$cut_file_number.fa" or die("$!:$prefix$cut_file_number.fa\n");
    push @{$cut_file_array},"$prefix$cut_file_number";
    open IN,"<$fasta_file" or die("$!:$fasta_file\n");
    $/='>';<IN>;$/="\n";
    while(my $line=<IN>){
        $line=~/\S+/;
        my $tag=$&;
        $/='>';
        my $seq=<IN>;
        chomp $seq;
        $/="\n";
        $number_counting++;
        print OUT ">$tag\n$seq";
        if($number_counting % $num_per_cut==0 and $number_counting<$total_num){
            $cut_file_number++;
            close OUT;
            open OUT,">$prefix$cut_file_number.fa" or die("$!:$prefix$cut_file_number.fa\n");
            push @{$cut_file_array},"$prefix$cut_file_number";
        }
    }
    close IN;
    close OUT;

    return 0;
}

###############  sub FilterPSL  ###############
####### &FilterPSL($psl_file,$filter_psl_file,$top_hit,$align_rate,$diff_rate);
sub FilterPSL{
    my ($psl,$output,$top_number,$min_coverage,$sub_diff_rate)=@_;
    my $last_query='';
    my $query_len=0;
    my $max_align_len=0;
    my @tmp_content;

    open OUT,">$output" or die("$!:$output\n");
    if($psl=~/\.gz$/){
        open IN,"gzip -dc $psl|" or die("$!:gzip -dc $psl|\n");
    }else{
        open IN,"<$psl" or die("$!:$psl\n");
    }
    while(my $line=<IN>){
        my @info=(split /\s+/,$line);
        ################# Clear data while END of <IN>
        if(eof(IN)){
            if($last_query ne $info[9]){
                print OUT $line;
            }else{
                $max_align_len=$info[0] if($info[0]>$max_align_len);
                push @tmp_content,[($info[0],$line)] if($info[0]/$query_len>=$min_coverage);
                $last_query='';
            }
        }
        ####################### Output data and renew an array while reading a new group
        if($last_query ne $info[9]){
            my $this_number=0;
            foreach my $e (sort { $b->[0]<=>$a->[0] } @tmp_content){
                $this_number++;
                last if($this_number>$top_number or 1-$e->[0]/$max_align_len>$sub_diff_rate);
                print OUT $e->[1];
            }
            ######### clear old data
            @tmp_content=();

            ######### creat new information
            $last_query=$info[9];
            $query_len=$info[10];
            $max_align_len=$info[0];
        }
        $max_align_len=$info[0] if($info[0]>$max_align_len);
        push @tmp_content,[($info[0],$line)] if($info[0]/$query_len>=$min_coverage);
    }
    close IN;
    close OUT;

    return 0;
}

###############  sub NonRedundantPSL  ###############
####### &NonRedundantPSL($psl_file,$nonredundant_psl_output,$overlap_cutoff,$stat_file,$origin_gff_file);
sub NonRedundantPSL{
    my ($psl,$output,$overlap,$stat_tmp,$ref_gff_file)=@_;
    my %stat;
    my %max_len;
    my %print_lines;

    ############## read reference gff file to stat exon number
    my %ref_gff_data;
    my %exon_num;
    if(defined $ref_gff_file){
        open IN,"<$ref_gff_file" or die("$!:$ref_gff_file\n");
        while(my $line=<IN>){
            my @info=(split /\s+/,$line);
            if($info[2] eq 'mRNA'){
                $info[8]=~/^ID=([^;\s]+);?/;
                $ref_gff_data{$1}{'gene_len'} = abs($info[4] - $info[3]) + 1;
                $ref_gff_data{$1}{'intron_len'} = $ref_gff_data{$1}{'gene_len'}; ## mRNA tag must appear before CDS.
            }elsif($info[2] eq 'CDS'){
                $info[8]=~/^Parent=([^;\s]+);?/;
                $ref_gff_data{$1}{'exon_num'}++;
                $ref_gff_data{$1}{'intron_len'} -= abs($info[4] - $info[3]) + 1;
            }
        }
        close IN;
    }
    ############# Read filter psl and nonredundant
    open IN,"<$psl" or die("$!:$psl\n");
    while(my $line=<IN>){
        my @info=(split /\s+/,$line);

        ######## Filter by intron length. Added in v1.4
        next if($max_total_intron_len and $info[7] > $max_total_intron_len);
        next if(defined $ref_gff_file and $ref_intron_len_restrict and $info[7] > $ref_gff_data{$info[9]}{'intron_len'}) * $ref_intron_len_restrict;

        $info[8]=~s/^\S//;
        my $mark=(defined $ref_gff_data{$info[9]}{'exon_num'} and $ref_gff_data{$info[9]}{'exon_num'} - 1 == $info[6]) ? 1 : 0;

        ################# [0]--start   [1]--end   [2]--identity   [3]--query   [4]--right_match   [5]--line number   [6]---max marker
        push @{$stat{$info[13]}{$info[8]}},[($info[15],$info[16],$info[0]/($info[0]+$info[1]),$info[9],$info[0],$.,$mark)];
        $max_len{$info[9]}=$info[0] if(!defined $max_len{$info[9]} or $info[0]>$max_len{$info[9]});
    }
    close IN;

    ############## Foreach queries
    foreach my $e (keys %stat){
        ############ Foreach strands (+ or -)
        foreach my $s (keys %{$stat{$e}}){
            my $address=$stat{$e}{$s};
            my @cluster;
            ################# Sort reagion by start position and end position, then get first cluster
            @{$address}=sort {$a->[0]<=>$b->[0] or $a->[1]<=>$b->[1]} @{$address};
            push @cluster,(shift @{$address});
            $cluster[0][6]=1 if($max_len{$cluster[0][3]}==$cluster[0][4]);
            my $end=$cluster[0][1];
            ################ Loop to find overlapping reagions
            foreach my $j (@{$address}){
                $j->[6]=1 if($max_len{$j->[3]}==$j->[4]);
                if($j->[1]<=$end){
                    $end=$j->[1] if($j->[1]>$end);
                }else{
                    ############### Mark if the match is the max one or only one
                    &FilterCluster(\@cluster,\%print_lines,$overlap);
                    $end=$j->[1];
                    @cluster=();
                }
                push @cluster,$j;
            }
            ################## Filter cluster at last
            &FilterCluster(\@cluster,\%print_lines,$overlap);
        }
    }
    ######################## Read filter psl again and output data
    open OUT,">$output" or die("$!:$output\n");
    open TMP,">$stat_tmp.tmp" or die("$!:$stat_tmp.tmp\n");
    open IN,"<$psl" or die("$!:$psl\n");
    while(my $line=<IN>){
        if(defined $print_lines{$.}){
            my @info=(split /\s+/,$line);
            print OUT $line;
            print TMP "$info[9]\t$info[10]\t$info[13]\t$info[14]";
            print TMP "\t$ref_gff_data{$info[9]}{'exon_num'}" if(defined $ref_gff_data{$info[9]}{'exon_num'}); ## TODO, add intron length and gene length column to stat log
            print TMP "\n";
        }
    }
    close IN;
    close OUT;
}

###############  sub FilterCluster  ###############
####### &FilterCluster(\@cluster,\%print_line,$overlap_cutoff);
sub FilterCluster{
    my ($cp,$pp,$overlap_cf)=@_;
    #################### Sort the cluster,mark=1 > coverage > identity
    @{$cp}=sort { $b->[6]<=>$a->[6] or $b->[4]<=>$a->[4] or $b->[2]<=>$a->[2] } @{$cp};
    ################## First will be the best
    my ($st,$ed)=($cp->[0][0],$cp->[0][1]);
    $pp->{$cp->[0][5]}=0;
    shift @{$cp};
    ################# Return if only one element in cluster
    return 0 unless(@{$cp});
    ################# Make a new cycle to find whose overlap to result < cutoff
    my @new_cluster;
    foreach my $e (@{$cp}){
        ####################### push to new cluster if mark=1 or no overlap
        if($e->[6]==1 or $e->[1]<$st or $e->[0]>$ed){
            push @new_cluster,$e;
        }else{
            ###################### push to new cluster if overlap region / both length < cutoff
            my @middle=sort {$a<=>$b} ($st,$ed,$e->[0],$e->[1]);
            push @new_cluster,$e if(($middle[2]-$middle[1]+1)/($ed-$st+1)<$overlap_cf and ($middle[2]-$middle[1]+1)/($e->[1]-$e->[0]+1)<$overlap_cf);
        }
    }
    ############### Start new loop again if exist elements in new cluster.
    &FilterCluster(\@new_cluster,$pp,$overlap_cf) if(@new_cluster);
    return 0;
}

###############  sub ReadPSL  ###############
####### &ReadPSL($psl,$extent_len,\@psl_info,\%pep_info,\%scaffold_info);
sub ReadPSL{
    my ($psl,$extent,$psl_info,$pep_info,$scaf_info)=@_;
    open IN,"<$psl" or die("$!:$psl\n");
    while(my $line=<IN>){
        my @info=(split /\s+/,$line);

        ############## extent subject seq
        my $subject_start=$info[15]+1-$extent;
        $subject_start=1 if($subject_start<1);
        my $subject_end=$info[16]+$extent;
        $subject_end=$info[14] if($subject_end>$info[14]);

        ############## change strand
        my $plus=$info[8]=~tr/+/+/;
        $plus=($plus==2 or $plus==0)?1:0;

        ############## @{$psl_info} [0]--query  [1]--query_len  [2]--subject  [3]--subject_start  [4]--subject_end  [5]--strand
        push @{$psl_info},[($info[9],$info[10],$info[13],$subject_start,$subject_end,$plus)];

        ############## %{$pep_info}   key: pep_id    value:1,prepare for seq;
        ${$pep_info}{$info[9]}=1;

        ############## %{scaf_info}   key: scaffold_id   value: array of position in @psl_info;
        push @{${$scaf_info}{$info[13]}},$#{$psl_info};
    }
    close IN;

    return 0;
}

###############  sub GenewiseConfig  ###############
####### &GenewiseConfig($pep_file,$ref_file,$genewise_dir,$shell,$genewise_query_per_dir,\@psl_information,\%pep_seq,\%scaffold_info);
sub GenewiseConfig{
    my ($pep,$ref,$dir,$shell,$num_per_dir,$psl_info,$pep_info,$scaf_info)=@_;
    my $all_query_number=@{$psl_info};
    my $counting=0;
    my $sub_dir_number="000";

    ####### Read pep_file
    open IN,"<$pep" or die("$!:$pep\n");
    $/='>';<IN>;$/="\n";
    while(my $line=<IN>){
        $line=~/\S+/;
        my $tag=$&;
        $/='>';
        my $seq=<IN>;
        chomp $seq;
        $/="\n";
        if(defined ${$pep_info}{$tag}){
            $seq=~s/\s//g;
            ############## %{$pep_info}   key: pep_id    value: sequence;
            ${$pep_info}{$tag}=$seq;
        }
    }
    close IN;

    ####### Read scaffold fasta and print genewise shell
    open SH,">$shell" or die("$!:$shell\n");
    open CLEAN,">$clean_dir/clean_step3_genewise.sh" or die("$!:$clean_dir/clean_step3_genewise.sh\n");
    open RESULT,">$genewise_result_list" or die("$!:$genewise_result_list\n");
    open IN,"<$ref" or die("$!:$ref\n");
    $/='>';<IN>;$/="\n";
    REF:while(my $line=<IN>){
        $line=~/\S+/;
        my $tag=$&;
        $/='>';
        my $seq=<IN>;
        chomp $seq;
        $/="\n";
        ################# If this scaffold have homolog alignment.
        if(defined ${$scaf_info}{$tag}){
            $seq=~s/\s//g;
            ################## Foreach alignment, %{scaf_info}   key: scaffold_id   value: array of position in @psl_info;
            foreach my $e (@{${$scaf_info}{$tag}}){
                $counting++;
                ############## Creat sub directory for origin data and result
                if($counting % $num_per_dir==1){
                    $sub_dir_number++;
                    mkdir("$dir/$sub_dir_number");
                }
                my $this=$psl_info->[$e];
                ################# Define filenames of origin data and result
                my $query_fa="$dir/$sub_dir_number/$this->[0].fa";
                my $subject_fa="$dir/$sub_dir_number/$this->[0]_$this->[2]_$this->[3]_$this->[4].fa";
                my $result_file="$dir/$sub_dir_number/$this->[0]_$this->[2]_$this->[3]_$this->[4].genewise";

                ################# Output origin data, query and subject
                open OUT,">$query_fa" or die("$!:$query_fa\n");
                print OUT ">$this->[0]\n${$pep_info}{$this->[0]}\n";
                close OUT;

                open OUT,">$subject_fa" or die("$!:$subject_fa\n");
                print OUT ">$this->[2]_$this->[3]_$this->[4]\n".substr($seq,$this->[3]-1,$this->[4]-$this->[3]+1)."\n";
                close OUT;

                ################# Choose paramester by different strand
                my $choose_strand=($this->[5]==1)?'-tfor':'-trev';
                print SH "$genewise $choose_strand -sum  -genesf  -gff $query_fa $subject_fa > $result_file 2> /dev/null;\n";

                print RESULT "$result_file\n";
                print CLEAN "rm $query_fa $subject_fa $result_file\n";

                ################# Last if all data were output.
                $all_query_number--;
                last REF if($all_query_number==0);
            }
        }
    }
    close IN;
    close SH;
    close RESULT;

    print CLEAN "rm -r $genewise_dir\nrm -r $genewise_shell.*\n";
    close CLEAN;
}

###############  sub ConvertFormat  ###############
####### &ConvertFormat($genewise_result,$gff3_file,$stat_file);
sub ConvertFormat{
    my ($genewise,$output,$gene_stat)=@_;
    my %len;
    my %id_num;

    ################## Read $stat_file.tmp
    my %stat;
    open IN,"<$gene_stat.tmp" or die("$!:$gene_stat.tmp\n");
    while(my $line=<IN>){
        my @info=(split /\s+/,$line);
        $stat{'p'}{$info[0]}=$info[1] if(!defined $stat{'p'}{$info[0]});
        $stat{'s'}{$info[2]}=$info[3] if(!defined $stat{'s'}{$info[2]});
        $stat{'e'}{$info[0]}=$info[4] if(defined $info[4] and !defined $stat{'e'}{$info[0]});
    }
    close IN;
    ################# Convert genewise format
    open OUT,">$output" or die("$!:$output\n");
    open STAT,">$gene_stat" or die("$!:$gene_stat\n");
    print STAT "#gene_ID\tBIT_score\tpep_st\tpep_ed\tref_pep_len\tthis_pep_len\tscaffold\tscaf_len\tstrand\tscaf_st\tscaf_ed\tShift";
    print STAT "\tref_exon_num\tthis_exon_num";
    print STAT "\n";
    if($genewise=~/\.gz$/){
        open IN,"gzip -dc $genewise |" or die ("$!: gzip -dc $genewise |\n");
    }else{
        open IN,"<$genewise" or die("$!:$genewise\n");
    }
    $/="//\n";
    while(my $line1=<IN>){
        next if ($line1!~/^Bits/);
        chomp $line1;
        my $line2=<IN>;
	    chomp $line2;
	    my $line3=<IN>;
	    chomp $line3;

        my ($pep_id,$pep_len,$scaf,$scaf_start,$bit_score,$pst,$ped);
        ######## eg ####### Bits   Query         start end Target      start end   idels introns
        ################### 99.46 scaffold810_unigene_0   19   61 scaffold217_2756383_2757511  629  501    0    0
        ###################  $1        $2                 $3   $4               $5             $6   $7
        ###################                                            $1       $2
        if($line1=~/\n(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/){
            ($bit_score,$pep_id,$pst,$ped)=($1,$2,$3,$4);
            $pep_len=$4-$3+1;
            my $subject=$5;
            if($subject=~/(\S+)_(\d+)_\d+$/){
                $scaf=$1;
                $scaf_start=$2;
            }else{
                die("Uncorrect line:\n$line1\n");
            }
        }else{
            die("Uncorrect file format:\n$line1\n");
        }
        my $cover=sprintf("%.2f",$pep_len/$stat{'p'}{$pep_id}*100);

        my $shift=-1;
	    my @cds;
        #### eg #### scaffold217_2756383_2757511     GeneWise        match   629     501     99.46   -       .       scaffold217_275638...
        ## line3 ### scaffold217_2756383_2757511     GeneWise        cds     629     501     0.00    -       0       scaffold217_275638...
        foreach my $e (split(/\n/,$line3)){
            my @info=(split /\s+/,$e);
            @info[3,4]=@info[4,3] if($info[3]>$info[4]);
            if($info[2] eq 'match'){
                $shift++;
            }elsif($info[2] eq 'cds'){
                $info[3]+=$scaf_start-1;
                $info[4]+=$scaf_start-1;
                push @cds,[@info[3,4,6,7]];
            }
        }
        @cds=sort {$a->[0]<=>$b->[0]} @cds;
        ###################### Change gene name if number>1
        $id_num{$pep_id}++;
        my $gene_id=$pep_id;
        $gene_id.='-D'.$id_num{$pep_id} if($id_num{$pep_id}>1);
        ###################### Output GFF format
        print STAT "$gene_id\t$bit_score\t$pst\t$ped\t$stat{'p'}{$pep_id}\t$pep_len\t$scaf\t$stat{'s'}{$scaf}\t$cds[0][2]\t$cds[0][0]\t$cds[-1][1]\t$shift";
        print STAT "\t$stat{'e'}{$pep_id}\t".(scalar @cds) if(defined $stat{'e'}{$pep_id});
        print STAT "\n";
        print OUT "$scaf\tGeneWise\tmRNA\t$cds[0][0]\t$cds[-1][1]\t$cover\t$cds[0][2]\t.\tID=$gene_id;Shift=$shift;\n";
        map {print OUT join("\t",$scaf,'GeneWise','CDS',$_->[0],$_->[1],'.',$_->[2],$_->[3],"Parent=$gene_id;")."\n" } @cds;
    }
    close IN;
    close OUT;
    close STAT;
    return 0;
}

###############  sub CheckWorkFinished  ###############
####### $finished_condition=&CheckWorkFinished($sge_qsub_dir);
sub CheckWorkFinished{
    my $qsub_dir=shift;
    my $finish_word='The-Work-is-completed!';
    foreach my $e (glob("$qsub_dir/*.sh")){
        my $finish_condition=0;
        ##################### Statistics total jobs
        my $total_jobs=0;
        $total_jobs=`grep -cw $finish_word "$e" `;
        chomp $total_jobs;
        ##################### Statistics echo *.o*
        foreach my $t (glob("$e.o*")){
            my $job_counting=0;
            $job_counting=`grep -cw $finish_word $t `;
            chomp $job_counting;
            if($job_counting==$total_jobs){
                my $err_log=$t;
                $err_log=~s/\.sh\.o/\.sh\.e/;
                if(-e $err_log and (stat $err_log)[7]==0){
                    $finish_condition=1;
                    last;
                }
            }
        }
        return 1 if($finish_condition==0);
    }
    return 0;
}

###############  sub GetTime  ###############
####### &GetTime();
sub GetTime{
    my $this_time=localtime;
    print $this_time."\n";
    return 0;
}

###############  sub GetHelp  ###############
####### &GetHelp($help_step);
sub GetHelp{
    my $step_num=shift;
    my $wrong;
    ####################### Check step 1
    if($step_num==1){
        $wrong.="Not exist pep_file: \"$pep_file\", please check.\n" unless(-e $pep_file);
        $wrong.="Not exist ref_file: \"$ref_file\", please check.\n" unless(-e $ref_file);
        $wrong.="Can't creat pep_cut_dir: \"$pep_cut_dir\", please check.\n" unless(-d $pep_cut_dir);
    }
    ###################### Check step 2
    if($step_num==2){
        $wrong.="Not exist blat_result_file:\"$blat_result_psl\", please check\n" unless(-e $blat_result_psl);
        $wrong.="Not correct int number of best hit para:\"$top_hit\", please set an int number\n" if($top_hit=~/\D/ or $top_hit-int($top_hit)!=0);
        $wrong.="Not correct float number of align coverage para:\"$align_rate\", please set a float number.\n" if($align_rate>1 or $align_rate<0);
        $wrong.="Not correct float number of divergence rate para:\"$diff_rate\", please set a float number.\n" if($diff_rate>1 or $diff_rate<0);
        $wrong.="Not correct float number of overlap cutoff para:\"$overlap_cutoff\", please set a float number.\n" if($overlap_cutoff>1 or $overlap_cutoff<0);
    }
    if($step_num==3){
        $wrong.="Not exist blat_nonredundant_file:\"$nonredundant_psl\", please check\n" unless(-e $nonredundant_psl);
        $wrong.="Not exist pep_file: \"$pep_file\", please check.\n" unless(-e $pep_file);
        $wrong.="Not exist ref_file: \"$ref_file\", please check.\n" unless(-e $ref_file);
        $wrong.="Can't creat pep_cut_dir: \"$genewise_dir\", please check.\n" unless(-d $genewise_dir);
        $wrong.="Not correct int number of extend length para:\"$extend_len\", please set an int number\n" if($extend_len=~/\D/ or $extend_len-int($extend_len)!=0);
    }
    if($step_num==4){
        $wrong.="Not exist stat_tmp_file: \"$stat_file.tmp\", please check.\n" unless(-e "$stat_file.tmp");
        $wrong.="Not exist genewise_result_list: \"$genewise_result_list\", please check.\n" unless(-e $genewise_result_list);
    }
    if(!defined $wrong){
        return 0;
    }else{
        print STDERR $wrong;
        exit;
    }
}
