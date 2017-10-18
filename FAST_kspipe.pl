#!/usr/bin/perl -w
####################################################################
# Calculating Ka, Ks Ka/Ks nuclear and protein alignment identity  #
# from DNA sequences (EST) 					   #
# Original by Kerr Wall 					   #
# Modified by Liying Cui, Severn Everett 			   #
# Modified by Michael R. McKain               			   #
# Date Sept 21 2005 						   #
# Revised Jun 29 2011						   #
# Removed ESTscan May 23 2012					   #
# Rewrite Feb 16 2013--M. R. McKain				#
####################################################################

use strict;
use Getopt::Long;
use Parallel::ForkManager;

my $usage = '--query QueryFastaFile --subj SubjectFastaFile --blast BlastFile --pairs PairsperSET --perid PercentID --alen AlignLength --processes MaxProcesses';

my ($infile, $hitfile, $blastfile);
my $perid = 40;
my $alen = 99;
my $proc=2;
my $setpairs = 500;
GetOptions ("query=s" => \$infile,
	       "subj=s" => \$hitfile,
           "blast=s" => \$blastfile,
	   "pairs=i" => \$setpairs,
           "perid=i" => \$perid,
           "alen=i" => \$alen,
	   "processes=i" => \$proc);

unless ($infile && $hitfile && $blastfile){
   die "$0 $usage\n";
}
my $pm = Parallel::ForkManager->new($proc);
my $qsho = substr($infile, 0, -5);
my $ssho = substr($hitfile, 0, -5);
my $bsho = substr($blastfile, 0, -5);
my $pair_id = 1;
my $count = 1;
if (-e "$qsho\_$ssho\_ks") {
	for (my $i=1; $i<=$count; $i++){
		if (-e "$qsho\_$ssho\_ks.$i"){
			$count++;
		}
		else { 
			system "mv $qsho\_$ssho\_ks $qsho\_$ssho\_ks.$i";
		}
	}
}
system "mkdir $qsho\_$ssho\_ks";
system "mkdir $qsho\_$ssho\_ks/data";


my $prequery = substr($infile, 0, index($infile, "."));
my $prehit = substr($hitfile, 0, index($hitfile, "."));


my %q_index = index_input($prequery, $infile);
my %s_index = index_input($prehit, $hitfile);
get_quality_pairs_ORIG($blastfile,$infile,$hitfile);


my $sets=0;
if($pair_id%$setpairs == 0){
        $sets= $pair_id/$setpairs;
}
else{
        $sets= (int($pair_id/$setpairs)+1);
}
my $j=1;
        for (my $i=1; $i <= $sets; $i++){
                open my $setfile, ">", "$i.set.txt";
                open my $paralogfile, "<", "$qsho\_$ssho\_ks/$qsho\_$ssho.paralogs.txt";
                for (my $k=0; $k<($i*$setpairs-$setpairs); $k++){
                        readline($paralogfile);
                }
                        while($j<=$pair_id && $j <= $i*$setpairs){
                                
				my $file_line = readline($paralogfile);
                                if($file_line){
					chomp $file_line;
                                	print $setfile "$j,$file_line\n";
                                }
				$j++;
                        }
                }
        
my @tasks = (1..$sets);
TASKS:
for my $task (@tasks){

	$pm->start and next TASKS;
	`perl FASTKs/bin/FAST_kspipe_part2.pl $qsho $ssho $task.set.txt $infile $hitfile`;

	$pm->finish;
}
$pm->wait_all_children;

`cat $qsho\_$ssho\_ks/$qsho\_$ssho.kaks.2.txt $qsho\_$ssho\_ks/set*/$qsho\_$ssho.kaks.txt > $qsho\_$ssho\_ks/$qsho\_$ssho.kaks.txt`;
`rm kspipealign* *.set.txt`;
#####################################
sub index_input {
	my %index;
	open my $index_out, ">", "$qsho\_$ssho\_ks/$_[0].index";
	open my $infile, "<", $_[1];
	my $count=0;
	
	while(<$infile>){
		chomp;
		if(/^>/){
			$count++;
			my $seqid = substr($_, 1);
			print $index_out "$_[0]" . "$count\t$seqid\n";
			$index{$seqid} = $_[0].$count;
		}
	}
	return(%index);
} 

#####################################
sub get_quality_pairs_ORIG {

	my %hits;
        my %gene_hits;
        my (%sseqs, %qseqs);
        my $temp_seqid;
        open my $qfile, "<", $_[1];
        while(<$qfile>){
                chomp;
                if(/>/){
			if(/\s/){
				/>(.*?)\s/;
				$temp_seqid = $1;
			}
			else{
                        	$temp_seqid = substr($_, 1);
			}
                }
                else{
                        $qseqs{$temp_seqid}.=$_;
                }
        }
        open my $sfile, "<", $_[2];
        while(<$sfile>){
                chomp;
                if(/>/){
                        if(/\s/){
				/>(.*?)\s/;
				$temp_seqid = $1;
			}
			else{
                        	$temp_seqid = substr($_, 1);
			}
                }
                else{
                        $sseqs{$temp_seqid}.=$_;
                }
        }

	open my $blast_file, "<", $_[0];
	open my $outfile, ">", "$qsho\_$ssho\_ks/$qsho\_$ssho.paralogs.txt";
	my $current_query = ""; 
	while (<$blast_file>) {
		chomp;
		my ($query,$hit,$perc_id,$aln_len, $mismatch, $indel, $qstart,$qend,$hstart,$hend) = split /\s+/;	
		next if ($perc_id >= 100.00);
                next if $query eq $hit;
                #my ($qgene, $sgene);
                #$query =~ /(.+c\d+_g\d+)/;
                #$qgene = $1;
                #$hit =~ /(.+c\d+_g\d+)/;
                #$sgene = $1;

                #next if $qgene eq $sgene;
                my $sublen = length($sseqs{$hit});
                my $qulen = length($qseqs{$query});
                next if (abs(($qend-$qstart)*3)/$qulen) < 0.75 || (abs(($hend-$hstart)*3)/$sublen) < 0.75;
		#next if exists $gene_hits{$sgene}{$qgene} || exists $gene_hits{$qgene}{$sgene};
                #next if exists $hits{$hit}{$query} || exists $hits{$query}{$hit};
                #next if exists $hits{$hit};
                #                #next if exists $hits{$query};
                #                                next if exists $gene_hits{$sgene}{$qgene} || exists $gene_hits{$qgene}{$sgene};
                #                                                next if exists $hits{$hit}{$query} || exists $hits{$query}{$hit};

	#	next if ($perc_id == 100.00);
	#	next if $query eq $hit;
		#next if exists $hits{$hit};
		#next if exists $hits{$query};
		next if exists $hits{$hit}{$query} || exists $hits{$query}{$hit};
#		print "$query,$hit,$perc_id,$aln_len,$qstart\t$qend\t$hstart\t$hend\n";
		if ($perc_id > $perid && $aln_len > $alen){
			print $outfile "$query\t$hit\t$qstart\t$qend\t$hstart\t$hend\n"  if (! exists $hits{$hit}{$query});
			$hits{$hit}{$query} = 1;
			$hits{$query}{$hit} = 1;
			#$gene_hits{$sgene}{$qgene}=1;
                        #$gene_hits{$qgene}{$sgene}=1;
			$pair_id++;
		}
	}
}
#####################################

	
