#!/usr/bin/perl -w
use strict;
my $qsho = $ARGV[0];
my $ssho = $ARGV[1];

my $setid = substr($ARGV[2], 0, index($ARGV[2], "."));
system "mkdir $qsho\_$ssho\_ks/set$setid";
my $prequery = substr($ARGV[3], 0, index($ARGV[3], "."));
my $prehit = substr($ARGV[4], 0, index($ARGV[4], "."));

my %qindex = index_input($prequery);
my %sindex = index_input($prehit);
my %rqindex = index_input_rev($prequery);
my %rsindex = index_input_rev($prehit);

my ($qseqs, $sseqs, $pepq, $peps) = subset_seqs($ARGV[2],%qindex, %sindex);
#my %qseqs = %$qseqs;
#my %sseqs = %$sseqs;
open my $setfile, "<", $ARGV[2];
while(<$setfile>){
        chomp;
        my ($pair_id, $line_info) = split /,/;
        create_fasta($pair_id, $line_info);
        run_muscle($pair_id);
        run_pal2nal($pair_id);
        create_phyml_aln($pair_id);
        run_paml($pair_id);
        get_kaks_values($pair_id, $line_info);
}
        
open my $newfile, ">", "$qsho\_$ssho\_ks/$qsho\_$ssho.kaks.2.txt";
print $newfile "Pair\tKa\tKs\tKa/Ks\tSeqQ\tSeqS\n";        
my $home = `pwd`;
print "$home\n";
        
#####################################
sub index_input {
	my %index;
	open my $infile, "<", "$qsho\_$ssho\_ks/$_[0].index";

	
	while(<$infile>){
		chomp;
		my ($new, $old) = split /\s+/;
		$index{$old}=$new;
	}
	return(%index);
} 

#####################################
sub index_input_rev {
        my %index;
        open my $infile, "<", "$qsho\_$ssho\_ks/$_[0].index";


        while(<$infile>){
                chomp;
                my ($new, $old) = split /\s+/;
                $index{$new}=$old;
        }
        return(%index);
}

#####################################
sub subset_seqs {
	my (%qseqs, %sseqs, %pepq, %peps, %fullq, %fulls, %fullpq, %fullps);
	open my $qin, "<", $ARGV[3];
		my $seqid;
		while(<$qin>){
			chomp;
			if(/^>/){
				$seqid = substr($_, 1);
			}
			else{
				$fullq{$seqid}.=$_;
			}
		}
	open my $sin, "<", $ARGV[4];
		$seqid=();
		while(<$sin>){
			chomp;
			if(/^>/){
				$seqid = substr($_, 1);
			}
			else{
				$fulls{$seqid}.=$_;
			}
		}
	 open my $pqin, "<", $prequery . ".pep";
		$seqid=();
		while(<$pqin>){
			chomp;
			if(/^>/){
				$seqid = substr($_, 1);
			}
			else{
				$fullpq{$seqid}.=$_;
			}
		}
		
	 open my $psin, "<", $prehit . ".pep";
		$seqid=();
		while(<$psin>){
			chomp;
			if(/^>/){
				$seqid = substr($_, 1);
			}
			else{
				$fullps{$seqid}.=$_;
			}
		}
		
	open my $infile, "<", $_[0];
	my $qindex = $_[1];
	my $sindex = $_[2];	
	while(<$infile>){
		chomp;
		my ($pair_id, $line_info) = split /,/;
		my ($qseq, $sseq) = split(/\s+/, $line_info);
		$qseqs{$qindex{$qseq}}=$fullq{$qseq};
		$sseqs{$sindex{$sseq}}=$fulls{$sseq};
		$pepq{$qindex{$qseq}}=$fullpq{$qseq};
		$peps{$sindex{$sseq}}=$fullps{$sseq};
	}
	undef(%fullq);
	undef(%fulls);
	undef(%fullpq);
	undef(%fullps);
	return(\%qseqs, \%sseqs, \%pepq, \%peps);
}
		
		
#####################################
sub create_fasta {
	open my $pepfile, ">", "$qsho\_$ssho\_ks/set$setid\/$_[0].pep";
	my ($pair1, $pair2) = split(/\s+/, $_[1]);
	print $pepfile ">$qindex{$pair1}\n${$pepq}{$qindex{$pair1}}\n>$sindex{$pair2}\n${$peps}{$sindex{$pair2}}\n";
	
	open my $cdnafile, ">", "$qsho\_$ssho\_ks/set$setid\/$_[0].cdna";
	print $cdnafile ">$qindex{$pair1}\n${$qseqs}{$qindex{$pair1}}\n>$sindex{$pair2}\n${$sseqs}{$sindex{$pair2}}\n";
}
#####################################
sub run_muscle {
	system "muscle -in $qsho\_$ssho\_ks/set$setid\/$_[0].pep -out $qsho\_$ssho\_ks/set$setid\/$_[0].pep.muscle.fsa";
	next if (-z "$qsho\_$ssho\_ks/set$setid\/$_[0].pep.muscle.fsa");
}
#####################################
sub run_pal2nal {
	system "perl ~/bin/pal2nal.v14/pal2nal.pl $qsho\_$ssho\_ks/set$setid\/$_[0].pep.muscle.fsa $qsho\_$ssho\_ks/set$setid\/$_[0].cdna -output fasta > $qsho\_$ssho\_ks/set$setid\/$_[0].p2n.fasta";
}
#####################################
sub create_phyml_aln {
	open my $infile, "<", "$qsho\_$ssho\_ks/set$setid\/$_[0].p2n.fasta";
	open my $outfile, ">", "$qsho\_$ssho\_ks/set$setid\/$_[0].p2n.phy";
	my ($len, $seqid, %tempseq);
	while(<$infile>){
		chomp;
		if(/^>/){
			if($seqid){
				$len = length($tempseq{$seqid});
			}
			$seqid = substr($_, 1);
			
		}
		else{
			$tempseq{$seqid} .= $_;
		}
	}
	
	print $outfile "2 $len\n";
	my $nsq;
	for my $seq (sort keys %tempseq){
		print $outfile "$seq  $tempseq{$seq}\n";
	}
}
#####################################
sub run_paml {
	chdir("$qsho\_$ssho\_ks/set$setid/");
	my $pair_id = $_[0];
	# run codeml with null
    my $file_args = qq`
    	seqfile = $_[0].p2n.phy 
    	outfile = $_[0].codeml.out
        `;
    my $codeml_args = `cat ~/bin/codeml.ctl.static`;
    open my $codeml, ">", "$_[0].codeml.ctl" or die "Can't open CODEML";
    print $codeml "$file_args";
    print $codeml "$codeml_args";
    close $codeml;
    system "~/bin/paml4.8/bin/codeml $_[0].codeml.ctl";
	chdir("../../");
}
#####################################
sub get_kaks_values {
        open my $ksout, ">>", "$qsho\_$ssho\_ks/set$setid/$qsho\_$ssho.kaks.txt";
       	open my $codemlfile, "<",  "$qsho\_$ssho\_ks/set$setid\/$_[0].codeml.out";
       	my ($pair1, $pair2) = split(/\s+/, $_[1]);
       	
                while (<$codemlfile>){
                        chomp;
                        my $line = $_;
                        if ($line =~ /^t=/){
                        #       if($line =~ /dN\/dS=\s*(\d+.\d+)\s+/){
                        #               print "$1 This works\n";
                        #               $dnds=$1;
                        #       }       
                                my ($dnds,$dn,$ds) = ($1,$2,$3) if ($line =~ /dN\/dS\s*=\s*(\d+\.\d+)\s+dN\s*=\s*(\d+\.\d+)\s+dS\s*=\s*(\d+\.\d+)/);

                                print $ksout "$_[0]\t$dn\t$ds\t$dnds\t$pair1\t$pair2\n";
                        }
                }

`rm "$qsho\_$ssho\_ks/set$setid\/$_[0].codeml.out" "$qsho\_$ssho\_ks/set$setid\/$_[0].cdna" "$qsho\_$ssho\_ks/set$setid\/$_[0].pep" "$qsho\_$ssho\_ks/set$setid\/$_[0].p2n.phy" "$qsho\_$ssho\_ks/set$setid\/$_[0].codeml.ctl" "$qsho\_$ssho\_ks/set$setid\/$_[0].p2n.fasta" "$qsho\_$ssho\_ks/set$setid\/$_[0].pep.muscle.fsa"`;
}


		
