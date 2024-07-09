#! /usr/bin/perl

my $anno=shift || die "usage: perl $0 anno \n";

my $out="$anno.cps_seq";
my $cps="$anno.cps";

open IN,"$anno";

$seq0=<IN>;

my ($chr0, $start0, $end0, $strand0, $monoID0, $groupID0,)=split(/\t/,$seq0);
my $label0="$chr0-$strand0-$groupID0";

open OUT,">$out";
open OUT2,">$cps";

print OUT ">$chr0\n";

my $count=1;

while (<IN>){
	chomp;
	($chr, $start, $end, $strand, $monoID, $groupID, $accID,)=split(/\t/,$_);
	$label="$chr-$strand-$groupID";
#	if($chr ne $chr0){
#		print OUT "\n>$chr\n";
#	}
	if($label eq $label0){
		$end0=$end;
		$count+=1;
	}
	else{
		print OUT2 "$chr0\t$start0\t$end0\t$strand0\t$groupID0*$count\n";
		print OUT "$groupID0";
		if($chr ne $chr0){
			print OUT "\n>$chr\n";
		}
		$chr0=$chr;
		$start0=$start;
		$end0=$end;
		$strand0=$strand;
		$groupID0=$groupID;
		$label0="$chr0-$strand0-$groupID0";
		$count=1;
	}
}
close IN;
close OUT;
close OUT2;
