#! /usr/bin/perl

my $cps=shift || die "usage: perl $0 cps \n";

my $out="$cps.cps_seq";
my $cps2="$cps.cps";

open IN,"$cps";

$seq0= <IN>;
chomp $seq0;
my ($chr0, $start0, $end0, $strand0, $mcps0)=split(/\t/, $seq0);
my $mgroup0=uc(substr($mcps0,0,1));
my $label0="$chr0-$strand0-$mgroup0";
($m0, $num0)=split(/\*/,$mcps0);
$count0=$num0;

open OUT,">$out";
open OUT2,">$cps2";

print OUT ">$chr0\n";

while (<IN>){
	chomp;
	($chr, $start, $end, $strand, $mcps)=split(/\t/,$_);
	$mgroup=uc(substr($mcps,0,1));
	$label="$chr-$strand-$mgroup";
	($m, $num)=split(/\*/,$mcps);
#	if($chr ne $chr0){
#		print OUT "\n>$chr\n";
#	}
	if($label eq $label0){
		$end0=$end;
		$count0+=$num;
	}
	else{
		print OUT2 "$chr0\t$start0\t$end0\t$strand0\t$mgroup0*$count0\n";
		print OUT "$mgroup0";
		if($chr ne $chr0){
			print OUT "\n>$chr\n";
		}
		$chr0=$chr;
		$start0=$start;
		$end0=$end;
		$strand0=$strand;
		$mgroup0=$mgroup;
		$label0="$chr0-$strand0-$mgroup0";
		$count0=$num;
	}
}
close IN;
close OUT;
close OUT2;
