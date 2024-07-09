#! /usr/bin/perl

my $cps2= shift || die "";

my $label=1;

open OUT0,">$cps2.temp";

open IN, "$cps2";
while (<IN>){
	chomp;
	($chr, $start, $end, $strand, $comp)=split(/\t/,$_);
	if($label % 2 == 1){
		$line1_chr=$chr;
		$line1_start=$start;
		$line1_end=$end;
		$line1_strand=$strand;
		($m1, $num1)=split(/\*/,$comp);
	}
	else {
		$line2_chr=$chr;
		$line2_start=$start;
		$line2_end=$end;
		$line2_strand=$strand;
		($m2, $num2)=split(/\*/,$comp);
		if(($line1_chr eq $line2_chr) and ($line1_strand eq $line2_strand) and ($line2_start-$line1_end < 10000)){
			my @array=($m1,$m2);
			my @sorted = sort { lc($a) cmp lc($b) } @array;
			if ($sorted[0] eq $m1){
				print OUT0 "$chr\t$line1_start\t$line2_end\t$strand\t$m1$m2\t$num1,$num2\t$m1$m2\n";
			}
			else{
				print OUT0 "$chr\t$line1_start\t$line2_end\t$strand\t$m2$m1\t$num2,$num1\t$m1$m2\n";
			}
		}
				
#		if($line1_chr eq $line2_chr){
#			$line1_chr=$line2_chr;
#			$line1_start=$line2_start;
#			$line1_end=$line2_end;
#			$line1_strand=$line2_strand;
#			$label=1;
#		}

	}
	$label+=1;
}

close OUT0;
close IN;


open IN1,"$cps2.temp";
open OUT1,">$cps2.cps";

open OUT21,">$cps2.cps_seq1";
open OUT2,">$cps2.cps_seq";
open OUT3,">$cps2.cps_seq2";

$seq0 = <IN1>;
chomp $seq0;
my ($chr0, $start0, $end0, $strand0, $dyad0, $dyadcom0, $dyadR0)=split(/\t/, $seq0);
my $label0 = "$chr0-$strand0-$dyad0";
my $count0 = 1;
my $note0 = "($dyadcom0)";

print OUT2 ">$chr0\n";
print OUT21 ">$chr0\n";
print OUT3 ">$chr0\n";

while(<IN1>){
	chomp;
	my ($chr, $start, $end, $strand, $dyad, $dyadcom, $dyadR)=split(/\t/, $_);
	$label="$chr-$strand-$dyad";

	if(($label eq $label0) and ($start-$end0 < 20000)){
		$end0=$end;
		$count0+=1;
		$note0.="($dyadcom)";
	}
	else {
		print OUT1 "$chr0\t$start0\t$end0\t$strand0\t$dyad0\t$count0\t$note0\t$dyadR0\n";

		print OUT2 "$dyadR0";
		print OUT21 "$dyad0-";
		print OUT3 "$dyadR0($count0)$strand0";
		if($chr ne $chr0){
			print OUT2 "\n>$chr\n";
			print OUT21 "\n>$chr\n";
			print OUT3 "\n>$chr\n";
		}
                $chr0=$chr;
                $start0=$start;
                $end0=$end;
                $strand0=$strand;
                $dyad0=$dyad;
		$dyadR0=$dyadR;
                $label0="$chr-$strand-$dyad";
                $count0=1;
		$note0 = "($dyadcom)";
	}
}






