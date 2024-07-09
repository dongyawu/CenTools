#!/usr/bin/perl

my $rfasta = shift || die "Usage: perl $0 Ref_fasta Que_fasta Output_Prefix \n";
my $qfasta = shift || die "";
my $prefix = shift || die "";

$threads=40;

system(qq(minimap2 -cx asm5 -t ${threads} $rfasta  $qfasta > ${prefix}.paf));

my $paf="${prefix}.paf";

#Read ref sequences
%refseq;
open IN00,"$rfasta";
while (<IN00>){
	chomp;
	if(/^>/){
		s/>//;
		$id=$_;
	}
	else {
		$refseq{$id}.=$_;
	}
}
close IN00;

#Read query sequences
%qseq;
open IN01,"$qfasta";
while (<IN01>){
	chomp;
	if(/^>/){
		s/>//;
		$id=$_;
	}
	else {
		$qseq{$id}.=$_;
	}
}
close IN01;

#Filter PAF
%quelen;
%reflen;
%allelism;

open OUT11,">$paf.fasta";
open OUT22,">$paf.log";

system(qq(cat $paf | sort -k1,1 -k3,3n > $paf.sort));
open IN1,"$paf.sort";
while (<IN1>){
	chomp;
	($qctg, $qlen, $qstart, $qend, $ori, $rctg, $rlen, $rstart, $rend, $alnbase, $alnlen, $quality,)=split(/\t/,$_);
	$quelen{$qctg}=$qlen;
	$reflen{$rctg}=$rlen;
	if(($alnbase > 20000) and ($quality > 55)){
		$lal=$qctg."-".$rctg;
		if(not exists $rctglab{$lal}){
			if(not exists $allelism{$qctg}){
				$allelism{$qctg}=$rctg;
				$count{$qctg}=1;
			}
			else {
				$allelism{$qctg}.=",".$rctg;
				$count{$qctg}+=1;
			}
			$rctglab{$lal}=1;
			$rmin{$lal}=$rstart;
			$rmax{$lal}=$rend;
			$qmin{$lal}=$qstart;
			$qmax{$lal}=$qend;
		}
		else {
			if($rend>$rmax{$lal}){
				$rmax{$lal}=$rend;
			}
			elsif($rstart<$rmin{$lal}){
				$rmin{$lal}=$rstart;
			}
			if($qend>$qmax{$lal}){
				$qmax{$lal}=$qend;
			}
			elsif($qstart<$qmin{$lal}){
				$qmin{$lal}=$qstart;
			}		
		}
		$qlabel=$qctg."-".$qstart."-".$qend;
		$rlabel=$ori."-".$rctg."-".$rstart."-".$rend;
		$hash{$qlabel}=$rlabel;
		
	}
}
close IN1;

open OUT1,">$paf.filt";
foreach $ctg (sort keys %allelism){
	print OUT1 $ctg."\t".$quelen{$ctg}."\t".$count{$ctg}."\t".$allelism{$ctg}."\n";}
close OUT1;

open IN2,"$paf.sort";
while (<IN2>){
	chomp;
	$line=$_;
	($qctg, $qlen, $qstart, $qend, $ori, $rctg, $rlen, $rstart, $rend, $alnbase, $alnlen, $quality,)=split(/\t/,$line);
	if(($count{$qctg}>1)&&($quality>55)&&($alnbase > 20000)){
		$lal=$qctg."-".$rctg;
		$rminR=sprintf("%.5f",$rmin{$lal}/$rlen);
		$rmaxR=sprintf("%.5f",$rmax{$lal}/$rlen);
		$rRrange=sprintf("%.5f",$rmaxR-$rminR);
		$qminR=sprintf("%.5f",$qmin{$lal}/$qlen);
		$qmaxR=sprintf("%.5f",$qmax{$lal}/$qlen);
		$rleft=$rmin{$lal};
		$rright=$rmax{$lal};
		$qleft=$qmin{$lal};
		$qright=$qmax{$lal};
		if(not exists $e{$lal}){
			$e{$lal}=1;

			###Left Ref-contig mapping against Que-contig
			if (not exists $hash3{$qctg}){
				$hash3{$qctg}=2;
				if(($ori eq "+")&&($rmaxR>0.95)){
					$conf="HC";}
				elsif(($ori eq "-")&&($rminR<0.05)){
					$conf="HC";}
				else {
					$conf="LC";}
				print OUT22 $conf."\t".$qctg."\t".$qlen."\t".$qmin{$lal}."-".$qmax{$lal}."\t".$qminR."-".$qmaxR."\t".$ori."\t".$rctg."\t".$rlen."\t".$rmin{$lal}."-".$rmax{$lal}."\t".$rminR."-".$rmaxR."\t".$rRrange."\n";

				##Sequence merge
				if($conf eq "HC"){
					$id=$ori.$rctg."+".$qctg;
					if (not exists $remove{$rctg}){
						$left=$qleft;
						if($ori eq "+"){
							$seq1=substr($refseq{$rctg}, 0, $rleft);
						}
						else {
							$revc_rctgseq = rev_and_com($refseq{$rctg});
							$seq1=substr($revc_rctgseq, 0, $rlen-$rright);
						}
						$leftctg=$rctg;
					}
					else{
						$seq1="";
						$left=$qright;
					}
					$leftrange=$qleft;
					$rightrange=$qright;
				}
				else{
					$left="";
				}
				print OUT22 "left: ".$left."\n";
			}

			###Right Ref-contig mapping against Que-contig
			elsif($hash3{$qctg} == $count{$qctg}){
				print OUT22 $qright."\t".$rightrange."\n";
				if($qright<$rightrange){
					$conf="LC";
				}
				else{
				if(($ori eq "+")&&($rminR<0.05)){
					$conf="HC";}
				elsif(($ori eq "-")&&($rmaxR>0.95)){
					$conf="HC";}
				else {
					$conf="LC";}
				}
				print OUT22 $conf."\t".$qctg."\t".$qlen."\t".$qmin{$lal}."-".$qmax{$lal}."\t".$qminR."-".$qmaxR."\t".$ori."\t".$rctg."\t".$rlen."\t".$rmin{$lal}."-".$rmax{$lal}."\t".$rminR."-".$rmaxR."\t".$rRrange."\n";

				##Sequence merge
				if ($conf eq "HC"){
					$id.="+".$ori.$rctg;
					if(not exists $remove{$rctg}){
						$right=$qright;
						if($ori eq "+"){
							$seq3=substr($refseq{$rctg}, $rright, $rlen-$rright);
						}
						else{
							$seq3= rev_and_com(substr($refseq{$rctg},0,$rleft));
						}
						$rightctg=$rctg;
					}
					else{
						$seq3="";
						$right=$qleft;
					}
				}
				else {
					$right="";
				}
				print OUT22 "right: ".$right."\n";

				##Output merged sequence
				if(($left ne "") && ($right ne "")){
					print OUT22 $id."\n";
					print OUT22 $leftctg."--".$left."--".$right."--".$rightctg."\n";
					$seq2=substr($qseq{$qctg}, $left, $right-$left-1);
					print OUT11 ">".$id."\n";
					print OUT11 $seq1.$seq2.$seq3."\n";
					print OUT22 "seq1:".length($seq1)."\t"."seq2:".length($seq2)."\t"."seq3:".length($seq3)."\n";
					$remove{$leftctg}=1;
					$remove{$rightctg}=1;
					foreach $aaa (@inrctg){
						$remove{$aaa}=1;
					}
				}
				$id="";
				$left="";
                                $leftctg="";
				$right="";
				$leftctg="";
				$rightctg="";
				$leftrange="";
				$rightrange="";
				@inrctg=();
			}

			##### Ref-contigs Within Query-contig
			else {
				if($qright>$rightrange){
					$rightrange=$qright;
				}
				if ($rRrange>0.8){
					print OUT22 "H#"."\t".$qctg."\t".$qlen."\t".$qmin{$lal}."-".$qmax{$lal}."\t".$qminR."-".$qmaxR."\t".$ori."\t".$rctg."\t".$rlen."\t".$rmin{$lal}."-".$rmax{$lal}."\t".$rminR."-".$rmaxR."\t".$rRrange."\n";
#					$remove{$rctg}=1;
					push @inrctg,$rctg;
					$id.="-".$rctg;
				}
				else {
					print OUT22 "L#"."\t".$qctg."\t".$qlen."\t".$qmin{$lal}."-".$qmax{$lal}."\t".$qminR."-".$qmaxR."\t".$ori."\t".$rctg."\t".$rlen."\t".$rmin{$lal}."-".$rmax{$lal}."\t".$rminR."-".$rmaxR."\t".$rRrange."\n";	
				}
				$hash3{$qctg}+=1;
			}
		}
	}
}


foreach $re (sort keys %refseq){
	if (not exists $remove{$re}){
		print OUT11 ">".$re."\n".$refseq{$re}."\n";
	}
}

close OUT11;
close OUT22;

sub rev_and_com {
	my $s = "";
	my $a = "";
	$s = shift ;
	$a = $s;
	$a =~ tr/atcgATCG/tagcTAGC/;
	return reverse($a);
}

system(qq(rm -f $paf.sort));
system(qq(N50 $paf.fasta $paf.n50 10000));

