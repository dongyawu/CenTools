#!/bin/bash

###########################
####Date 2023 July
###########################
####By DW
###########################
####Dependencies: meryl, minimap2, winnowmap2, RagTag, seqkit, NextPolish2, N50...
###########################

echo "WARNING: Edit this Script Before Using!!!"
echo ""
echo "A user-friendly version coming soon..."
echo ""
echo "^-^"
echo ""

if [[ "$#" -lt 3 ]]; then
	echo "Usage: ./gap_filling_polising.sh <query fasta> <subspecies XI/GJ> <prefix> <other assembly e.g. verkko> "
	echo ""
	echo "Note: all related files should be named uniformly!!!"
	echo ""
	exit 0	
fi

threads=10

query=$1
sub=$2
id=$3
prefix=$3
verkko=$4
hifi="/PATH/${id}/*fastq.gz"
echo "HiFi reads: ${hifi}"
ngsdir="/PATH/ngs"
echo "NGS reads: ${ngsdir}/${id}/*clean.f*.gz"

gfapl="/PATH/gap_filling_by_assembly.pl"
dotplot="/PATH/run_unimap_dotplot.sh"

###Synteny Control
###GJ and XI assemblies aligned to corresponding Refs...

refGJ="/PATH/NIP_Shang_MP.fasta"
refXI="/PATH/MH63_Song_MP.fasta"

if [[ $sub == "GJ" ]]; then
	ref=$refGJ
elif [[ $sub == "XI" ]]; then 
	ref=$refXI
else
	ref=$refXI
fi

echo ">>>Scaffold contigs!"
echo "ragtag.py scaffold ${ref} ${query} -f 2000 --remove-small -t ${threads} --aligner minimap2 -o ${id}_ragtag"
ragtag.py scaffold ${ref} ${query} -f 2000 --remove-small -t ${threads} --aligner minimap2 -o ${id}_ragtag
mv ${id}_ragtag/ragtag.scaffold.fasta ${id}_ragtag/${id}_scaffold.fasta
mv ${id}_ragtag/ragtag.scaffold.stats ${id}_ragtag/${id}_scaffold.stats
mv ${id}_ragtag/ragtag.scaffold.agp ${id}_ragtag/${id}_scaffold.agp
mv ${id}_ragtag/ragtag.scaffold.confidence.txt ${id}_ragtag/${id}_scaffold.confidence.txt
mv ${id}_ragtag/ragtag.scaffold.asm.paf ${id}_ragtag/${id}_scaffold.asm.paf
rm -f ${id}_ragtag/ragtag.scaffold.asm.paf.log
rm -f ${id}_ragtag/ragtag.scaffold.err

echo ">>>Keep super-scaffold!" 
cat ${id}_ragtag/${id}_scaffold.fasta | grep ">" | grep "RagTag" | sed "s/>//" > ${id}_ragtag/${id}_spscaf.list
seqkit grep -f ${id}_ragtag/${id}_spscaf.list ${id}_ragtag/${id}_scaffold.fasta > ${id}_ragtag/${id}_chr.fasta
echo ">>>Split super-scaffold into contigs!"
cat ${id}_ragtag/${id}_chr.fasta | seqkit fx2tab | cut -f 2 | sed -r 's/n+/\n/gi'  | cat -n | seqkit tab2fx | seqkit replace -p "(.+)" -r "Contig{nr}" > ${id}_ctg.fasta

###gap-filling-by-assembly
echo ">>>Gap filling by other assembly (e.g. verkko)!"
perl ${gfapl} ${id}_ctg.fasta $verkko ${id}_ctg_gfa
sh ${dotplot} ${id}_ctg_gfa.paf.fasta ${ref} ${id}_ctg_gfa_${sub}
rm -f out.*

###rescaffold into scaffolds
echo ">>>Rescaffold GFA-contigs into scaffolds!"
ragtag.py scaffold ${ref} ${id}_ctg_gfa.paf.fasta -f 2000 --remove-small -t ${threads} --aligner minimap2 -o ${id}_ragtag2
mv ${id}_ragtag2/ragtag.scaffold.fasta ${id}_ragtag2/${id}_gfa_scaffold.fasta
mv ${id}_ragtag2/ragtag.scaffold.stats ${id}_ragtag2/${id}_gfa_scaffold.stats
mv ${id}_ragtag2/ragtag.scaffold.agp ${id}_ragtag2/${id}_gfa_scaffold.agp

seqkit grep -f ${id}_ragtag/${id}_spscaf.list ${id}_ragtag2/${id}_gfa_scaffold.fasta > ${id}_ragtag2/${id}_gfa_chr.fasta
ln -s ${id}_ragtag2/${id}_gfa_chr.fasta ./
cat ${id}_ragtag2/${id}_gfa_scaffold.agp | awk '$5=="U"{print $1"\t"$2"\t"$3}' > ${id}_gfa_chr.gap
N50 ${id}_gfa_chr.fasta ${id}_gfa_chr.n50 10000


########################
###Polish by nextPolish2
########################
echo ">>>Polish using nextPolish2!"
echo ""

assembly="${id}_gfa_chr.fasta"
echo "${assembly} will be polished using HiFi reads!"

###Map hifi reads
echo "  >>>Map HiFi reads using winnowmap......"
meryl count k=15 output ${id}_merylDB ${assembly}
meryl print greater-than distinct=0.9998 ${id}_merylDB > repetitive_k15.txt

winnowmap -t ${threads} -W repetitive_k15.txt -ax map-pb ${assembly} ${hifi} | samtools sort -@ ${threads} -o ${id}_hifi_sort.bam
samtools index -@ ${threads} ${id}_hifi_sort.bam

###Prepare k-mer dataset files
echo "  >>>Prepare k-mer dataset from NGS short reads......"
if [ ! -f "${ngsdir}/${id}/${id}.k21.yak" ]; then
	echo "no prepared Yak library (k=21)"
	echo "Preparing Yak library (k=21)"
	yak count -t ${threads} -k 21 -b 37 -o ${id}.k21.yak <(cat ${ngsdir}/${id}/${id}*.clean.f*gz) <(cat ${ngsdir}/${id}/${id}*.clean.f*gz)
	k21="${id}.k21.yak"
else
	echo "found ${ngsdir}/${id}/${id}.k21.yak"
	k21="${ngsdir}/${id}/${id}.k21.yak"
fi

if [ ! -f "${ngsdir}/${id}/${id}.k31.yak" ]; then
	echo "no prepared Yak library (k=31)"
	echo "Preparing Yak library (k=31)"
	yak count -t ${threads} -k 31 -b 37 -o ${id}.k31.yak <(cat ${ngsdir}/${id}/${id}*.clean.f*gz) <(cat ${ngsdir}/${id}/${id}*.clean.f*gz)
	k31="${id}.k31.yak"
else
	echo "found ${ngsdir}/${id}/${id}.k31.yak"
	k31="${ngsdir}/${id}/${id}.k31.yak"
fi

###Run nextPolish2
echo "  >>>Run nextPolish2...... "
nextPolish2 -t 5 -r ${id}_hifi_sort.bam ${assembly} ${k21} ${k31} > ${id}_gfa_chr_np2.fa

##########
####QV
###########

echo ">>>QV"
###QV using hifi kmers
echo "  >>>QV using HiFi library......"

meryl count k=19 threads=${threads} ${hifi} output ${id}_hifi_reads.meryl
meryl greater-than 5 output ${id}_hifi_reads_gt5.meryl ${id}_hifi_reads.meryl

sh $MERQURY/merqury.sh ${id}_hifi_reads_gt5.meryl ${assembly} ${id}_gfa_chr_np2.fa QV_hifi_out

###QV using ngs kmers
echo "  >>>QV using NGS library......"

if [! -d ${ngsdir}/${id}/${id}_ngs_reads.meryl]; then
	echo "Preparing ngs meryl library"
	meryl count k=19 threads=${threads} ${ngsdir}/${id}/*.clean.f*gz output ${id}_ngs_reads.meryl
	meryl greater-than 5 output ${id}_ngs_reads_gt5.meryl ${id}_ngs_reads.meryl
	ngsmyl="${id}_ngs_reads.meryl"
	ngsmylgt5="${id}_ngs_reads_gt5.meryl"
else 
	echo "found ${ngsdir}/${id}/${id}_ngs_reads.meryl"
	ngsmyl="${ngsdir}/${id}/${id}_ngs_reads.meryl"
	ngsmylgt5="${ngsdir}/${id}/${id}_ngs_reads_gt5.meryl"
fi

sh $MERQURY/merqury.sh ${ngsmylgt5} ${assembly} ${id}_gfa_chr_np2.fa QV_ngs_out

###QV using combined kmers
echo "  >>>QV using HiFi+NGS library......"
meryl union-sum output ${id}_combine_reads.meryl ${id}_hifi_reads.meryl ${ngsmyl}
meryl greater-than 5 output ${id}_combine_reads_gt5.meryl ${id}_combine_reads.meryl

sh $MERQURY/merqury.sh ${id}_combine_reads_gt5.meryl ${assembly} ${id}_gfa_chr_np2.fa QV_combine_out

rm -rf logs

echo ""
echo ">>>Polish FINISHED!"
echo ""
echo "bye!"
