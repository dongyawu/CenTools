FILE=./CENH3_ChIPseq
gatk_tmpfile=./CENH3_ChIPseq/gatk_tmp
ref_genome=./0_all_assembly/NH051.fasta

#####  set index
cd $FILE/P051
ln -s $ref_genome NH051.fasta
bowtie2-build NH051.fasta NH051

#####  quality control
fastp -i $FILE/P051/E200013966_L01_76_1.fq.gz -o NH051_1.clean.fastq.gz -I $FILE/P051/E200013966_L01_76_2.fq.gz -O NH051_2.clean.fastq.gz -5 --cut_front_window_size 4 --cut_front_mean_quality 20 -3 --cut_tail_window_size 4 --cut_tail_mean_quality 20 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 -w 10 -R NH051

##### clean reads map to the reference
bowtie2 -p 35 -x $FILE/P051/NH051 -1 NH051_1.clean.fastq.gz  -2  NH051_2.clean.fastq.gz 
              --very-sensitive --no-mixed --no-discordant --maxins 800 --rg-id NH051  
              --rg "PL:ILLUMINA" --rg "SM:NH051" -S NH051.sam

#####  sorted and get bam files
/public/software/lib/miniconda3/bin/samtools view -@ 35 -bS NH051.sam  > NH051.bam
/public/software/lib/miniconda3/bin/samtools sort -@ 35 NH051.bam -o NH051.sorted.bam
/public/software/lib/miniconda3/bin/samtools index -@ 35  NH051.sorted.bam
/public/software/lib/miniconda3/bin/samtools flagstat -@ 35  NH051.sorted.bam > NH051.sorted.flagstat
# rm NH051.sam NH051.bam

#####  Keep only the single best comparison result
samtools view -h ./P051/NH051.sorted.bam > NH051_temp.sam
samtools view NH051_temp.sam | awk '{if ($5 == 255 || $5 > max[$1]) {max[$1] = $5; best_line[$1] = $0}} END {for (key in best_line) {print best_line[key]}}' > NH051.best_alignment.sam
head -n 100 NH051_temp.sam | grep "^@" > head_NH051
cat head_NH051 NH051.best_alignment.sam >> NH051.best_alignment_head.sam
samtools view -bS NH051.best_alignment_head.sam | samtools sort -o NH051.filtered_sorted.bam
mv NH051.filtered_sorted.bam NH051-CK.filtered_sorted.bam
# rm NH051_temp.sam NH051.best_alignment_head.sam head_NH051 NH051.best_alignment.sam


#####  calculate standardized coverage by bamCompare
/public/software/lib/miniconda3/bin/samtools index -@ 5 NH051-CK.filtered_sorted.bam
/public/software/lib/miniconda3/bin/samtools index -@ 5 NH051-B1.filtered_sorted.bam
/public/software/lib/miniconda3/bin/samtools index -@ 5 NH051-B2.filtered_sorted.bam
bamCompare -b1 NH051-B1.filtered_sorted.bam -b2 NH051-CK.filtered_sorted.bam --outFileFormat bedgraph -o NH051-B1_log2ratio_1k.bdg --binSize 1000 --operation log2 -p 5 --extendReads
bamCompare -b1 NH051-B2.filtered_sorted.bam -b2 NH051-CK.filtered_sorted.bam --outFileFormat bedgraph -o NH051-B2_log2ratio_1k.bdg --binSize 1000 --operation log2 -p 5 --extendReads

