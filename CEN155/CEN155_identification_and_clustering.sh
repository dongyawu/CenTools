#####  1\\  de novo TRASH repeat finding 
sh ./software/TRASH/TRASH_run.sh ${prefix}.fasta --simpleplot --horclass CEN155 --frep 100 --par 5 --o ./${prefix}_repeat


#####  2\\  get consensus for rice accessions
mkdir ./${prefix}_HOR
python 02_consensus.py -i ./${prefix}_repeat/Summary.of.repetitive.regions.${prefix}.fasta.csv -o ./${prefix}_HOR/TRASH_${prefix}_consensus.txt


#####  3\\  TRASH CEN155 confirmation
sh ./software/TRASH/TRASH_run.sh ${prefix}.fasta 
							--seqt /${prefix}_HOR/TRASH_${prefix}_consensus.txt 
							--simpleplot --horclass `less ./${prefix}_HOR/TRASH_${prefix}_consensus.txt | cut -d "," -f 1 | sed -n "2p"`  
							--par 5 
							--o ./${prefix}_HOR


#####  4\\  TRASH confirmation by BLAST
makeblastdb -in ${prefix}.fa -dbtype nucl -out ${prefix}
blastn -query ../CEN155_consensus.fasta -db ${prefix} -evalue 1e-6 -outfmt 6 -num_threads 6 -out ${prefix}_consensus.blast
python 04_consensus_blast_region.py -i ${prefix} -o ${prefix}_BLAST_CR_multi.txt



#####  5\\ selection of representative monomers
## step1 select CENtype monomers ##
python 05_1_select_CENtype.py -i ./all_repeats -o ./allrepeats_CENtype

## step2 statistic in one material  ##
python 05_2_TRASH_monomer_dict_count.py -i1 TRASH_allrepeat_unalign_70m_rmdup.monomerID -i2 ${prefix} -o ${each_chr}_${prefix}_count.txt

## step3 get count table  ##
python 05_3_merge_dict_count.py -i1 material.list -i2 ${each_chr} -o ${each_chr}_70m_counttable.txt

## step4 get pro table  ##
python 05_4_count2pro_repeat.py -i1 ${each_chr}_61m_counttable.txt -i2 material.list -o ${each_chr}_61m_protable.txt

## step5 select_representatives  ##
python 05_5_count2pro_all.py -i material.list -o1 ALL_70m_counttable.txt -o2 ALL_70m_protable.txt
python 05_5_select_representatives.py -i1 ALL_70m_protable.txt 
									  -i2 TRASH_70_allrepeat_unalign_rmdup_line.monomrID 
									  -o1 Selected_70m_representative_protable.txt 
									  -o2 Selected_70m_representative_protable.monomrID

## step6 build phylogenetic tree ##
less Selected_70m_representative_protable.monomrID | while read CID seq length
do
	echo ">${CID}" >> Selected_70m_representative.fasta
	echo "${seq}" >> Selected_70m_representative.fasta
done
muscle -align Selected_70m_representative.fasta -output Selected_70m_representative_align.fa
FastTreeMP -nt Selected_70m_representative_align.fa > Selected_70m_representative_align.tree


#####  6\\ manually split monomers into 15 superfamilies


#####  7\\ remian monomers assgin to the superfamilies
python 07_match_all2group_byEdit.py -i1 TRASH_allrepeat_unalign_70m_rmdup.monomerID -i2 know_6624_group_info.txt -o TRASH_allrepeat_match2groups_all.txt


#####  8\\ annotations for 70 rice accessions
## step1 annotate monomers ##
python 08_1_monomer_anno.py -i1 ${prefix} -o ${prefix}_monomer_anno.txt
for prefix in `less material.list`
do
    awk -v OFS='\t' '{if(($5 == "else") && ($3-$2+1 <140)) {$6 = "group_short"}; print}' ${prefix}_monomer_anno.txt > ${prefix}_monomer_anno_short.txt
    awk -v OFS='\t' '{if(($5 == "else") && ($3-$2+1 > 170)) {$6 = "group_long"}; print}' ${prefix}_monomer_anno_short.txt > ${prefix}_monomer_anno_long.txt
    sed "s/group100/group_rare/g" ${prefix}_monomer_anno_long.txt > ${prefix}_monomer_anno_final.txt
    sed -i '1i Chr\tStart\tEnd\tStrand\tMonomerID\tGroup\tSample\tType' ${prefix}_monomer_anno_final.txt
    rm ${prefix}_monomer_anno_short.txt ${prefix}_monomer_anno_long.txt ${prefix}_monomer_anno.txt
done

## step2 reassign monomer superfamilies according to dected ED in 7 ##
python 08_2_monomer_anno_add.py -i1 ${i} -i2 TRASH_allrepeat_match2groups_all.txt -o ${i}_monomer_anno_add.txt


#####  9\\ add intervals to the annotation file
for each_chr in {01..12}; do python 09_monomer_add_interval.py -i ${prefix} -c Chr${each_chr} -o ${prefix}_full_annotation.rmdup.txt; done
