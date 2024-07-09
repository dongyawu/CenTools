###define PBS arguments
#PBS -N verkko
#PBS -l nodes=1:ppn=40,walltime=30:00:00:00
#PBS -q high
#PBS -j oe
#PBS -o vk_log

###job starting reminder
#echo "Starting job at:"
date
hostname

###job dir
cd $PBS_O_WORKDIR

###job command

threads=30

for i in `sample.list`;

do

hifi="/PATH/${i}_ccs.fastq.gz"
ont="/PATH/${i}_ont.fastq.gz"

if [ -f ${ont} ]; then
	nano="--nano ${ont}"
else 
	nano=""
fi

time verkko -d ${i}_verkko \
--hifi $hifi $nano \
--local-memory 300 \
--sto-run ${threads} 300 24 \
--mer-run ${threads} 300 24 \
--ovb-run ${threads} 300 24 \
--red-run ${threads} 300 24 \
--mbg-run ${threads} 300 24 \
--utg-run ${threads} 300 24 \
--spl-run ${threads} 300 24 \
--ali-run ${threads} 300 24 \
--pop-run ${threads} 300 24 \
-utp-run ${threads} 300 24 \
-lay-run ${threads} 300 24 \
-sub-run ${threads} 300 24 \
-par-run ${threads} 300 24 \
-cns-run ${threads} 300 24 \
--threads ${threads}

cat ${i_verkko}/assembly.fasta | seqkit fx2tab | cut -f 2 | sed -r 's/n+/\n/gi' | cat -n | seqkit tab2fx | seqkit replace -p "(.+)" -r "Contig{nr}" > ${id}_verkko_ctg.fasta

done

echo "Finished!"

date
