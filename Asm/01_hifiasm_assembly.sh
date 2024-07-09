###define PBS arguments
#PBS -N hifiasm
#PBS -l nodes=1:ppn=40,walltime=30:00:00:00
#PBS -q high
#PBS -j oe
#PBS -o hifiasm_log

###job starting reminder
echo "Starting job at:"
date
hostname

###job dir
cd $PBS_O_WORKDIR

###job command

threads=40

for i in `cat sample.list`;
do

hifi="/PATH/${i}_ccs.fastq.gz"
ont="/PATH/${i}_ont.fastq.gz"

hifiasm --version

if [ -f $ont ]; then
	echo "ONT reads found: ${i}"
	time hifiasm -o ${i}_hifiasm -t ${threads} ${hifi} --ul ${ont}
else 
	echo "ONT reads NOT found: ${i}"
	time hifiasm -o ${i}_hifiasm -t ${threads} ${hifi}
fi

echo "hifiasm finished ${i}"

###NOTE: Primary assembly is used
awk '/^S/{print ">"$2"\n"$3}' ${i}_hifiasm.p_ctg.gfa | fold > ${i}_hifiasm_ctg.fa

#Brief statistic
N50 ${i}_hifiasm_ctg.fa ${i}_hifiasm_ctg.n50 10000

done

echo "Finished!"

date
