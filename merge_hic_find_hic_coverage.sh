#!/bin/bash
#SBATCH -J find_hic_hifi_coverage
#SBATCH -A naiss2025-5-531
#SBATCH -p shared
#SBATCH --mem=200GB
#SBATCH -t 24:00:00
#SBATCH -o %x.out
#SBATCH -e %x.err

# source and actiavte envs
source $HOME/.bash_profile

# define inputs
dirofbams=~/LyrenosaPangenome/NewHiCBams/Diploid
dirofhic=~/LyrenosaPangenome/NewHiC/Diploid
oldhicdir=~/LyrenosaPangenome/NewHiC/Diploid
topuphicdir=/cfs/klemming/projects/supr/yantlab_storage/Denzel/LyrenosaPangenome/HiCTopup/Diploid

libpath=~/Lineages
outputdir=~/LyrenosaPangenome/FindHiCStats/Diploid
nchrs=32

# make outdir
mkdir -p $outputdir
cd $outputdir

# export blast db for hifiadapterfilt
export PATH=$PATH:/cfs/klemming/projects/supr/yantlab_storage/Denzel/Tool_Bin/HiFiAdapterFilt/DB
export NUMEXPR_MAX_THREADS=60

# place information
echo STARTED `date`

# loop through bam files and execute code
for file in $dirofbams/*.bam; do
    echo $file # when looping through a dir you just get the file name like sample1.bam, no absolute oaths
    if [ -f $file ]; then
        
        # store bam file
        BamFile=$(basename "${file}")
        BamPrefix=$(basename "${file}" .bam) # name of the bam/fasta file without the extension

        population="${BamPrefix%%.*}"

	# only run for VEL
	if [[ $population != "VEL3171" ]]; then
		continue
	fi

        # set old hic files
        naive_oldhic1=$oldhicdir/$population/*_1.fq.gz
        naive_oldhic2=$oldhicdir/$population/*_2.fq.gz
	
	# set topup hic files
        naive_newhic1=$topuphicdir/$population/*_1.fq.gz
        naive_newhic2=$topuphicdir/$population/*_2.fq.gz

        echo "original files for old hic1 ${naive_oldhic1},original files for new hic1 ${naive_newhic1}"
        echo "original files for old hic2 ${naive_oldhic2}, original files for new hic2 ${naive_newhic2}"

	# activate env
	conda activate hifiadapterfilt

	# get filtered fastq from bam
	mkdir -p $outputdir
        cd $outputdir
	mkdir -p $BamPrefix # make directory for this specific sample
        cd $BamPrefix
        mkdir -p s1_filter_adapters
        cd s1_filter_adapters
        cp $file ./
        hifiadapterfilt.sh \
                -p $BamPrefix \
                -l 44 \
                -m 97 \
                -t 56 \
                -o .

        # store hifi reads, need to have finished hifiasm pipe, in future maybe let s1 be done before hand
        hifi=$outputdir/$BamPrefix/s1_filter_adapters/"${BamPrefix}.filt.fastq.gz"
	
	# activate seqkit
	conda deactivate
	conda activate seqkit

	# ONLY OLD
	# create and cd into correct dirs
	mkdir -p $outputdir
	cd $outputdir
	mkdir -p $BamPrefix # make directory for this specific sample
        cd $BamPrefix
        mkdir -p s0_query_data
        cd s0_query_data

	echo "Combining hic files incase there is more than 1 in the input data dir"
	cat $naive_oldhic1 > "${population}.merged_1.fq.gz"
        cat $naive_oldhic2 > "${population}.merged_2.fq.gz"
	
	# set new variables for fastqs
	oldhic1=$outputdir/$BamPrefix/s0_query_data/"${population}.merged_1.fq.gz"
	oldhic2=$outputdir/$BamPrefix/s0_query_data/"${population}.merged_2.fq.gz"
        
	# find stats
	echo "Getting old hic stats with seqkit"

        seqkit stats $oldhic1 --threads 30 -Ta -o "${population}.hic1.tsv"
        seqkit stats $oldhic2 --threads 30 -Ta -o "${population}.hic2.tsv"
        seqkit stats $hifi --threads 30 -Ta -o "${population}.hifi.tsv"


	# ONLY TOPUP
	# run stats for new topup only
        cd ../
        mkdir -p s0_query_data_topup
        cd s0_query_data_topup

	echo "Combining hic files incase there is more than 1 in the input data dir"
        cat $naive_newhic1 > "${population}.merged_1.fq.gz"
        cat $naive_newhic2 > "${population}.merged_2.fq.gz"

        # set new variables for fastqs
        newhic1=$outputdir/$BamPrefix/s0_query_data_topup/"${population}.merged_1.fq.gz"
        newhic2=$outputdir/$BamPrefix/s0_query_data_topup/"${population}.merged_2.fq.gz"

	# get stats
        echo "Getting topup hic stats with seqkit"

	seqkit stats $newhic1 --threads 30 -Ta -o "${population}.hic1.tsv"
        seqkit stats $newhic2 --threads 30 -Ta -o "${population}.hic2.tsv"
        seqkit stats $hifi --threads 30 -Ta -o "${population}.hifi.tsv"


	# MERGED TOPUP and OLD
	cd ../
	mkdir -p s0_query_data_merged
	cd s0_query_data_merged

	
	echo "Combining hic files incase there is more than 1 in the input data dir"
        cat $oldhic1 $newhic1 > "${population}.merged_1.fq.gz"
        cat $oldhic2 $newhic2 > "${population}.merged_2.fq.gz"

        # set new variables for fastqs
        concathic1=$outputdir/$BamPrefix/s0_query_data_merged/"${population}.merged_1.fq.gz"
        concathic2=$outputdir/$BamPrefix/s0_query_data_merged/"${population}.merged_2.fq.gz"

        # get stats
        echo "Getting topup hic stats with seqkit"

        seqkit stats $concathic1 --threads 30 -Ta -o "${population}.hic1.tsv"
        seqkit stats $concathic2 --threads 30 -Ta -o "${population}.hic2.tsv"
        seqkit stats $hifi --threads 30 -Ta -o "${population}.hifi.tsv"

	conda deactivate

	# create summary of data
	cd ../
	mkdir -p s0_query_data_summary
	cd s0_query_data_summary

	conda activate pd-env
	
	find_hic_hifi_coverage.py -i $outputdir -od ./
	
	conda deactivate
    fi
done

echo "Done"
