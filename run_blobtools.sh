#!/bin/bash
#SBATCH -J run_blobtools
#SBATCH -A naiss2025-5-531
#SBATCH -p shared 
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -t 40:00:00
#SBATCH -o %x.out 
#SBATCH -e %x.err

# Load necessary modules
source $HOME/.bash_profile

### read before running ###
# below are the conda envs used in the pipe, their names, the package installed in them, and the version, and the line their activation is found on in the script
 
# conda environment name : package installed : version : activation line
# seqkit : seqkit : v2.10.0 : line 206
# minimap2_samtools : minimap2 : v2.30-r1287 : line 157
# compleasm : compleasm : v0.2.7 : line 192
# diamond : diamond : v2.1.12 : line 300
# btk : blobtools : v4.4.5 : line 313
# pd env : python : v3.9 : line 353

# with diamond there is a very important caveat, the version of diamond you used to make the uniportdb/taxdump has to be exactly the same as the version in the conda environment
# you will use in this pipe. it doesnt have to be the same as the version above, it just needs to match with what you used for the dbs

### now set all variables
### universal variables (used by different steps) ###

# account to charge cpu hours to, if this changes you will need to change the account at the top of the script, but this variable controls those in the subscripts below
account="naiss2025-5-531"

# this is where you want the whole analysis to be stored
output_dir=~/LyrenosaPangenome/DiploidYahs_VEL3171/VEL3171.hifi_reads/s3_quality_control/Yahs/VEL3171.hifi_reads.hic.hap1/s3_quality_control/Blobtools
assembly=~/LyrenosaPangenome/DiploidYahs_VEL3171/VEL3171.hifi_reads/s3_quality_control/Yahs/VEL3171.hifi_reads.hic.hap1/s1_run_yahs_read_qc/yahs.out_scaffolds_final.fa



### minimap variables ###

# input fastq file, NOTE the minimap step assumes hifi reads are being used
fastq=/cfs/klemming/projects/supr/yantlab_storage/Denzel/LyrenosaPangenome/DiploidYahs_VEL3171/VEL3171.hifi_reads/s1_filter_adapters/VEL3171.hifi_reads.filt.fastq.gz

# amount of memory/RAM assigned to the minimap script, if you run into OOM error, consider changing the partition from shared
# whatever RAM you use, divide by 2 to get a rough estimate of the number of cores you will be getting if the partition is shared
minimap_mem="120GB"

# set minimap partition
minimap_partition="shared"

# number of cpus used for the actual minimap step
minimap_threads=48

# number of time given to the minimap script, this doesnt really matter much because we are only charged for what we used so I set it to one day
minimap_time="24:00:00"

# name for the output bam file
minimap_bam=VEL3171.hap1.vs.reads.bam



### compleasm variables ###

# ram/memory assigned to the compleasm step, divide value by two to know how many cores you are giving the script (if using shared partition)
compleasm_mem="100GB"

# set partition for compleasm script
compleasm_partition="shared"

# number of cpus used for the actual compleasm step
compleasm_threads=32

# number of time given to the compleasm script, this doesnt really matter much because we are only charged for what we used so I set it to one day
compleasm_time="4:00:00"

# lineage to use in compleasm
lineage="brassicales"

# odb used in compleasm
odb="odb12"



### diamond variables ###

# ram/memory assigned to the diamond step, divide value by two to know how many cores you are giving the script (if using shared partition)
diamond_mem="40GB"

# set partition for diamond script
diamond_partition="shared"

# number of cpus used for the actual diamond step
diamond_threads=16

# number of time given to the diamond script, this doesnt really matter much because we are only charged for what we used so I set it to one day
diamond_time="4:00:00"

# scaffolds with a bp length greater than this are considered large scaffolds and are processed accordingly; currently 1Mb
export large_scaffold_length=1000000

# how big the subsampling should be from the large scaffolds, currently 100Kb
subsample_size=100000

# uniprot db used in diamond step
uniprotdb=~/Tool_Bin/blobtools/Databases/UNIPROT/reference_proteomes.dmnd



### blobtools variables ###
taxdumpdir=~/Tool_Bin/blobtools/Databases/taxdump

# this is the name for the final blobdir created
blobdir_outname=VEL3171.bp.hap1


# before running blobtools we need to first:
# 1. run minimap on the assembly
# 2. run compleasm on the assembly
# 3. run diamond on the different scaffolds of the assembly

# create directories
mkdir -p $output_dir
cd $output_dir
mkdir -p Scripts/OnE

# set var for scriptsdir
scriptsdir=$output_dir/Scripts

# step1: run minimap

echo 'Creating minimap script'

# Generate the SLURM job script
job_script="${scriptsdir}/s1_minimap.sh"
cat <<EOF > $job_script
#!/bin/bash
#SBATCH -J s1_minimap
#SBATCH -A ${account}
#SBATCH -p ${minimap_partition}
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --mem=${minimap_mem}
#SBATCH -t ${minimap_time}
#SBATCH -o ${scriptsdir}/OnE/s1_minimap.out
#SBATCH -e ${scriptsdir}/OnE/s1_minimap.err

source \$HOME/.bash_profile

mkdir -p $output_dir
cd $output_dir

mkdir -p Mapping_Depth
cd Mapping_Depth

conda activate minimap2_samtools

minimap2 -ax map-pb -t ${minimap_threads} ${assembly} ${fastq} | samtools sort -@${minimap_threads} -O BAM -o ${minimap_bam} -

conda deactivate

echo 'Done'
EOF

# step2: run compleasm

echo 'Creating Compleasm script'

# Generate the SLURM job script
job_script="${scriptsdir}/s2_compleasm.sh"
cat <<EOF > $job_script
#!/bin/bash
#SBATCH -J s2_compleasm
#SBATCH -A ${account}
#SBATCH -p ${compleasm_partition}
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --mem=${compleasm_mem}
#SBATCH -t ${compleasm_time}
#SBATCH -o ${scriptsdir}/OnE/s2_compleasm.out
#SBATCH -e ${scriptsdir}/OnE/s2_compleasm.err

source \$HOME/.bash_profile

mkdir -p $output_dir
cd $output_dir

mkdir -p Compleasm
cd Compleasm

conda activate compleasm
compleasm run -a ${assembly} -o ./ -l ${lineage} --odb ${odb} -t ${compleasm_threads}
conda deactivate

echo 'Compleasm completed'

echo 'Done'
EOF

# step3: running diamond

echo 'Starting Diamond Pipe'

# activate env
conda activate seqkit

# Ensure output directory exists
mkdir -p $output_dir
cd $output_dir
mkdir -p Diamond
cd Diamond
mkdir -p split_scaffolds
mkdir -p diamond_out

split_scaffold_dir=$output_dir/Diamond/split_scaffolds
diamond_output_dir=$output_dir/Diamond/diamond_out

cd $split_scaffold_dir

# Step 1: Split the assembly into individual scaffolds
seqkit split2 --by-size 1 --out-dir $split_scaffold_dir $assembly

# Step 2: Create dir for large scaffolds
cd $split_scaffold_dir
mkdir -p large_scaffolds

# step 3: find files with large scaffolds and move to large scaffold directory. large scaffold >= xMb
find . -type f -name "*.fa" -exec bash -c '
for file do
    fasta_length=$(grep -v "^>" "$file" | wc -m)
    if (( "$fasta_length" > "$large_scaffold_length" )); then
        mv "$file" large_scaffolds/
    fi
done
' _ {} +

# deactivate env
conda deactivate

# if large scaffold is not empty
if [ -n "$( ls -A large_scaffolds )" ]; then

    echo "large scaffolds found. Beginning fragmenting pipe...."

    # step 4: break large scaffolds into xbp fragments

    # move into large scaffolds dir
    cd large_scaffolds

    # go through all files in large scaffold dir
    for file in *;do

	# echo file name
	echo "Subsampling from ${file}"
	
	# move into dir with large scaffolds
	cd $split_scaffold_dir/large_scaffolds
	
	# get name of scaffold
	scaffold_name=$(basename "$file" .fa)

	# randomly subsample 100kb from the assembly
	randomise_and_subsample_fasta.py -f $file \
			--fragment_size $subsample_size -od ../ \
			-op $scaffold_name


    done
fi


# Step 5: Identify scaffold files and sort them by size (descending)
scaffold_files=($(ls -S $split_scaffold_dir/*.fa))

# Step 6: Loop through each scaffold and create appropriate SLURM jobs
for i in "${!scaffold_files[@]}"; do
    scaffold="${scaffold_files[$i]}"
    scaffold_name=$(basename "$scaffold" .fa)

    # Generate the SLURM job script
    job_script="${scriptsdir}/${scaffold_name}_diamond_job.sh"
    cat <<EOF > $job_script
#!/bin/bash
#SBATCH -J diamond_${scaffold_name}
#SBATCH -A ${account}
#SBATCH -p ${diamond_partition}
#SBATCH -n 1
#SBATCH --ntasks=1
#SBATCH --mem=${diamond_mem} 
#SBATCH -t ${diamond_time}
#SBATCH -o ${scriptsdir}/OnE/${scaffold_name}.out 
#SBATCH -e ${scriptsdir}/OnE/${scaffold_name}.err

source \$HOME/.bash_profile

mkdir -p $diamond_output_dir
cd $diamond_output_dir

conda activate diamond

run_diamond_blobtools.py -f $scaffold --pre $scaffold_name --db $uniprotdb -t ${diamond_threads} --ssize $subsample_size 

echo 'Done'
EOF

done

# run sbatch jobs, and prevent job limit error
run_blobtools.py --sbatch_directory $scriptsdir

# create blobdir
conda activate btk

mkdir -p $output_dir
cd $output_dir

# check if viewer dir exists, if it does then remove it
if [ -d 'viewer' ]; then
    rm -r 'viewer'
fi

# make viewer dir
mkdir -p viewer
cd viewer

# start creating symbolic links and files

# get hits from diamond
cat $output_dir/Diamond/diamond_out/*.out > diamond.hits

# create symbplic link for taxdump
ln -s $taxdumpdir ./

# symbolic link for bam file
ln -s $output_dir/Mapping_Depth/*.bam ./

# symbplic link for assmebly
ln -s $assembly ./

# symbolic link to busco file
ln -s $output_dir/Compleasm/brassicales_odb12/full_table_busco_format.tsv ./

# make adjustments to busco table
sed -i '1s/^/# The lineage dataset is: brassicales_odb12 (Creation date: 2024-01-08, number of genomes: 1, number of BUSCOs: 4311)\n/' full_table_busco_format.tsv
sed -i '1s/^/# Compleasm version is: 5.2.7\n/' full_table_busco_format.tsv

# run blobtools
blobtools create --fasta *.fa --taxdump taxdump --hits *.hits --cov *.bam --busco *.tsv --hits-cols 1=qseqid,2=staxids,3=bitscore,5=sseqid,10=qstart,11=qend,14=evalue --threads 2 $blobdir_outname

conda deactivate

conda activate pd-env

filter_contaminated_scaffolds.py -f $assembly -b $blobdir_outname -od ./FilteredScaffolds

conda deactivate

echo 'Done'

