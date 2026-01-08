#!/bin/bash
#SBATCH -J running_yahs_wholepipe_VEL3171
#SBATCH -A naiss2025-5-531
#SBATCH -p memory
#SBATCH --mem=500GB
#SBATCH -t 144:00:00
#SBATCH -o %x.out
#SBATCH -e %x.err

# source and actiavte envs
source $HOME/.bash_profile

# define inputs
#ref=~/LyrenosaPangenome/ReferenceGenomes/GCA_026151155.1/ncbi_dataset/data/GCA_026151155.1/GCA_026151155.1_UiO_Aaren_v1.0_genomic.fna
outputdir=~/LyrenosaPangenome/DiploidYahs_VEL3171
nhap=2
dirofbams=~/LyrenosaPangenome/NewHiCBams/Diploid
dirofhic=~/LyrenosaPangenome/FindHiCStats/Diploid
libpath=~/Lineages

# make outdir
mkdir -p $outputdir
cd $outputdir

# place information
echo STARTED `date`

# export blast db for hifiadapterfilt
export PATH=$PATH:/cfs/klemming/projects/supr/yantlab_storage/Denzel/Tool_Bin/HiFiAdapterFilt/DB
export NUMEXPR_MAX_THREADS=60

# loop through bam files and execute code
for file in $dirofbams/VEL*.bam; do
    echo $file # when looping through a dir you just get the file name like sample1.bam, no absolute oaths
    if [ -f $file ]; then

	cd $outputdir

	#Hifiadapterfilt
	conda activate hifiadapterfilt
        echo "Beginning hifiadpterfilt for ${file} .."

	BamFile=$(basename "${file}")
	BamPrefix=$(basename "${file}" .bam) # name of the bam/fasta file without the extension

	population="${BamPrefix%%.*}"

	# merge old and new hic
	

        # set hic
        hic1=$dirofhic/$BamPrefix/s0_query_data_merged/*_1.fq.gz
        hic2=$dirofhic/$BamPrefix/s0_query_data_merged/*_2.fq.gz

        echo "hic1 ${hic1}"
        echo "hic2 ${hic2}"

	mkdir -p $BamPrefix # make directory for this specific sample
	cd $BamPrefix
	mkdir -p s1_filter_adapters
	cd s1_filter_adapters
	mv $file ./
	hifiadapterfilt.sh \
	        -p $BamPrefix \
        	-l 44 \
        	-m 97 \
        	-t 56 \
        	-o .

	# set variable for filtered reads
	FilteredHifiReads="${BamPrefix}.filt.fastq.gz"
	mv $BamFile $dirofbams # move bam file back to directory of bams
	conda deactivate

	# hifiasm
	conda activate hifiasm
	echo "Beginning hifiasm for ${file} .."
	cd ../
	mkdir -p s2_run_hifiasm
	cd s2_run_hifiasm
	hifiasm -v # confirm hifiasm version
	hifiasm -o $BamPrefix \
        	--n-hap $nhap --h1 $hic1 --h2 $hic2 \
        	-t 96 ../s1_filter_adapters/$FilteredHifiReads
	conda deactivate

	# combine gfa files for hap1 and 2 so we can try to get dual haplome maps
	cat *hic.hap1.p_ctg.gfa *hic.hap2.p_ctg.gfa > "${BamPrefix}.hic.hap_comb.p_ctg.gfa"

	# set hifiasm output var
	HifiasmOutDir=$outputdir/$BamPrefix/s2_run_hifiasm

	cd $outputdir
        cd $BamPrefix
        mkdir -p s3_quality_control
        cd s3_quality_control

	# loop through files and execute code
	# loop through files and execute code
        for gfa in $HifiasmOutDir/*hic.hap*.p_ctg.gfa; do
            echo $gfa # when looping through a dir you just get the file name like sample1.bam, no absolute oaths
            if [ -f $gfa ]; then
		conda activate compleasm
		echo "Found file: ${gfa}"
		GfaFile=$(basename "${gfa}")
                GfaPrefix=$(basename "${gfa}" .p_ctg.gfa)
		# GfaPrefixNoBadChars=${GfaPrefix//./_}
		
		# convert to fasta file
		# awk '$1 == "S" {print ">" $2 "\n" $3}' $gfa > $HifiasmOutDir/"${GfaPrefix}.fa"
		awk '/^S/{print ">"$2;print $3}' $gfa > $HifiasmOutDir/"${GfaPrefix}.fa"
		

		echo "Using prefix ${GfaPrefix}"
		cd $outputdir
		cd $BamPrefix
		cd s3_quality_control
		mkdir -p Compleasm
		cd Compleasm
		mkdir -p $GfaPrefix
		cd $GfaPrefix

		echo "Now for compleasm"
                compleasm run -a $HifiasmOutDir/"${GfaPrefix}.fa" \
                               -o ./ -l brassicales --odb odb12 -L $libpath -t 60
                conda deactivate

                # run gfa stats
                echo "Running gfa stats"
                cd ../../
                mkdir -p Gfastats
                cd Gfastats
                mkdir -p $GfaPrefix
                cd $GfaPrefix

                conda activate gfastats
                gfastats -f $HifiasmOutDir/"${GfaPrefix}.fa" -j 30 > stats.txt
                conda deactivate

                # run qc summary scripts
                echo "Now for qc summary metrics"
                cd ../../
                mkdir -p QC_Summary
                cd QC_Summary
                mkdir -p $GfaPrefix
                cd $GfaPrefix

                conda activate pd-env
                create_qc_summary_table.py --fasta $HifiasmOutDir/"${GfaPrefix}.fa" \
                       --gfastats ../../Gfastats/$GfaPrefix/stats.txt \
                       --compleasm ../../Compleasm/$GfaPrefix/summary.txt \
                       -o output --assem_name $GfaPrefix
                conda deactivate
		
		cd ../../

		# make dir for yahs quast, dont worry about it now
		mkdir -p QuastYahs
                mkdir -p Yahs
                cd Yahs
                mkdir -p $GfaPrefix
                cd $GfaPrefix
		mkdir -p s1_run_yahs_read_qc
		cd s1_run_yahs_read_qc
		mkdir -p tmp_pairtools
		
		# define yahs dir with scaffolds
		yahsdir=$outputdir/$BamPrefix/s3_quality_control/Yahs/$GfaPrefix
		conda activate minimap2_samtools
		
		# index fasta
		samtools faidx $HifiasmOutDir/"${GfaPrefix}.fa"

		# create tsv of name and size
		cut -f1,2 $HifiasmOutDir/"${GfaPrefix}.fa.fai" > $HifiasmOutDir/"${GfaPrefix}.contigsizes"
	
		# index with bwa
		bwa index $HifiasmOutDir/"${GfaPrefix}.fa"

		# map hic reads to fasta from hifiasm
		bwa mem -5SP -T0 -t 96 $HifiasmOutDir/"${GfaPrefix}.fa" $hic1 $hic2 | \
		pairtools parse --min-mapq 40 --walks-policy 5unique \
		--max-inter-align-gap 30 --nproc-in 96 --nproc-out 96 --chroms-path $HifiasmOutDir/"${GfaPrefix}.contigsizes" | \
		pairtools sort --tmpdir=tmp_pairtools | pairtools dedup --mark-dups --output-stats stats.txt | \
		pairtools split --output-pairs mapped.pairs --output-sam -|samtools view -bS -@ 96 | \
		samtools sort -@ 96 -o mapped.PT.bam

		conda deactivate
		conda activate yahs

		yahs $HifiasmOutDir/"${GfaPrefix}.fa" mapped.PT.bam
		
		# copy scaffold to quast dir for later
		cp yahs.out_scaffolds_final.fa $outputdir/$BamPrefix/s3_quality_control/QuastYahs/"${GfaPrefix}.yahs.out_scaffolds_final.fa"

		yahsfasta=$yahsdir/s1_run_yahs_read_qc/yahs.out_scaffolds_final.fa
		
		conda deactivate

		conda activate minimap2_samtools

		cd $yahsdir
		mkdir -p s2_pretext
		cd s2_pretext
		mkdir -p pretextmap
		cd pretextmap

		echo 'bwa pipe'

		bwa index $yahsfasta && \
		bwa mem -t 96 $yahsfasta $hic1 $hic2 | \
		samtools sort -@96 -o hic_sorted.bam && \
		samtools view -@96 -h hic_sorted.bam | PretextMap -o map.pretext --mapq 0

		echo 'hifi coverage'

		# align hifi data to scaffolds
		samtools faidx $yahsfasta

		echo 'minimap mapping hifi reads to scaffold assem'

		minimap2 -ax map-pb -t 48 $yahsfasta $outputdir/$BamPrefix/s1_filter_adapters/$FilteredHifiReads | samtools sort -@48 -O BAM -o hifi_sorted.bam -

		# set var for pretextmap
		premap=$yahsdir/s2_pretext/pretextmap/map.pretext
		hicsortbam=$yahsdir/s2_pretext/pretextmap/hic_sorted.bam
		hifisortbam=$yahsdir/s2_pretext/pretextmap/hifi_sorted.bam
		agp=$yahsdir/s1_run_yahs_read_qc/yahs.out_scaffolds_final.agp


		# tracks
		conda deactivate
		conda activate tracks

		cd ../
		mkdir -p tracks
		cd tracks  

		echo 'tidk telomeres'
		tidk search --fasta $yahsfasta --string TTTAGGG --output telomeres --dir telomeres --extension bedgraph

		cd telomeres

		echo 'telomeres graph'

		cat *.bedgraph | PretextGraph -i $premap -n "telomeres"  
		
		# go back to tracks dir
		cd ../

		mkdir -p coverage
		cd coverage

		echo 'coverage graph hic'

		bedtools genomecov -ibam $hicsortbam -bga > coverage_output_hic.bedgraph
		cat coverage_output_hic.bedgraph | PretextGraph -i $premap -n "coveragehic"

		echo 'coverage graph hifi'

		bedtools genomecov -ibam $hifisortbam -bga > coverage_output_hifi.bedgraph
		cat coverage_output_hifi.bedgraph | PretextGraph -i $premap -n "coveragehifi"
		
		# move back to tracks dir
		cd ../
		mkdir -p gaps
		cd gaps

		echo 'gaps graph'
		grep -w 0$ ../coverage/coverage_output_hic.bedgraph | sed 's/0$/200/g' > gaps.bedgraph 
		cat gaps.bedgraph | PretextGraph -i $premap -n "gaps"
 
		conda deactivate
		
		# rereun qc pipe on yash assembly
		echo "Using prefix ${GfaPrefix}, starting qc for yahs scaffolds"
                cd $yahsdir
		mkdir -p s3_quality_control
                cd s3_quality_control
                mkdir -p Compleasm
                cd Compleasm
		
		conda activate compleasm

                echo "Now for compleasm"
                compleasm run -a $yahsfasta -o ./ -l brassicales -L $libpath --odb odb12 -t 60
                conda deactivate

                # run gfa stats
                echo "Running gfa stats"
                cd ../
                mkdir -p Gfastats
                cd Gfastats

                conda activate gfastats
                gfastats -f $yahsfasta -j 30 > stats.txt
                conda deactivate

                # run qc summary scripts
                echo "Now for qc summary metrics"
                cd ../
                mkdir -p QC_Summary
                cd QC_Summary

                conda activate pd-env
                create_qc_summary_table.py --fasta $yahsfasta \
                        --gfastats ../Gfastats/stats.txt \
                        --compleasm ../Compleasm/summary.txt \
                        -o output --assem_name $GfaPrefix
                conda deactivate	
	

            fi
        done

	# run quast for hifiasm
        echo "Now for running quast for yahs"
        cd $outputdir
        cd $BamPrefix
        mkdir -p s3_quality_control
        cd s3_quality_control

        mkdir -p Quast
        cd Quast

        conda activate quast
        quast -t 60 "${HifiasmOutDir}"/*.fa
        conda deactivate

	# run quast for yahs
        cd $outputdir
        cd $BamPrefix
        mkdir -p s3_quality_control
        cd s3_quality_control

        mkdir -p QuastYahs
        cd QuastYahs

        conda activate quast
        quast -t 60 *.fa
        conda deactivate


    fi
done

echo "Done"
