#!/bin/bash
#SBATCH -J running_pretextmap_post_blobtools_VEL3171_hap2
#SBATCH -A naiss2025-5-531
#SBATCH -p memory
#SBATCH --mem=500GB
#SBATCH -t 36:00:00
#SBATCH -o %x.out
#SBATCH -e %x.err

# after running blobtools on hap1 yahs scaffolds, we will use the fasta with contaminated scaffolds filtered out and rereun pretext map pipe

# source and actiavte envs
source $HOME/.bash_profile

# define inputs
#ref=~/LyrenosaPangenome/ReferenceGenomes/GCA_026151155.1/ncbi_dataset/data/GCA_026151155.1/GCA_026151155.1_UiO_Aaren_v1.0_genomic.fna
outputdir=~/LyrenosaPangenome/DiploidYahs_VEL3171
blob_fasta=~/LyrenosaPangenome/DiploidYahs_VEL3171/VEL3171.hifi_reads/s3_quality_control/Yahs/VEL3171.hifi_reads.hic.hap2/s3_quality_control/Blobtools/viewer/FilteredScaffolds/yahs.out_scaffolds_final.wquery.fa
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

        # set variable for filtered reads
        FilteredHifiReads="${BamPrefix}.filt.fastq.gz"
        conda deactivate

        # hifiasm
        conda activate hifiasm
        echo "Beginning hifiasm for ${file} .."
        cd ../
        mkdir -p s2_run_hifiasm
        cd s2_run_hifiasm
        conda deactivate

        cd $outputdir
        cd $BamPrefix
        mkdir -p s3_quality_control
        cd s3_quality_control

        GfaPrefix=VEL3171.hap2.blobtools.filt.scaffolds_final.wquery
            
        echo "Using prefix ${GfaPrefix}"
            
        mkdir -p Yahs
        cd Yahs
        mkdir -p $GfaPrefix
        cd $GfaPrefix
        mkdir -p s1_run_yahs_read_qc
        cd s1_run_yahs_read_qc
        mkdir -p tmp_pairtools
            
        # define yahs dir with scaffolds
        yahsdir=$outputdir/$BamPrefix/s3_quality_control/Yahs/$GfaPrefix
        
        yahsfasta=$blob_fasta
        
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

echo "Done"
