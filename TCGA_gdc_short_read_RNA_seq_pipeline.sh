python star_align.py \
--genomeDir <star_index_path> \
--FastqFileIn <input_fastq_path> \
--workDir <work_dir> \
--out <output_bam> \
--genomeFastaFiles <reference> \
--runThreadN 8 \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--limitBAMsortRAM 0 \
--alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--twopass1readsN -1 \
--sjdbOverhang 100 \
--outSAMstrandField intronMotif \
--outSAMunmapped Within

### For users without access to the ICGC pipeline:

### Step 1: Building the STAR index.*

STAR
--runMode genomeGenerate
--genomeDir <star_index_path>
--genomeFastaFiles <reference>
--sjdbOverhang 100
--sjdbGTFfile <gencode.v22.annotation.gtf>
--runThreadN 8

### Step 2: Alignment 1st Pass.

STAR
--genomeDir <star_index_path>
--readFilesIn <fastq_left_1>,<fastq_left2>,... <fastq_right_1>,<fastq_right_2>,...
--runThreadN <runThreadN>
--outFilterMultimapScoreRange 1
--outFilterMultimapNmax 20
--outFilterMismatchNmax 10
--alignIntronMax 500000
--alignMatesGapMax 1000000
--sjdbScore 2
--alignSJDBoverhangMin 1
--genomeLoad NoSharedMemory
--readFilesCommand <bzcat|cat|zcat>
--outFilterMatchNminOverLread 0.33
--outFilterScoreMinOverLread 0.33
--sjdbOverhang 100
--outSAMstrandField intronMotif
--outSAMtype None
--outSAMmode None

### Step 3: Intermediate Index Generation.

STAR
--runMode genomeGenerate
--genomeDir <output_path>
--genomeFastaFiles <reference>
--sjdbOverhang 100
--runThreadN <runThreadN>
--sjdbFileChrStartEnd <SJ.out.tab from previous step>

### Step 4: Alignment 2nd Pass.

STAR
--genomeDir <output_path from previous step>
--readFilesIn <fastq_left_1>,<fastq_left2>,... <fastq_right_1>,<fastq_right_2>,...
--runThreadN <runThreadN>
--outFilterMultimapScoreRange 1
--outFilterMultimapNmax 20
--outFilterMismatchNmax 10
--alignIntronMax 500000
--alignMatesGapMax 1000000
--sjdbScore 2
--alignSJDBoverhangMin 1
--genomeLoad NoSharedMemory
--limitBAMsortRAM 0
--readFilesCommand <bzcat|cat|zcat>
--outFilterMatchNminOverLread 0.33
--outFilterScoreMinOverLread 0.33
--sjdbOverhang 100
--outSAMstrandField intronMotif
--outSAMattributes NH HI NM MD AS XS
--outSAMunmapped Within
--outSAMtype BAM SortedByCoordinate
--outSAMheaderHD @HD VN:1.4
--outSAMattrRGline <formatted RG line provided by wrapper>


######################################
######################################
######################################
######################################

# install dependencies
#adapter trimming script: http://www.bcgsc.ca/platform/bioinfo/software/adapter-trimming-for-small-rna-sequencing
#BWA 0.5.7
#samtools 0.1.7
#miRNA profiling 0.2.8




## miRNA pipeline quality metrics
samtools view <abc.bam> | awk '{arr[length($10)]+=1} END {for (i in arr) {print i" "arr[i]}}' | sort -t " " -k1n

#$infile is a miRNA fastq file
#$ref is a genomic reference fasta file
{cluster-host}~> bwa aln $ref $infile > $sai
{cluster-host}~> bwa samse -n 10 $ref $sai $infile > $sam



# project library

{cluster-host}~> mkdir -p <PROJECT>/<LIBID1> 
{cluster-host}~> mkdir <PROJECT>/<LIBID2>
{cluster-host}~> mkdir <PROJECT>/<LIBID3>


# file renaming

{cluster-host}~> cp <MY_LIBID>.bam <PROJECT>/<LIBID1>/<LIBID1>_<INDEX1>.bam
{cluster-host}~> cp <MY_LIBID>_adapter.report <PROJECT>/<LIBID1>/<LIBID1>_<INDEX1>_adapter.report
{cluster-host}~> cp <MY_LIBID>.bam <PROJECT>/<LIBID1>/<LIBID1>_<INDEX2>.bam
{cluster-host}~> cp <MY_LIBID>_adapter.report <PROJECT>/<LIBID1>/<LIBID1>_<INDEX2>_adapter.report
{cluster-host}~> cp <MY_LIBID>.bam <PROJECT>/<LIBID2>/<LIBID2>_<INDEX1>.bam
{cluster-host}~> cp <MY_LIBID>_adapter.report <PROJECT>/<LIBID2>/<LIBID2>_<INDEX1>_adapter.report
{cluster-host}~> cp <MY_LIBID>.bam <PROJECT>/<LIBID1>/<LIBID2>_<INDEX2>.bam
{cluster-host}~> cp <MY_LIBID>_adapter.report <PROJECT>/<LIBID2>/<LIBID2>_<INDEX2>_adapter.report

# bam annotation

{xhost or cluster-host}~> samtools view -h <PROJECT><LIBID1>/<LIBID1>_<INDEX1>.bam > <PROJECT><LIBID1>/<LIBID1>_<INDEX1>.sam


# annotate 

<PROJECT>/<LIBID>/<LIBID>_<INDEX>.sam
<PROJECT>/<LIBID>/<LIBID>_<INDEX>_adapter.report


# annotation.pl

{cluster-host}~> perl <BASEDIR>/v0.2.8/code/annotation/annotate.pl -m <mirbase> -u <ucsc_database> -o <species_code> -p <PROJECT>


# alignment stats

{dbhost or xhost}~> perl <BASEDIR>/v0.2.8/code/library_stats/alignment_stats.pl -p <PROJECT>


# TCGA style results

{dbhost or xhost}~> perl <BASEDIR>/v0.2.8/code/custom\_output/tcga.pl -p <PROJECT>

# expression matricies

{dbhost or xhost}~> perl <BASEDIR>/v0.2.8/code/library\_stats/expression\_matrix.pl -m <mirbase> -o <mirbase\_species\_code> -p <PROJECT>

# library stats / min map

{dbhost}~> perl <BASEDIR>/v0.2.8/code/library\_stats/tcga/expression\_matrix\_mimat.pl -m <mirbase> -o <mirbase\_species\_code> -p <PROJECT>


#visualize results 
 
{xhost}~> perl <BASEDIR>/v0.2.8/code/library_stats/graph_libs.pl -p <PROJECT>









