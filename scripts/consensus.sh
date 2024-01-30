FAST5_DIR="/import/home2/yhchenmath/Dataset/ONT/LifeScienceData/20230306_BECKI_1/no_sample/20230306_1656_MN41990_FAT91015_0b9a9d99/fast5_pass/"
FASTA="/import/home2/yhchenmath/Code/NanoporeData/config/reference_T.fasta"
FASTQ="$HOME/guppy/guppy_out_2"
LOG="/import/home2/yhchenmath/Log/Temp/"
# Software
nanopolish="/import/home2/yhchenmath/software/nanopolish/nanopolish"
canu="/import/home2/yhchenmath/software/canu/build/bin/canu"
minimap="/import/home2/yhchenmath/software/minimap2/minimap2"
samtools="/import/home2/yhchenmath/software/samtools/bin/samtools"


# $nanopolish index -d $FAST5_DIR $FASTQ/reads_T.fasta
$canu -p test -d $LOG/test_canu genomeSize=100k -nanopore-raw $FASTQ/reads_T.fasta
$canu -correct -p test -d $LOG/test_canu genomeSize=100k -nanopore-raw $FASTQ/reads_T.fasta

