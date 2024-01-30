FAST5_DIR="/import/home2/yhchenmath/Dataset/ONT/LifeScienceData/1_mM/20240116_1747_X2_FAY40021_3dd7153f/fast5_pass/"
CONFIG="rna_r9.4.1_70bps_hac.cfg"
FASTA="/import/home2/yhchenmath/Code/NanoporeData/configs/reference_T.fasta"
FASTQ="$HOME/guppy/guppy_out_1mM"
# Software
nanopolish="/import/home2/yhchenmath/software/nanopolish/nanopolish"
minimap="/import/home2/yhchenmath/software/minimap2/minimap2"
samtools="/import/home2/yhchenmath/software/samtools/bin/samtools"
export HDF5_PLUGIN_PATH=/home/yhchenmath/.local/hdf5/lib/plugin

# Basecall
/import/home2/yhchenmath/guppy/ont-guppy/bin/guppy_basecaller -i $FAST5_DIR -s $FASTQ -c $CONFIG -x 'cuda:2'
# Merged the different guppy fastq outputs
cat $FASTQ/pass/*.fastq > $FASTQ/guppy-merged.fastq
# Converted the merged fastq to fasta format
paste - - - - < $FASTQ/guppy-merged.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $FASTQ/reads.fasta
# Create reads_T.fasta
# awk '/^[^>]/{ gsub("U","T") }1' reads.fasta > reads_T.fasta
# Minimap2 and samtools
$minimap -ax map-ont -uf -t 8 --secondary=no  $FASTA $FASTQ/reads_T.fasta | $samtools sort -o $FASTQ/reads_T-ref.sorted.bam
$samtools index $FASTQ/reads_T-ref.sorted.bam


# Nanopolish
$nanopolish index -d $FAST5_DIR $FASTQ/reads_T.fasta
$nanopolish eventalign --reads $FASTQ/reads_T.fasta \
--bam $FASTQ/reads_T-ref.sorted.bam \
--genome $FASTA \
--signal-index \
--scale-events \
--summary $FASTQ/summary.txt \
--threads 20 > $FASTQ/eventalign.txt


