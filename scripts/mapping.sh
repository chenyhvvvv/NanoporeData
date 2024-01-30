FASTQ_DIR="/import/home2/yhchenmath/guppy/guppy_out_2mM/pass/"
# gunzip $FASTQ_DIR/*.gz
OUT_DIR="/import/home2/yhchenmath/guppy/guppy_out_2mM/Analysis"

cat $FASTQ_DIR/*.fastq > $OUT_DIR/guppy-merged.fastq
paste - - - - < $OUT_DIR/guppy-merged.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > $OUT_DIR/reads.fasta
awk '/^[^>]/{ gsub("U","T") }1' $OUT_DIR/reads.fasta > $OUT_DIR/reads_T.fasta


# /import/home2/yhchenmath/Dataset/ONT/LifeScienceData/20230306_BECKI_2/Analysis
# /import/home2/yhchenmath/Dataset/ONT/LifeScienceData/20230306_BECKI_1/Analysis
# /import/home2/yhchenmath/Dataset/ONT/LifeScienceData/2_mM/Analysis
# /import/home2/yhchenmath/Dataset/ONT/LifeScienceData/1_mM/Analysis