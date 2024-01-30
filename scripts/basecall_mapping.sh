FAST5_DIR="$HOME/Dataset/ONT/LifeScienceData/2_mM/20240116_1747_X3_FAY41870_3c86da9a/fast5_pass/"
CONFIG="rna_r9.4.1_70bps_hac.cfg"
FASTQ="$HOME/guppy/guppy_out_2mM"
/import/home2/yhchenmath/guppy/ont-guppy/bin/guppy_basecaller -i $FAST5_DIR -s $FASTQ -c $CONFIG -x 'cuda:2'
