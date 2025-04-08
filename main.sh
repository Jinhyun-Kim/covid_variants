BASE_DIR=/nfs_share/students/jinhyun/covid
REFERENCE=$BASE_DIR/reference/NC_045512.2.fasta
cat $BASE_DIR/wildtype/*.fasta $BASE_DIR/standard/*.fasta $BASE_DIR/alpha/*.fasta $BASE_DIR/beta/*.fasta $BASE_DIR/delta/*.fasta $BASE_DIR/omicron/*.fasta > results/all_samples.fasta

mafft --auto --thread 20 --keeplength --addfragments results/all_samples.fasta $REFERENCE > results/aligned.fasta

python extract_vcf.py 