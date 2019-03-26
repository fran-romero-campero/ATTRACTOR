#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j yes

Rscript multiplicity_motifs.R $1 $2 $3
