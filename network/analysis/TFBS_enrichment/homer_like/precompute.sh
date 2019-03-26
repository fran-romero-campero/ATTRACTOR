#! /bin/bash

promoter_up=( 2000 1500 1000 500 )
promoter_down=( 500 200 0 )
scores=( 100 95 90 85 80 )

for i in `seq 0 3`;
do
   for j in `seq 0 2`;
   do
      for k in `seq 0 4`;
      do
         qsub -N p_${promoter_up[$i]}_${promoter_down[$j]}_${scores[$k]} -o p_${promoter_up[$i]}_${promoter_down[$j]}_${scores[$k]} multiplicity_motifs.sh ${promoter_up[$i]} ${promoter_down[$j]} ${scores[$k]}
      done
   done
done
