runs=( 29 )
#types=( 'MI_sim_reco_acts' 'MI_sim_reco_genfit' 'MI_sim_reco_truth' )
types=( 'MI_sim_reco_truth' )
echo ${#runs[@]}
echo ${#types[@]}
for ((k=0; k<${#runs[@]}; k++))
do
  for ((j=0; j<${#types[@]}; j++))
  do
    echo run ${runs[$k]} ${types[$j]}
    root -b -q DistortionCorrectionMatrixInversion.C"(${runs[$k]},\"${types[$j]}\")" 1>log_${runs[$k]}_${types[$j]} 2>&1
  done
done
