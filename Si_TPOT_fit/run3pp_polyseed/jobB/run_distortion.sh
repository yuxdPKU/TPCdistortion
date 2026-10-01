runs=( 79516 )
types=(
	'MI_data_reco_FITacts_USEMMStrue'
	'MI_data_reco_FITacts_USEMMSfalse'
)
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
