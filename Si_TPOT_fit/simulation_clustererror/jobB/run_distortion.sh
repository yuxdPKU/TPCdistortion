runs=( 29 )
types=(
	'MI_sim_reco_FITacts_DISTORTIONINPUTtrue_TRUTHSEEDINGtrue'
	'MI_sim_reco_FITacts_DISTORTIONINPUTtrue_TRUTHSEEDINGfalse'
	'MI_sim_reco_FITacts_DISTORTIONINPUTfalse_TRUTHSEEDINGtrue'
	'MI_sim_reco_FITacts_DISTORTIONINPUTfalse_TRUTHSEEDINGfalse'
	'MI_sim_reco_FITgenfit_DISTORTIONINPUTtrue_TRUTHSEEDINGtrue'
	'MI_sim_reco_FITgenfit_DISTORTIONINPUTtrue_TRUTHSEEDINGfalse'
	'MI_sim_reco_FITgenfit_DISTORTIONINPUTfalse_TRUTHSEEDINGtrue'
	'MI_sim_reco_FITgenfit_DISTORTIONINPUTfalse_TRUTHSEEDINGfalse'
	'MI_sim_reco_FITtruth_DISTORTIONINPUTtrue_TRUTHSEEDINGtrue'
	'MI_sim_reco_FITtruth_DISTORTIONINPUTtrue_TRUTHSEEDINGfalse'
	'MI_sim_reco_FITtruth_DISTORTIONINPUTfalse_TRUTHSEEDINGtrue'
	'MI_sim_reco_FITtruth_DISTORTIONINPUTfalse_TRUTHSEEDINGfalse'
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
