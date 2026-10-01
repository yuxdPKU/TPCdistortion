runs=( 29 )
#types=( 'MI_sim_reco_acts' 'MI_sim_reco_genfit' 'MI_sim_reco_truth' )
#types=( 'MI_sim_reco_truth_notpot' )
#types=( 'MI_sim_reco_acts_truthseeding' 'MI_sim_reco_genfit_truthseeding' )
#types=( 'MI_sim_reco_truth_extrapolate' )
#types=( 'MI_sim_reco_acts_truthseeding_includesecondaries' 'MI_sim_reco_genfit_truthseeding_includesecondaries' )
types=(
	#'MI_sim_reco_FITacts_EXTRAdefault_CLUSERRraw_DISTORTIONINPUTfalse'
	#'MI_sim_reco_FITacts_EXTRAforward_CLUSERRraw_DISTORTIONINPUTfalse'
	#'MI_sim_reco_FITacts_EXTRAbackward_CLUSERRraw_DISTORTIONINPUTfalse'
	#'MI_sim_reco_FITacts_EXTRAbidirectional_CLUSERRraw_DISTORTIONINPUTfalse'
	#'MI_sim_reco_FITacts_EXTRAbidirectional_CLUSERRsim_DISTORTIONINPUTfalse'
	#'MI_sim_reco_FITacts_EXTRAbidirectional_CLUSERRraw_DISTORTIONINPUTtrue'
	#'MI_sim_reco_FITgenfit_DISTORTIONINPUTtrue'
	#'MI_sim_reco_FITgenfit_DISTORTIONINPUTfalse'
	#'MI_sim_reco_FITtruth_DISTORTIONINPUTtrue'
	#'MI_sim_reco_FITtruth_DISTORTIONINPUTfalse'
	#'MI_sim_reco_FITacts_EXTRAbidirectional_CLUSERRraw_DISTORTIONINPUTfalse_nominalSeeding'
	#'MI_sim_reco_FITacts_EXTRAbidirectional_CLUSERRraw_DISTORTIONINPUTtrue_nominalSeeding'
	'MI_sim_reco_FITgenfit_DISTORTIONINPUTfalse_nominalSeeding'
	'MI_sim_reco_FITgenfit_DISTORTIONINPUTtrue_nominalSeeding'
)
#runs=( 0 )
#types=( 'MI_sim_reco_acts_truthseeding_includesecondaries_constBField' 'MI_sim_reco_genfit_truthseeding_includesecondaries_constBField' 'MI_sim_reco_truth_extrapolate_constBField' )
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
