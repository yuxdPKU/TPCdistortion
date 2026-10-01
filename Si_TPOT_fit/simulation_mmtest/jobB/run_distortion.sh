runs=( 29 )
#types=( 'MI_sim_reco_acts' 'MI_sim_reco_genfit' 'MI_sim_reco_truth' )
#types=( 'MI_sim_reco_truth_notpot' )
#types=( 'MI_sim_reco_acts_truthseeding' 'MI_sim_reco_genfit_truthseeding' )
#types=( 'MI_sim_reco_truth_extrapolate' )
#types=( 'MI_sim_reco_acts_truthseeding_includesecondaries' 'MI_sim_reco_genfit_truthseeding_includesecondaries' )
runs=( 0 )
types=( 'MI_sim_reco_acts_truthseeding_includesecondaries_chargedgeantino_emptymm')
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
