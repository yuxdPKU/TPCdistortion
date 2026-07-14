phiType=( 'central' 'west' 'east' )
types=( 'MI_sim_reco_acts' 'MI_sim_reco_genfit' 'MI_sim_reco_truth' 'MI_sim_reco_acts_truthseeding' 'MI_sim_reco_genfit_truthseeding' )
echo ${#phiType[@]}
for ((k=0; k<${#phiType[@]}; k++))
do
  echo run ${phiType[$k]}
  root -b -q -l draw1D_r_from3D_sim_closure.C"(\"${phiType[$k]}\")"
  root -b -q -l draw1D_z_from3D_sim_closure.C"(\"${phiType[$k]}\")"
  for ((j=0; j<${#types[@]}; j++))
  do
    echo run ${phiType[$k]} ${types[$j]}
    root -b -q -l draw2D_rz_from3D_sim_closure.C"(29,\"${types[$j]}\",\"${phiType[$k]}\")"
  done
done

#selectZ=( 0 10 30 60 )
#echo ${#selectZ[@]}
#for ((k=0; k<${#selectZ[@]}; k++))
#do
#  echo run ${selectZ[$k]}
#  root -b -q -l draw2D_rphi_from3D_sim_closure.C"(29,\"MI_sim_reco_truth_notpot\",${selectZ[$k]})"
#done

#selectR=( 40 50 60 70 )
#echo ${#selectR[@]}
#for ((k=0; k<${#selectR[@]}; k++))
#do
#  echo run ${selectR[$k#]}
#  root -b -q -l draw2D_phiz_from3D_sim_closure.C"(29,\"MI_sim_reco_truth_notpot\",${selectR[$k]})"
#done
