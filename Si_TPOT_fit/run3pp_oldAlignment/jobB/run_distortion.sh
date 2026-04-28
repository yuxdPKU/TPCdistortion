runs=( 79507 79516 )
#runs=( 79507 79509 79510 79511 79512 79513 79514 79515 79516 )
echo ${#runs[@]}
for ((k=0; k<${#runs[@]}; k++))
do
  echo run ${runs[$k]}
  root -b -q DistortionCorrectionMatrixInversion.C"(${runs[$k]})" 1>log_${runs[$k]} 2>&1
done
