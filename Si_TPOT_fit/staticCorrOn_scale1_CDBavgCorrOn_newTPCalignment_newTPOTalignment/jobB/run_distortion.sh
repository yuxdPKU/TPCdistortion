#runs=( 53534 53630 53756 53877 53744 53876)
runs=( 53877 )
echo ${#runs[@]}
for ((k=0; k<${#runs[@]}; k++))
do
  echo run ${runs[$k]}
  root -b -q DistortionCorrectionMatrixInversion.C"(${runs[$k]})" 1>log_${runs[$k]} 2>&1
done
