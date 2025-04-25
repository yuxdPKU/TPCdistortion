#runs=( 53285 53630 53756 53877 53534 53744 53876)
runs=( 53534 53630 53756 53877 53744 53876)
#runs=( 53534 53630 53744 )
for ((k=0; k<${#runs[@]}; k++))
do
  echo run ${runs[$k]}
  root -b -q DistortionCorrectionMatrixInversion.C"(${runs[$k]})" 1>log_${runs[$k]} 2>&1
done
