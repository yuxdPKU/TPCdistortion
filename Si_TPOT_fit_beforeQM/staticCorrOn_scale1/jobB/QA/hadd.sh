#runs=( 53285 53630 53756 53877 53534 53744 53876)
#runs=( 53534 53630 53756 53877 53744 53876)
runs=( 53534 )
for ((k=0; k<${#runs[@]}; k++))
do
  echo hadd ${runs[$k]}
  hadd -k -f allqa_${runs[$k]}.root ../../Reconstructed/${runs[$k]}/*qa.root
done
