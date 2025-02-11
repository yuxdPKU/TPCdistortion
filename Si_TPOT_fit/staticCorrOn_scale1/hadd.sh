runs=(53285)

for ((k=0; k<${#runs[@]}; k++))
do
  hadd -f -k -v clusters_seeds_${runs[$k]}_qa_all.root root/clusters_seeds_${runs[$k]}*qa.root
done
