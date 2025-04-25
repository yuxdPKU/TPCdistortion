#!/bin/bash

get_closest_numbers() {
    number=$1
    lower=$(( (number / 100) * 100 ))
    upper=$(( lower + 100 ))
    echo "$lower $upper"
}

path=/sphenix/u/xyu3/hftg01/DST_FOR_DISTORTION/Reconstructed

#53877 - 400khz
#53876 - 430khz
#53756 - 380khz
#53744 - 300khz
#53630 - 550khz
#53534 - 250khz
#53285 - 70khz
#runs=(53285 53534 53630 53744 53756 53876 53877)
runs=(53534 53630 53744 53756 53876 53877)

nseg=100
for ((k=0; k<${#runs[@]}; k++))
do

  if [ -d "seedfilelist/${runs[$k]}" ]; then
    echo "Dir seedfilelist/${runs[$k]} exists. Delete!"
    /usr/bin/rm -rf seedfilelist/${runs[$k]}
  fi
  mkdir -p seedfilelist/${runs[$k]}

  for ((j=0; j<${nseg}; j++))
  do
    segment=$(printf "%05d" $j)
    out=seedfilelist/${runs[$k]}/seed_${segment}.list
    > ${out}

    ls ${path}/${runs[$k]}/clusters_seeds_${runs[$k]}-${j}-*.root_dst.root > ${out}
  done

  out=${runs[$k]}_seg.txt
  > ${out}

  for ((i=0; i<${nseg}; i++))
  do
    echo "$i" >> ${out}
  done


done

echo "Done"
