#!/bin/bash

get_closest_numbers() {
    number=$1
    lower=$(( (number / 100) * 100 ))
    upper=$(( lower + 100 ))
    echo "$lower $upper"
}

runs=(54966 54967 54968 54969)

path=/sphenix/lustre01/sphnxpro/production/run2auau/physics/ana466_2024p012_v001/DST_TRKR_CLUSTER

NumEvtPerDst=10000
NumEvtPerJob=100
NumJobPerDst=$((${NumEvtPerDst} / ${NumEvtPerJob}))

for ((k=0; k<${#runs[@]}; k++))
do
  out=${runs[$k]}_seg_id.txt
  > ${out}

  result=$(get_closest_numbers ${runs[$k]})
  result_arr=($result)
  numlower=${result_arr[0]}
  numupper=${result_arr[1]}

  nseg=`ls ${path}/run_000${numlower}_000${numupper}/dst/*${runs[$k]}* | wc -l`

  echo run ${runs[$k]} : ${nseg} dsts, $((${nseg} * ${NumJobPerDst})) jobs

  for ((i=0; i<${nseg}; i++))
  do
    for ((j=0; j<${NumJobPerDst}; j++))
    do
      echo "$i $j" >> ${out}
    done
  done

done
