#!/bin/bash

#53877 - 400khz
#53876 - 430khz
#53756 - 380khz
#53744 - 300khz
#53630 - 550khz
#53534 - 250khz
#53285 - 70khz

get_closest_numbers() {
    number=$1
    lower=$(( (number / 100) * 100 ))
    upper=$(( lower + 100 ))
    echo "$lower $upper"
}

runs=(53877 53876 53756 53744 53630 53534 53285)

path=/sphenix/lustre01/sphnxpro/physics/slurp/streaming/physics/ana441_2024p007

NumEvtPerDst=10000
NumEvtPerJob=500
NumJobPerDst=$((${NumEvtPerDst} / ${NumEvtPerJob}))

for ((k=0; k<${#runs[@]}; k++))
do
  out=${runs[$k]}_seg_id.txt
  > ${out}

  result=$(get_closest_numbers ${runs[$k]})
  result_arr=($result)
  numlower=${result_arr[0]}
  numupper=${result_arr[1]}

  nseg=`ls ${path}/run_000${numlower}_000${numupper}/*${runs[$k]}* | wc -l`

  echo run ${runs[$k]} : ${nseg} dsts, $((${nseg} * ${NumJobPerDst})) jobs

  for ((i=0; i<${nseg}; i++))
  do
    for ((j=0; j<${NumJobPerDst}; j++))
    do
      echo "$i $j" >> ${out}
    done
  done

done
