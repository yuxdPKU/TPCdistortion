#!/bin/bash

get_closest_numbers() {
    number=$1
    lower=$(( (number / 100) * 100 ))
    upper=$(( lower + 100 ))
    echo "$lower $upper"
}

runspecies='run3pp'
runtype='physics'
anabuild='ana532'
cdbtag='2025p009'
version='v001'

path=/sphenix/lustre01/sphnxpro/production/${runspecies}/${runtype}/${anabuild}_${cdbtag}_${version}
echo Searching in ${path}

runs=(79516)

for ((k=0; k<${#runs[@]}; k++))
do
  echo Processing run ${runs[$k]}

  result=$(get_closest_numbers ${runs[$k]})
  result_arr=($result)
  numlower=${result_arr[0]}
  numupper=${result_arr[1]}

  out=lamination_dst_list_${runs[$k]}
  > ${out}

  NumDstPerJob=100
  #nseg=`ls ${path}/DST_TRKR_CLUSTER/run_000${numlower}_000${numupper}/dst/*${runs[$k]}* | wc -l`
  nseg=`ls ${path}/DST_TRKR_CLUSTER/run_000${numlower}_000${numupper}/*${runs[$k]}* | wc -l`
  if [ "$NumDstPerJob" -gt "$nseg" ]; then
    echo Total available number of segments less than ${NumDstPerJob}, set to $nseg
    NumDstPerJob=${nseg}
  fi

  for ((j=0; j<${NumDstPerJob}; j++))
  do
    formatted_segment=$(printf "%05d" "$j")
    #echo "${path}/DST_TRKR_CLUSTER/run_000${numlower}_000${numupper}/dst/DST_TRKR_CLUSTER_${runspecies}_${anabuild}_${cdbtag}_${version}-000${runs[$k]}-${formatted_segment}.root" >> ${out}
    echo "${path}/DST_TRKR_CLUSTER/run_000${numlower}_000${numupper}/DST_TRKR_CLUSTER_${runspecies}_${anabuild}_${cdbtag}_${version}-000${runs[$k]}-${formatted_segment}.root" >> ${out}
  done

done 

echo "Done"
