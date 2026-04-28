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
cdbtag='nocdbtag'
version='v001'

path=/sphenix/lustre01/sphnxpro/production/${runspecies}/${runtype}/${anabuild}_${cdbtag}_${version}
echo Searching in ${path}

runs=(79516)

NumEvtPerDst=1000
NumEvtPerJob=1000
NumJobPerDst=$((${NumEvtPerDst} / ${NumEvtPerJob}))

for ((k=0; k<${#runs[@]}; k++))
do

  run=${runs[$k]}
  if [ -d "filelist/${runs[$k]}" ]; then
    echo "Dir filelist/${runs[$k]} exists. Delete!"
    /usr/bin/rm -rf filelist/${runs[$k]}
  fi
  mkdir -p filelist/${runs[$k]}

  result=$(get_closest_numbers ${runs[$k]})
  result_arr=($result)
  numlower=${result_arr[0]}
  numupper=${result_arr[1]}

  nseg=100000000

  if [[ "${runspecies}" == "run2pp" ]]; then
    subsystem='TPC'
    for ((i=0; i<24; i++))
    do
      ebdc=$(printf "%02d" $i)
      nsegtemp=`ls ${path}/DST_STREAMING_EVENT_${subsystem}${ebdc}/run_000${numlower}_000${numupper}/dst/*${run}* | wc -l`
      echo run ${run} ${subsystem}${ebdc} raw hit DST ${nsegtemp} segments
      if [ $nseg -gt $nsegtemp ]; then
        nseg=$nsegtemp
      fi
    done
  fi

  if [[ "${runspecies}" == "run3auau" || "${runspecies}" == "run3pp" ]]; then
    subsystem='ebdc'
    for ((i=0; i<24; i++))
    do
      ebdc=$(printf "%02d" $i)
      for ((ii=0; ii<2; ii++))
      do
        nsegtemp=`ls ${path}/DST_STREAMING_EVENT_${subsystem}${ebdc}_${ii}/run_000${numlower}_000${numupper}/*${run}* | wc -l`
        echo run ${run} ${subsystem}${ebdc}_${ii} raw hit DST ${nsegtemp} segments
        if [ $nseg -gt $nsegtemp ]; then
          nseg=$nsegtemp
        fi
      done
    done
  fi

  for ((j=0; j<${nseg}; j++))
  do
    segment=$(printf "%05d" $j)
    out=filelist/${runs[$k]}/rawhit_${segment}.list
    > ${out}

    if [[ "${runspecies}" == "run2pp" ]]; then
      subsystem='TPC'
      for ((i=0; i<24; i++))
      do
        ebdc=$(printf "%02d" $i)
        #echo "${path}/DST_STREAMING_EVENT_${subsystem}${ebdc}/run_000${numlower}_000${numupper}/dst/DST_STREAMING_EVENT_${subsystem}${ebdc}_${runspecies}_${anabuild}_${cdbtag}_${version}-000${run}-${segment}.root" >> ${out}
        echo "DST_STREAMING_EVENT_${subsystem}${ebdc}_${runspecies}_${anabuild}_${cdbtag}_${version}-000${run}-${segment}.root" >> ${out}
      done
    fi

    if [[ "${runspecies}" == "run3auau" || "${runspecies}" == "run3pp" ]]; then
      subsystem='ebdc'
      for ((i=0; i<24; i++))
      do
        ebdc=$(printf "%02d" $i)
        for ((ii=0; ii<2; ii++))
        do
          #echo "${path}/DST_STREAMING_EVENT_${subsystem}${ebdc}_${ii}/run_000${numlower}_000${numupper}/DST_STREAMING_EVENT_${subsystem}${ebdc}_${ii}_${runspecies}_${anabuild}_${cdbtag}_${version}-000${run}-${segment}.root" >> ${out}
          echo "DST_STREAMING_EVENT_${subsystem}${ebdc}_${ii}_${runspecies}_${anabuild}_${cdbtag}_${version}-000${run}-${segment}.root" >> ${out}
        done
      done
    fi

  done

  out=${runs[$k]}_seg_id.txt
  > ${out}

  nseg=`ls filelist/${runs[$k]}/rawhit_*.list | wc -l`
  echo run ${runs[$k]} : ${nseg} dsts, $((${nseg} * ${NumJobPerDst})) jobs

  for ((i=0; i<${nseg}; i++))
  do
    for ((j=0; j<${NumJobPerDst}; j++))
    do
      echo "$i $j" >> ${out}
    done
  done

done

echo "Done"
