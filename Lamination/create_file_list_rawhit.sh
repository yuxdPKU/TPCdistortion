#!/bin/bash

get_closest_numbers() {
    number=$1
    lower=$(( (number / 100) * 100 ))
    upper=$(( lower + 100 ))
    echo "$lower $upper"
}

runspecies='run2auau'
runtype='physics'
anabuild='ana464'
cdbtag='nocdbtag'
version='v001'

path=/sphenix/lustre01/sphnxpro/production/${runspecies}/${runtype}/${anabuild}_${cdbtag}_${version}
echo Searching in ${path}

runs=(54966 54967 54968 54969)

NumEvtPerDst=10000
NumEvtPerJob=10000
NumJobPerDst=$((${NumEvtPerDst} / ${NumEvtPerJob}))

for ((k=0; k<${#runs[@]}; k++))
do

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

#  subsystem='INTT'
#  for ((server=0; server<8; server++))
#  do
#    nsegtemp=`ls ${path}/DST_STREAMING_EVENT_${subsystem}${server}/run_000${numlower}_000${numupper}/dst/*${runs[$k]}* | wc -l`
#    echo run ${runs[$k]} ${subsystem}${server} raw hit DST ${nsegtemp} segments
#    if [ $nseg -gt $nsegtemp ]; then
#      nseg=$nsegtemp
#    fi
#  done
#
#  subsystem='MVTX'
#  for ((felix=0; felix<6; felix++))
#  do
#    nsegtemp=`ls ${path}/DST_STREAMING_EVENT_${subsystem}${felix}/run_000${numlower}_000${numupper}/dst/*${runs[$k]}* | wc -l`
#    echo run ${runs[$k]} ${subsystem}${felix} raw hit DST ${nsegtemp} segments
#    if [ $nseg -gt $nsegtemp ]; then
#      nseg=$nsegtemp
#    fi
#  done

  subsystem='TPC'
  for ((i=0; i<24; i++))
  do
    ebdc=$(printf "%02d" $i)
    nsegtemp=`ls ${path}/DST_STREAMING_EVENT_${subsystem}${ebdc}/run_000${numlower}_000${numupper}/dst/*${runs[$k]}* | wc -l`
    echo run ${runs[$k]} ${subsystem}${ebdc} raw hit DST ${nsegtemp} segments
    if [ $nseg -gt $nsegtemp ]; then
      nseg=$nsegtemp
    fi
  done

#  subsystem='TPOT'
#  nsegtemp=`ls ${path}/DST_STREAMING_EVENT_${subsystem}/run_000${numlower}_000${numupper}/dst/*${runs[$k]}* | wc -l`
#  echo run ${runs[$k]} ${subsystem} raw hit DST ${nsegtemp} segments
#  if [ $nseg -gt $nsegtemp ]; then
#    nseg=$nsegtemp
#  fi

  for ((j=0; j<${nseg}; j++))
  do
    segment=$(printf "%05d" $j)
    out=filelist/${runs[$k]}/rawhit_${segment}.list
    > ${out}

#    subsystem='INTT'
#    for ((server=0; server<8; server++))
#    do
#      echo "${path}/DST_STREAMING_EVENT_${subsystem}${server}/run_000${numlower}_000${numupper}/dst/DST_STREAMING_EVENT_${subsystem}${server}_${runspecies}_${anabuild}_${cdbtag}_${version}-000${runs[$k]}-${segment}.root" >> ${out}
#    done
#
#    subsystem='MVTX'
#    for ((felix=0; felix<6; felix++))
#    do
#      echo "${path}/DST_STREAMING_EVENT_${subsystem}${felix}/run_000${numlower}_000${numupper}/dst/DST_STREAMING_EVENT_${subsystem}${felix}_${runspecies}_${anabuild}_${cdbtag}_${version}-000${runs[$k]}-${segment}.root" >> ${out}
#    done

    subsystem='TPC'
    for ((i=0; i<24; i++))
    do
      ebdc=$(printf "%02d" $i)
      echo "${path}/DST_STREAMING_EVENT_${subsystem}${ebdc}/run_000${numlower}_000${numupper}/dst/DST_STREAMING_EVENT_${subsystem}${ebdc}_${runspecies}_${anabuild}_${cdbtag}_${version}-000${runs[$k]}-${segment}.root" >> ${out}
    done

#    subsystem='TPOT'
#    echo "${path}/DST_STREAMING_EVENT_${subsystem}/run_000${numlower}_000${numupper}/dst/DST_STREAMING_EVENT_${subsystem}_${runspecies}_${anabuild}_${cdbtag}_${version}-000${runs[$k]}-${segment}.root" >> ${out}

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
