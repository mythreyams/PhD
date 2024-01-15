#!/bin/bash

input="$PWD/Proximity_Zone_33/Proximity_Zone_33000.tif"

seedX="194"
seedY="409"
intidist="1"
alivedist="1"
outsidedist="1"
sigma="1"
alpha="0.5"
beta="0.6"
propscaling=( "-10" "-1" "0.1" "0" "0.1" "1" "10" )
advecscaling=( "-10" "-1" "0.1" "0" "0.1" "1" "10" )
curvscaling=( "-10" "-1" "0.1" "0" "0.1" "1" "10" )
maxrmserror="0.02"
maxiterations="10000"
geoOrShape="1"
trialoffset="3"
outsideoffset="50"
i="1"
j="0.1"
k="-1"
#for i in "${propscaling[@]}"
#do

#for j in "${advecscaling[@]}"
#do

#for k in "${curvscaling[@]}"
#do

output="$PWD/outputpics/output_p${i}_a${j}_c${k}.png"
echo "$input $output $seedX $seedY $intidist $sigma $alpha $beta ${i} ${j} ${k} $maxrmserror $maxiterations $geoOrShape $alivedist $outsidedist"
$PWD/Surface $input $output $seedX $seedY $intidist $sigma $alpha $beta ${i} ${j} ${k} $maxrmserror $maxiterations $geoOrShape $alivedist $outsidedist $trialoffset $outsideoffset

#done
#done
#done

