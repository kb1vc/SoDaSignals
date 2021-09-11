#!/bin/bash

i=64
j=2

echo '#size trials elapsed seconds-per-pt' > res.dat

for ((p2 = 4; p2 < 18; p2++))
do
    for ((p3 = 0; p3 < 3; p3++))
    do
	for ((p5 = 0; p5 < 1; p5++))
	do
	    size=$((2**${p2} * 3**${p3} * 5**${p5}))
	    if [ ${size} -lt 100000 ]
	    then
		s=$((${p2} + ${p3} + ${p5}))
		p=$((${p2} * ${p3} * ${p5}))
	        echo -n "${p2} ${p3} ${p5} $s $p " >> res.dat
  	        ./FFTTiming 100 ${size} >> res.dat
	        echo ${size}
	    fi
	done
    done
done

sort -n res.dat > sres.dat
