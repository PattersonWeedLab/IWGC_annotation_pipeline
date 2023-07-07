#!/bin/bash 

list="$1"

wc -l ${list} > ${list}.log 
counter=1
for x in $( sed 's/ /#/g' ${list} )
do 
echo $( echo ${x} | sed 's/#/ /g' )
$( echo ${x} | sed 's/#/ /g' ) 
echo finished ${counter}
counter=$(( ${counter} + 1 )) 
done

