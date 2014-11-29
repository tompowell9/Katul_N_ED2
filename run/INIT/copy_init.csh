#!/bin/bash
#generate monoculture css and pps files based on New_PFT css and pss files 

#define the PFT
PFT=(23 24 26 27 29)

for i_PFT in $(seq 0 4)
do
    #define prefix
    prefix="PFT_${PFT[$i_PFT]}_PV_2009_patch"
	suffix=".lat10.4lon-85.4"

	# loop over the patches
	for i_patch in $(seq 1 6)
	do
		#copy pss files
		cp "./New_PFT_PV_2009_patch$i_patch$suffix.pss" "./$prefix$i_patch$suffix.pss"

		#copy css files
		cp "./New_PFT_PV_2009_patch$i_patch$suffix.css" "./$prefix$i_patch$suffix.css"

		#change the PFT
		for org_PFT in $(seq 0 4)
		do
			sed -i "s: ${PFT[$org_PFT]} : ${PFT[$i_PFT]} :" "./$prefix$i_patch$suffix.css" 
		done
		
		#copy ED2IN files

	done
#echo ${PFT[$i_PFT]}
done
