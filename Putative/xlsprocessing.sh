#! /bin/bash
#landmarkfilename="/home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/manual_overlaps_sure_David.landmarkAscii"
landmarkfilename=$2
OutputFolderName=$1
OutputFileName=$4
#OutputName="/home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122"
local_line_num=0

while read line
do 
#echo "line=$line"
    local_line_num=0
    putative_centroidX=$( echo $line | cut -d$' ' -f1 )
    putative_centroidX=`echo ${putative_centroidX} | sed -e 's/[eE]+*/\\*10\\^/'`
    putative_centroidY=$( echo $line | cut -d$' ' -f2 )
    putative_centroidY=`echo ${putative_centroidY} | sed -e 's/[eE]+*/\\*10\\^/'`
    putative_centroidZ=$( echo $line | cut -d$' ' -f3 )
    putative_centroidZ=`echo ${putative_centroidZ} | sed -e 's/[eE]+*/\\*10\\^/'`
    #echo $putative_centroidX $putative_centroidY $putative_centroidZ
    min_dist=1000
    linenumerselected=0

    # putatives xls
    while read local_line
    do
    #echo $local_line
	  local_line_num=`expr $local_line_num + 1 `
	  prox_centroidX=$( echo $local_line | cut -d$',' -f2 )
	  prox_centroidY=$( echo $local_line | cut -d$',' -f3 )
	  prox_centroidZ=$( echo $local_line | cut -d$',' -f4 )

	  #echo $prox_centroidX $prox_centroidY $prox_centroidZ
	  
	  dist=$(echo "sqrt ( ( $putative_centroidX - $prox_centroidX ) ^2 + ( $putative_centroidY - $prox_centroidY ) ^2 + ( $putative_centroidZ - $prox_centroidZ ) ^2 )" | bc -l )
	  
	  echo $dist
	  if [ 1 -eq "$(echo "$dist < $min_dist" | bc -l)" ]
	  then
	    
	    
	      if [ 1 -eq "$(echo "$dist < 20" | bc -l)" ]
	      then
	      
	    
		min_dist=$dist

		
		linenumerselected=$local_line_num

		minX=$prox_centroidX
		minY=$prox_centroidY
		minZ=$prox_centroidZ
		#sed "/$minX $minY $minZ/s/$/ 1 $putative_centroidX $putative_centroidY $putative_centroidZ /" "${OutputName}/Proximities1.csv"
	      fi
	  
	  fi
    done < ${OutputFolderName}/Proximities1.csv

local_line_num=0



# this is for adding to the found line number (this is so lame)
  # putatives xls
    while read local_line
    do
    #echo $local_line
	  local_line_num=`expr $local_line_num + 1 `
	  
	  if [ 1 -eq "$(echo "$local_line_num == $linenumerselected" | bc -l)" ]
	  then
	    
	      echo "${local_line},$3,$putative_centroidX,$putative_centroidY,$putative_centroidZ,"  >> ${OutputFileName}
	 
	  
	  fi
    done < ${OutputFolderName}/Proximities1.csv 
   
done < $landmarkfilename 
