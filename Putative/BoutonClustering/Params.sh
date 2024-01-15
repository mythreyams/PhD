#! /bin/bash

pathName=$1
cellName="DS_20140122"
proxName="ProximityDetection"
dataName="${pathName}/${cellName}"
outputpath=$4

for i in $(seq -f "%03g" $2 $3)
do

#sectionFolders=find ${dataName}/S${i}/ProximityDetection/ -name "proximity_roi*" -type d
#($(ls ${dataName} | grep S${i}))
#echo $sectionFolders
 #find ${dataName}/S${i}/ProximityDetection/ -name "proximity_roi*" -type d 
#"${sectionFolders[@]}"
        #cd ${dataName}/S${i}_*
        #pwd
        len=`expr length ${dataName}`
	echo $len
        
        for j in `find ${dataName}/S${i}_* -type d`
		    
	do
		cd $j
		pwd
		
		
		# count the number of images within this region before invoking the program
		file=0
		for k in `find $j -name "roi*.png"`
			
		do
		    file=`expr $file + 1`
		    
		done
		
		file=`expr $file - 1`
		
		echo "$file"
		#roinum="0"
		#roinum= {printf $j | cut -c `expr $len + 7`-`expr $len + 8`}
		#roinum=`expr $roinum + 1`
		#echo "$roinum"
		#echo $roinum

		xz=`echo "$j" | cut -c 58- `
		echo $xz
		# Call the program
		/home/mythreya/projects/BoutonParameterEstimation/BoutonFinder "/home/mythreya/ResampledSpatialgraph.am" "${j}" 0 ${file} "${j}/roi_load.hx" ${outputpath} "${i}" "${xz}"  
		
		    
		    
		    
	done
       
       
        #find ${dataName}/S${i}/ProximityDetection/ -name "proximity_roi*" -type d
        #proxFolders=$(ls ${dataName}/S${i}/ProximityDetection)
        #echo result
           

done


#cell_spatial_graph="/network_home/mythreya/Documents/output/S05_S28/DS_20140122/S005_S028_merged_finished_test2.am"
#proximity_image_folder="/network_home/mythreya/Documents/output/S05_S28/DS_20140122/S016/ProximityDetection/proximity_roi3/"
#hxfile="/network_home/mythreya/Documents/output/S05_S28/DS_20140122/S016/ProximityDetection/proximity_roi3/loadROI.hx"
#outputpath="/network_home/mythreya/Documents/output/S05_S28/"
#tmp="${cell_spatial_graph} ${proximity_image_folder} 0 48 ${hxfile} ${outputpath} "
#network_home/mythreya/projects/BoutonParameterEstimation/BoutonFinder $tmp;




