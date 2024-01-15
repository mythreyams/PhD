#! /bin/bash

pathName="/network_home/mythreya/Documents/output/S05_S28"
cellName="DS_20140122"
proxName="ProximityDetection"
dataName="${pathName}/${cellName}"
outputpath="/network_home/mythreya/Documents/output/S05_S28/"

for i in $(seq -f "%03g" $1 $2)
do

#sectionFolders=find ${dataName}/S${i}/ProximityDetection/ -name "proximity_roi*" -type d
#($(ls ${dataName} | grep S${i}))
#echo $sectionFolders
 #find ${dataName}/S${i}/ProximityDetection/ -name "proximity_roi*" -type d 
#"${sectionFolders[@]}"
        cd ${dataName}/S${i}/ProximityDetection/
        pwd
       
        roinum=-1
        for j in `find ${dataName}/S${i}/ProximityDetection/ -name "proximity_roi*" -type d`
		    
	do
		cd $j
		pwd
		
		roinum=`expr $roinum + 1`
		
		# count the number of images within this region before invoking the program
		file=0
		for k in `find $j -name "proximity*.tif"`
			
		do
		    file=`expr $file + 1`
		    
		done
		
		file=`expr $file - 1`
		
		echo "$file"
		
		# Call the program
		/network_home/mythreya/projects/BoutonParameterEstimation/BoutonFinder "${dataName}/S005_S028_merged_finished_test2.am" "${j}" 0 ${file} "${j}/loadROI.hx" ${outputpath} "${i}" "$roinum"  
		
		    
		    
		    
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




