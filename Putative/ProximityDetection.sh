#! /bin/bash

cellName="DS_20140122"
mergeddir="Merged"
transformdir="Ready_to_merge"
mergedfile="final.am"
pathName="/katz/NeuroMorph/Final_Katz_Confocal"

#0: test without loading images
#1: normal
noImageTest="1"

#mergedfile=$4;
#name=$2
#dataName="${pathName}/${cellName}"
dataName=$1
OutputName="Proximities"
OutputPath="${dataName}/${OutputName}"
echo OutputPath

exit
# for l in 0;
# do
#   for k in 4;
#   do
OutputName="/home/mythreya/Documents/outputprox/${name}/${cellName}"
mergedFilename="/katz/NeuroMorph/Final_Katz_Confocal/${cellName}/${mergeddir}/${mergedfile}"
if [ ! -e $mergedFilename ]
    then
    echo "File $mergedFilename does not exist!"
    exit 0
fi
#z_scale_dir="/katz/NeuroMorph/Final_Katz_Confocal/${cellName}/scale_offset/"


mkdir "${dataName}/${OutputName}"
chmod 777 "${dataName}/${OutputName}"
mkdir "/home/mythreya/Documents/outputprox/${name}/"
chmod 777 "/home/mythreya/Documents/outputprox/${name}/"
mkdir ${OutputName}
chmod 777 ${OutputName}
cp -v $mergedFilename $OutputName

tmp="${mergedFilename} ${OutputName} 4"
echo ${tmp}
/home/mythreya/projects/AxonDendriteProximityDetector/AxonDendriteProximityFinder $tmp;



for i in $(seq -f "%03g" $3 $4)
do	

	sectionFolders=($(ls ${dataName} | grep S${i}))
	for secDir in "${sectionFolders[@]}"
	do
	    echo $secDir
	    transformated_file="${dataName}/${transformdir}/S${i}*.am"
	    echo ${transformated_file}
	    if [ ! -e $transformated_file ]
		then
		echo "File $transformated_file does not exist!"
		exit 0
	    fi

	    
	    #mkdir "${OutputName}/S${i}"
	    #chmod -R 777 "${OutputName}/S${i}"
	    
    #	cd "${z_scale_dir}"
	    
	    
    # 	if `find -name "S"${i}"_Zscaled_*_Trans.am"`
    # 	
    # 	then
	    #z_scale_filename="${z_scale_dir}/S${i}.txt"
	    #z_scale_filename= basename "${z_scale_dir}S${i}_Zscaled_*.am"
	    #echo ${z_scale_filename}
	    #fi
	    
	    
	    
	    #cd "${OutputName}"
	    #chmod -R 777 "S${i}/"
	    
	    cd "${dataName}/$secDir"
	    
	    
	    
	    
    # 	

	    file=0
	    for j in `find $dataName/$secDir/brick01/ -name "deconvolved_*.png"`
		    
	    do
		file=`expr $file + 1`
		
	    done
	    
	    file=`expr $file - 1`
	    
	    echo "$file"

	    
	    
	    
	    
	    
	    
	    #mkdir "${OutputName}/S${i}/ProximityDetection"
	    #chmod 777 "${OutputName}/S${i}/ProximityDetection"

	    
	    #cd "${dataName}/$secDir"
	    #echo "changing directory to ${dataName}/$secDir "


    
	    #Proper format for whole cell:  ./ProximityFinder 'amira_file_with_transformations.am' 'deconvolved_image%03d.png' 'output_name' startIndex endIndex includeUnknownSegments  'amira_file_whole_cell.am'  'nondecon_images.png'
	    #tmp=" ${dataName}/$secDir ${i} 0 ${file} ${OutputName}/S${i}/Rois ${mergedFilename} ${transformated_file}"
	    tmp="${transformated_file} deconvolved_z%03d.png ${OutputName} 0 ${file} 0 ${mergedFilename} $noImageTest ${OutputName}/Proximity_Landmarks.landmarkAscii ${OutputName}/Proximities.csv ${OutputName}/Proximities1.csv ${i}"
	    echo ${tmp}
	    /home/mythreya/projects/AxonDendriteProximityImagePublisher/AxonDendriteProximityImages $tmp;
	    \cp ${OutputName}/Proximities1.csv ${OutputName}/Proximities.csv
	    
	    
# 	    # The output of the prox detection is arranged in each prox zone 
# 	    cd "${OutputName}/S${i}/ProximityDetection"
# 	    echo "changing directory to ${OutputName}/S${i}/ProximityDetection "
# 	    
# 	    roicnt=0
# 	    
# 	    # go inside each prox folder and call the image count thingy
# 	    proxFolders=($(ls -d */ ))
# 	    for proxDir in "${proxFolders[@]}"
# 	    do
# 	      echo ${OutputName}/S${i}/ProximityDetection/$proxDir
# 	      
# 	      cd ${OutputName}/S${i}/ProximityDetection/$proxDir
# 	      
# 	      prox_file=0
# 	      axon_num=${proxDir:5:1}
# 	      dend_num=${proxDir:12:1}
# 	      echo "$axon_num"
# 	      echo "$dend_num"
# 	      
# 	      for k in `find ${OutputName}/S${i}/ProximityDetection/$proxDir -name "proximity*.tif"`
# 		      
# 	      do
# 		  prox_file=`expr $prox_file + 1`
# 		  
# 	      done
# 	      
# 	      prox_file=`expr $prox_file - 1`
# 	      
# 	      echo "$prox_file"
# 	      
# 	      #proxtmp="$mergedFilename ${OutputName}/S${i}/ProximityDetection/$proxDir 0 $prox_file ${OutputName}/S${i}/ProximityDetection/$proxDirloadROI.hx ${OutputName}/S${i}/ProximityDetection/$proxDir ${OutputName}/S${i}/ProximityDetection/$proxDirLandmarks.landmarkAscii"
# 	      proxtmp="$mergedFilename ${OutputName}/S${i}/ProximityDetection/$proxDir 0 $prox_file ${OutputName}/S${i}/ProximityDetection/${proxDir}loadROI.hx \
# 	      ${OutputName}/S${i}/ProximityDetection/$proxDir ${OutputName} \
# 	      ${i} ${roicnt} ${axon_num} ${dend_num}  
# 	      ${OutputName}/S${i}/ProximityDetection/${proxDir}manual_bouton_landmarks.landmarkAscii ${OutputName}/S${i}/ProximityDetection/${proxDir}manual_overlap_landmarks.landmarkAscii"
# 	      echo $proxtmp	
# 	      
# 	      #/home/mythreya/projects/BoutonDetector/BoutonFinder $proxtmp;
# 	      
# 	      roicnt=`expr $proicnt + 1`
# 	      
# 	    
# 	    done


#done

# done
done

done