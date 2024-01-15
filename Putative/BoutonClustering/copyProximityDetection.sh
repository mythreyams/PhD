#! /bin/bash
cellName="DS_20140122"
mergeddir="Merged"
transformdir="ZScaled_Trans"
mergedfile="S005_S028_merged_finished_test2.am"
pathName="/katz/NeuroMorph/Final_Katz_Confocal"

#0: test without loading images
#1: normal
noImageTest=$4

#mergedfile=$4;
name=$1
dataName="${pathName}/${cellName}"

# for l in 0;
# do
#   for k in 4;
#   do
OutputName="/nas1/Data_NeuroMorph/confocal_tracing/SynapseDetection/${name}/${cellName}/"
mergedFilename="/katz/NeuroMorph/Final_Katz_Confocal/${cellName}/${mergeddir}/${mergedfile}"
if [ ! -e $mergedFilename ]
    then
    echo "File $mergedFilename does not exist!"
    exit 0
fi
#z_scale_dir="/katz/NeuroMorph/Final_Katz_Confocal/${cellName}/scale_offset/"



mkdir "/nas1/Data_NeuroMorph/confocal_tracing/SynapseDetection/${name}/"
chmod 777 "/nas1/Data_NeuroMorph/confocal_tracing/SynapseDetection/${name}/"
mkdir ${OutputName}
chmod 777 ${OutputName}
cp -v $mergedFilename $OutputName





for i in $(seq -f "%03g" $2 $3)
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

	    
	    mkdir "${OutputName}/S${i}"
	    chmod -R 777 "${OutputName}/S${i}"
	    
    #	cd "${z_scale_dir}"
	    
	    
    # 	if `find -name "S"${i}"_Zscaled_*_Trans.am"`
    # 	
    # 	then
	    #z_scale_filename="${z_scale_dir}/S${i}.txt"
	    #z_scale_filename= basename "${z_scale_dir}S${i}_Zscaled_*.am"
	    #echo ${z_scale_filename}
	    #fi
	    
	    
	    
	    cd "${OutputName}"
	    chmod -R 777 "S${i}/"
	    
	    cd "${dataName}/$secDir"
	    
	    
	    
	    
    # 	

	    file=0
	    for j in `find $dataName/$secDir/brick01/ -name "deconvolved_*.png"`
		    
	    do
		file=`expr $file + 1`
		
	    done
	    
	    file=`expr $file - 1`
	    
	    echo "$file"

	    
	    
	    
	    
	    
	    
	    mkdir "${OutputName}/S${i}/ProximityDetection"
	    chmod 777 "${OutputName}/S${i}/ProximityDetection"

	    
	    cd "${dataName}/$secDir"
	    echo "changing directory to ${dataName}/$secDir "


    
#Proper format for whole cell:  ./ProximityFinder 'amira_file_with_transformations.am' 'deconvolved_image%03d.png' 'output_name' startIndex endIndex includeUnknownSegments  'amira_file_whole_cell.am'  'nondecon_images.png'

tmp="${transformated_file} deconvolved_z%03d.png ${OutputName}/S${i}/ProximityDetection/proximity 0 ${file} 1 ${mergedFilename} $noImageTest"
echo ${tmp}
/home/mythreya/Documents/AxonDendriteProximityDetector/AxonDendriteProximityFinder $tmp;

#done

# done
done

done