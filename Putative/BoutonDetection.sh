#! /bin/bash


#cell_spatial_graph="/nas1/Data_NeuroMorph/confocal_tracing/SynapseDetection/S016_test_240415_6/DS_20140122/S005_S028_merged_finished_test2.am"
#cell_spatial_graph="/nas1/Data_NeuroMorph/confocal_tracing/SynapseDetection/S016_test_240415_6/DS_20140122/S016/ProximityDetection/proximitycell_3.am"
#proximity_image_folder="/nas1/Data_NeuroMorph/confocal_tracing/SynapseDetection/04_29_15_1/DS_20140122/S016/ProximityDetection/proximity_CELL#1-->CELL#3_roi3"
#hxfile="/nas1/Data_NeuroMorph/confocal_tracing/SynapseDetection/04_29_15_1/DS_20140122/S016/ProximityDetection/proximity_CELL#1-->CELL#3_roi3/loadROI.hx"
cell_spatial_graph="/home/neuromorph/projects/input_data/04_29_15_1/DS_20140122/S016/ProximityDetection/proximitycell_3.am"
proximity_image_folder="/home/neuromorph/projects/input_data/04_29_15_1/DS_20140122/S016/ProximityDetection/proximity_CELL#1-->CELL#3_roi3"
hxfile="/home/neuromorph/projects/input_data/04_29_15_1/DS_20140122/S016/ProximityDetection/proximity_CELL#1->CELL#3_roi3/loadROI.hx"
outputname="/home/neuromorph/projects/output_data/$1"
tmp="${cell_spatial_graph} ${proximity_image_folder} 0 46 ${hxfile} $outputname"
/home/neuromorph/projects/BoutonFinderLocal/BoutonFinder $tmp;

