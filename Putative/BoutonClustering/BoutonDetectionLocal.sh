#! /bin/bash


#cell_spatial_graph="/nas1/Data_NeuroMorph/confocal_tracing/SynapseDetection/S016_test_240415_6/DS_20140122/S005_S028_merged_finished_test2.am"
cell_spatial_graph="/network_home/mythreya/Documents/output/S05_S28/DS_20140122/S005_S028_merged_finished_test2.am"
proximity_image_folder="/network_home/mythreya/Documents/output/S05_S28/DS_20140122/S016/ProximityDetection/proximity_roi3/"
hxfile="/network_home/mythreya/Documents/output/S05_S28/DS_20140122/S016/ProximityDetection/proximity_roi3/loadROI.hx"
outputpath="/network_home/mythreya/Documents/output/S05_S28/"
tmp="${cell_spatial_graph} ${proximity_image_folder} 0 48 ${hxfile} ${outputpath} "
/network_home/mythreya/projects/BoutonParameterEstimation/BoutonFinder $tmp;

