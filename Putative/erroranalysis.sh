
OutputName="/home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122"


mergetemp="${OutputName}/results/ProximitiesPlusOverlapPlusBoutonManual1Landmarks12 ${OutputName}/results/ProximitiesPlusOverlapPlusBoutonManual2Landmarks12 /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors stage2man1vs2"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;
mergetemp="${OutputName}/results/ProximitiesPlusOverlapPlusBoutonManual1Landmarks12 ${OutputName}/results/ProximitiesPlusOverlapPlusBoutonAutoLandmarks12 /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors stage2man1vsAuto"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;
mergetemp="${OutputName}/results/ProximitiesPlusOverlapPlusBoutonManual2Landmarks12 ${OutputName}/results/ProximitiesPlusOverlapPlusBoutonAutoLandmarks12 /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors stage2man2vsAuto"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

exit

mergetemp="${OutputName}/results/ProximitiesPlusOverlapManual1Landmarks12 ${OutputName}/results/ProximitiesPlusOverlapManual2Landmarks12 /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors stage1man1vs2"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;
mergetemp="${OutputName}/results/ProximitiesPlusOverlapManual1Landmarks12 ${OutputName}/results/ProximitiesPlusOverlapAutoLandmarks12 /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors stage1man1vsAuto"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;
mergetemp="${OutputName}/results/ProximitiesPlusOverlapManual2Landmarks12 ${OutputName}/results/ProximitiesPlusOverlapAutoLandmarks12 /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/results/errors stage1man2vsAuto"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

exit


# ground truth 1 - mythreya
mergetemp="${OutputName}/manual_overlaps_sure.landmarkAscii ${OutputName}/FixedFinalmerged_crossing_centroids_list.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff crossing"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

mergetemp="${OutputName}/manual_overlaps_unsure.landmarkAscii ${OutputName}/FixedFinalmerged_contact_centroids_list.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff contact"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

mergetemp="${OutputName}/manual_boutons_sure.landmarkAscii ${OutputName}/FixedFinalmerged_bouton_centroids_sure.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff bouton_sure"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

mergetemp="${OutputName}/manual_boutons_unsure.landmarkAscii ${OutputName}/FixedFinalmerged_bouton_centroids_unsure.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff bouton_unsure"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;


# ground truth 2 - david
mergetemp="${OutputName}/manual_overlaps_sure_David.landmarkAscii ${OutputName}/FixedFinalmerged_crossing_centroids_list.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff crossing_david"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

mergetemp="${OutputName}/manual_overlaps_unsure_David.landmarkAscii ${OutputName}/FixedFinalmerged_contact_centroids_list.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff contact_david"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

mergetemp="${OutputName}/manual_boutons_sure_David.landmarkAscii ${OutputName}/FixedFinalmerged_bouton_centroids_sure.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff bouton_sure_david"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

mergetemp="${OutputName}/manual_boutons_unsure_David.landmarkAscii ${OutputName}/FixedFinalmerged_bouton_centroids_unsure.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff bouton_unsure_david"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

# ground truth 1 vs 2 

mergetemp="${OutputName}/manual_overlaps_sure.landmarkAscii ${OutputName}/manual_overlaps_sure_David.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff crossing_1_vs_2"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

mergetemp="${OutputName}/manual_overlaps_unsure.landmarkAscii ${OutputName}/manual_overlaps_unsure_David.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff contact_1_vs_2"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

mergetemp="${OutputName}/manual_boutons_sure.landmarkAscii ${OutputName}/manual_boutons_sure_David.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff bouton_sure_1_vs_2"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;

mergetemp="${OutputName}/manual_boutons_unsure.landmarkAscii ${OutputName}/manual_boutons_unsure_David.landmarkAscii /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff /home/mythreya/Documents/outputprox/unsorted_dist10/DS_20140122/error_stuff bouton_unsure_1_vs_2"
echo $mergetemp
/home/mythreya/projects/ErrorCalulator/ErrorCalculator $mergetemp;


# sure + unsure merged