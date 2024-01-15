%% This is the main registration script and it calls relevant routines to 
% perform  the following  
% 1) Read FN contours and cell locations to be registered from the Amira
% files
% 2) Register these to the Avg FN reference frame
% 3) Write the registered FN along with axes, registred cells and 
% related stats.

% This scripts reads the input contour, cell location and morphology files 
% in the Amira format. User must change them to use the input and ouput file 
% formats specific to the project. Output of the registration routine is a 
% n*3 matrix (where n is the number of data points) and does not depend on 
% the depend on the visualization tool used such as Amira. This script reads
% inputs into raw data and provides it for the registration routine and
% susequently converts the output of the registration routine for
% visualizatoin (in this case Amira).

%%%%%%%%%%%%%%%%%%%%%% Caution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order for this script to produce results consistent with the paper
% the following directions need to be maintained for inputs. In cases where
% the input directions are not as expected the output can be converted
% appropriately (see code below)

% X Axis : Lateral - Medial
% Y Axis : Rostral - Caudal
% Z Axis : Ventral - Dorsal
%%%%%%%%%%%%%%%%%%%%%% Caution End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Input and Output Paths                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up the paths
FN_Contour_File = strcat(pwd, '\..\..\TestInputs\FN_Contours.am');
FN_Cells_to_be_registered = strcat(pwd, '\..\..\TestInputs\FN_VMCells_C3.landmarkAscii');

output_folder = strcat(pwd, '\..\..\TestOutputs\');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Avg Axes and radii of the FN reference frame definitions            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% canonical coordinate directions
% In Amira visualization software, the Y coordinate goes from top to bottom.
% therefore flip Y axis
canonical_coordinate_directions = [1,-1,1];

% canonical axes directions
% flip Z axis to maintain the preferred positive directions in Amira
canonical_axes_directions = [1,1,-1];

% avg FN radii
avg_fn_radii = [806.4708058 634.2044403 430.7335172];
% avg FN axes
avg_fn_axes =  [ -0.7499    0.5258    0.4015;...
                  0.6137    0.7795    0.1254;...
                  0.2470   -0.3404    0.9072 ];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coordinate Axes Conversions                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XAxisCanonical = ConvertAmiraToCanonicalCoordinates([1 0 0]);
YAxisCanonical = ConvertAmiraToCanonicalCoordinates([0 1 0]);
ZAxisCanonical = ConvertAmiraToCanonicalCoordinates([0 0 1]);

% Write axes in the Canonical coordinate system (invert Y axis since it points downwards in Amira)
writeSpatialGraphLineAmira([0,0,0],XAxisCanonical*avg_fn_radii(1), strcat(output_folder, '0_X-Canonical-LM.am'),0, 0,canonical_axes_directions);
writeSpatialGraphLineAmira([0,0,0],YAxisCanonical*avg_fn_radii(1), strcat(output_folder, '0_Z-Canonical-VD.am'),2, 0,canonical_axes_directions);
writeSpatialGraphLineAmira([0,0,0],ZAxisCanonical*avg_fn_radii(1), strcat(output_folder, '0_Y-Canonical-RC.am'),5, 0,canonical_axes_directions);

% Write axes of the Avg FN
writeSpatialGraphLineAmira([0,0,0],avg_fn_axes(:,1)*avg_fn_radii(1), strcat(output_folder, '4_AvgFN_Y-Prime.am'),0, 0,canonical_coordinate_directions);
writeSpatialGraphLineAmira([0,0,0],avg_fn_axes(:,2)*avg_fn_radii(2), strcat(output_folder, '4_AvgFN_X-Prime.am'),2, 0,canonical_coordinate_directions);
writeSpatialGraphLineAmira([0,0,0],avg_fn_axes(:,3)*avg_fn_radii(3), strcat(output_folder, '4_AvgFN_Z-Prime.am'),5, 0,canonical_coordinate_directions);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read FN data from different experiments                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read the FN contours and cells coordinates
% there is a tag in the Amira spatial graph file which indicates different 
% sections of the file indicated by '@'. Points begin at @4
[spatialgraph_points_original] =  readPtsSpatialGraphAmira(FN_Contour_File);

[cells_original] = readLandmarkFileAmira( FN_Cells_to_be_registered );

% Change axis and direction: Interchange Y and Z axis 
% Axes in the original data: 
% X - Lateral -> Medial, 
% Y - Dorsal -> Ventral,
% Z - Caudal -> Rostral
% 
% Required Directions:
% X - Lateral -> Medial, 
% Y - Rostral -> Caudal,
% Z - Ventral -> Dorsal.
%
% Note : In Amira Y axis goes from top to bottom thereofore alwars flip y
% axis to be bottom to top
spatialgraph_points_canonical = ConvertAmiraToCanonicalCoordinates(spatialgraph_points_original);
spatialgraph_points_centered_canonical = [spatialgraph_points_canonical(:,1)-mean(spatialgraph_points_canonical(:,1))...
                                          spatialgraph_points_canonical(:,2)-mean(spatialgraph_points_canonical(:,2))...
                                          spatialgraph_points_canonical(:,3)-mean(spatialgraph_points_canonical(:,3))];
%spatialgraph_coeff = pca(spatialgraph_points_canonical);
cells_canonical = ConvertAmiraToCanonicalCoordinates(cells_original);
cells_centered_canonical = [(cells_canonical(:,1)-mean(spatialgraph_points_canonical(:,1)))...
                            (cells_canonical(:,2)-mean(spatialgraph_points_canonical(:,2)))...
                            (cells_canonical(:,3)-mean(spatialgraph_points_canonical(:,3))) ];

% write data in canonical form
writeLandmarkFileAmira( spatialgraph_points_centered_canonical ,...
    strcat(output_folder, '1_Spatialgraph_canonical.landmarkAscii'), canonical_coordinate_directions);
writeLandmarkFileAmira( cells_centered_canonical ,...
    strcat(output_folder, '1_Centered_Cells_canonical.landmarkAscii'),canonical_coordinate_directions);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fit an ellipsoid and fetch the its axes and radii                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit an Ellipsoid to the contours
[center, fn_radii, fn_axes, FN_surface_points] = fitEllipsoid(spatialgraph_points_centered_canonical);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Register the FN ellipsoid and the cells into Avg FN reference frame %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotate the ellipsoid into its local reference frame by inverting rotation
% and scale it along the Avg FN radii; thus registering it.
[centered_FN_surface_local,rotated_FN_surface_local,scaled_FN_surface_local,fn_evecs] = ...
    rigidTransformation(FN_surface_points,center, fn_axes, avg_fn_radii, fn_radii, 0,1,1,0);

% Rotate the registered FN to the avg_fn axes so that it can be viewed in
% canonical form
[centered_FN_surface_canonical,rotated_FN_surface_canonical,scaled_FN_surface_canonical,fn_evecs] = ...
    rigidTransformation(scaled_FN_surface_local,center, avg_fn_axes, avg_fn_radii, fn_radii, 0,0,0,0);

% Rotate the cells to the local axis of its FN ellipsoid
[centered_cells_local,rotated_cells_local,scaled_cells_local,fn_evecs] = ...
    rigidTransformation(cells_centered_canonical,center, fn_axes, avg_fn_radii, fn_radii, 0,1,1,0);

% Rotate the registered cells to the avg_fn axes so that it can be viewed in
% canonical form
[centered_cells_canonical,rotated_cells_canonical,scaled_cells_canonical,fn_evecs] = ...
    rigidTransformation(scaled_cells_local,center, avg_fn_axes, avg_fn_radii, fn_radii, 0,0,0,0);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write the registered FN(along with axes) and cells                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write the Centered and Registered FN and Cells as landmark files
% Centered FN ellipsoid
writeLandmarkFileAmira( centered_FN_surface_local ,...
    strcat(output_folder, '2_Centered_FN_Canonical.landmarkAscii'),canonical_coordinate_directions);

% Registered FN surface in Avg FN coordinate system
% Lets visualize the longest axis as Y in the Local cordinate system
% Therefore interchange X and Y coordinates
scaled_FN_surface_local_vertical = [scaled_FN_surface_local(:,2) scaled_FN_surface_local(:,1) scaled_FN_surface_local(:,3)]; 
writeLandmarkFileAmira( scaled_FN_surface_local_vertical ,...
    strcat(output_folder, '3_Registered_FN_surface_local.landmarkAscii'),canonical_coordinate_directions);

% Registered Cells in Avg FN coordinate system
% Lets visualize the longest axis as Y in the Local cordinate system
% Therefore interchange X and Y coordinates
scaled_cells_local_vertical = [scaled_cells_local(:,2) scaled_cells_local(:,1) scaled_cells_local(:,3)]; 
writeLandmarkFileAmira( scaled_cells_local_vertical ,...
    strcat(output_folder, '3_Registered_Cells_local.landmarkAscii'),canonical_coordinate_directions);

% Registered FN surface in Canonical coordinate system
writeLandmarkFileAmira( rotated_FN_surface_canonical ,...
    strcat(output_folder, '4_Registered_FN_surface_Canonical.landmarkAscii'),canonical_coordinate_directions);

% Registered Cells in Canonical coordinate system
writeLandmarkFileAmira( rotated_cells_canonical ,...
    strcat(output_folder, '4_Registered_Cells_Canonical.landmarkAscii'),canonical_coordinate_directions);


% Write axes in the Canonical coordinate system 
writeSpatialGraphLineAmira([0,0,0],fn_axes(:,1)*fn_radii(1), strcat(output_folder, '2_FN_Y-Prime.am'),0, 0,canonical_coordinate_directions);
writeSpatialGraphLineAmira([0,0,0],fn_axes(:,2)*fn_radii(2), strcat(output_folder, '2_FN_X-Prime.am'),2, 0,canonical_coordinate_directions);
writeSpatialGraphLineAmira([0,0,0],fn_axes(:,3)*fn_radii(3), strcat(output_folder, '2_FN_Z-Prime.am'),5, 0,canonical_coordinate_directions);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write the center, axes direction and axes extent information as a csv%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% write stats about this ellipsoid fit
writeFNEllipsoidStats(scaled_FN_surface_local,output_folder,'01' ,center,fn_axes,fn_radii);
    
% write stats about the registered cells
writeRegisteredCellStats(scaled_cells_local,output_folder,'01','C3');





    



