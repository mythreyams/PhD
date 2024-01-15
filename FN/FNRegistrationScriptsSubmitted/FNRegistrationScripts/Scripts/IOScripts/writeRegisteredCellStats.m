function [cell_center,cell_evecs,cells_radii] = writeRegisteredCellStats(registered_cells,stat_path,expname,whiskername)
%
% writeRegisteredCellStats takes registered cell coordinates as input
% and write an xls file with the following information.
% Number of cells, Volume of the convex hull occupied by the cells,
% distance from centroid, centroid coordinates, principal component extent
% and directions (in degrees)
%
    % Open and write the header of the xls file
    fileID1 = fopen(strcat(stat_path,expname,'_',whiskername,'_stat_registered_cells.csv'),'w');
    fprintf(fileID1,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
    'Exp Name','Whisker Name','CellCount','Slab Volume in Cubic Microns','Distance from FN centroid in Microns',...
    'Centroid Y (Rostral-Caudal axis)','Centroid X (Lateral-Medial axis)','Centroid Z (Ventral-Dorsal axis)',...
    'Eigen Vector 1 Y component','Eigen Vector 1 X component','Eigen Vector 1 Z component',...
    'Eigen Vector 2 Y component','Eigen Vector 2 X component','Eigen Vector 2 Z component',...
    'Eigen Vector 3 Y component','Eigen Vector 3 X component','Eigen Vector 3 Z component',...
    'Eigen Vector 1 Extent in Microns','Eigen Vector 1 angle w.r.t Y in degrees', 'Eigen Vector 1 angle w.r.t X in degrees', 'Eigen Vector 1 elevation w.r.t XY plane in degrees',...
    'Eigen Vector 2 Extent in Microns','Eigen Vector 2 angle w.r.t Y in degrees', 'Eigen Vector 2 angle w.r.t X in degrees', 'Eigen Vector 2 elevation w.r.t XY plane in degrees',...
    'Eigen Vector 3 Extent in Microns','Eigen Vector 3 angle w.r.t Y in degrees', 'Eigen Vector 3 angle w.r.t X in degrees', 'Eigen Vector 3 elevation w.r.t XY plane in degrees');

    % Get cell count and centroid
    i = 1;
    raw_reg(i).landmarks = registered_cells;
    raw_reg(i).cellcount = size(raw_reg(i).landmarks,1);
    raw_reg(i).centroid = mean(raw_reg(i).landmarks,1);
    raw_reg(i).centroidX = raw_reg(i).centroid(1);
    raw_reg(i).centroidY = raw_reg(i).centroid(2);
    raw_reg(i).centroidZ = raw_reg(i).centroid(3);
        
    cell_center = raw_reg(i).centroid;
    % compute distance from center
    raw_reg(i).distFromFNCenter = pdist2(raw_reg(i).centroid, [0,0,0] );
    
    % find pca of cells
    [raw_reg(i).coeff,raw_reg(i).score,raw_reg(i).latent,raw_reg(i).tsquared,raw_reg(i).explained,raw_reg(i).mu] = pca(raw_reg(i).landmarks);
    cell_evecs = raw_reg(i).coeff;
    
    % The desired Axes directions of Average Row slab are as follows
    % PCA1 : Along +ve Y axis
    % PCA2: Along +ve Z axis
    % PCA3: Along +ve X axis
    % Therefore negate if in the oppsote direction to canonical directions
    % mentioned above
    if(cell_evecs(1,1)<0)
        cell_evecs(:,1) = -cell_evecs(:,1);
    end
    if(cell_evecs(3,2)<0)
        cell_evecs(:,2) = -cell_evecs(:,2);
    end
    if(cell_evecs(2,3)<0)
        cell_evecs(:,3) = -cell_evecs(:,3);
    end
        
    pca_extent_reg = range(raw_reg(i).score);
    raw_reg(i).pca_extentX = pca_extent_reg(1);
    raw_reg(i).pca_extentY = pca_extent_reg(2);
    raw_reg(i).pca_extentZ = pca_extent_reg(3);
    
    cells_radii = [pca_extent_reg(1)/2 pca_extent_reg(2)/2 pca_extent_reg(3)/2];
    
    % convert the direction cosines into angles
    [raw_reg(i).theta1, raw_reg(i).z1, raw_reg(i).rho1] = cart2sph(cell_evecs(1,1),cell_evecs(2,1),cell_evecs(3,1));
    [raw_reg(i).theta2, raw_reg(i).z2, raw_reg(i).rho2] = cart2sph(cell_evecs(1,2),cell_evecs(2,2),cell_evecs(3,2));
    [raw_reg(i).theta3, raw_reg(i).z3, raw_reg(i).rho3] = cart2sph(cell_evecs(1,3),cell_evecs(2,3),cell_evecs(3,3));
    
    % error check principal axes angle to stay within 90 degrees
    raw_reg(i).theta1 = abs(rad2deg(raw_reg(i).theta1));
    raw_reg(i).theta2 = abs(rad2deg(raw_reg(i).theta2));
    raw_reg(i).theta3 = abs(rad2deg(raw_reg(i).theta3));
    
    raw_reg(i).z1 = abs(rad2deg(raw_reg(i).z1));
    raw_reg(i).z2 = abs(rad2deg(raw_reg(i).z2));
    raw_reg(i).z3 = abs(rad2deg(raw_reg(i).z3));
        
    if( raw_reg(i).theta1 > 90)
        raw_reg(i).theta1 = 180-raw_reg(i).theta1;
    end
    if( raw_reg(i).theta2 > 90)
        raw_reg(i).theta2 = 180-raw_reg(i).theta2;
    end
    if( raw_reg(i).theta3 > 90)
        raw_reg(i).theta3 = 180-raw_reg(i).theta3;
    end
    
    if( raw_reg(i).z1 > 90)
        raw_reg(i).z1 = 180-raw_reg(i).z1;
    end
    if( raw_reg(i).z2 > 90)
        raw_reg(i).z2 = 180-raw_reg(i).z2;
    end
    if( raw_reg(i).z3 > 90)
        raw_reg(i).z3 = 180-raw_reg(i).z3;
    end
    
    % generate convex hull around the cells and find the volume enclosed by
    % it
    raw_reg(i).dt = delaunayTriangulation(raw_reg(i).landmarks(:,1),raw_reg(i).landmarks(:,2),raw_reg(i).landmarks(:,3));
    [raw_reg(i).vertices, raw_reg(i).volume] = convexHull(raw_reg(i).dt);
    
    % write all the paremeters into the xls sheet
    fprintf(fileID1,'%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',...
        expname,whiskername,...
        raw_reg(i).cellcount,raw_reg(i).volume,raw_reg(i).distFromFNCenter,...
        raw_reg(i).centroid(1),raw_reg(i).centroid(2),raw_reg(i).centroid(3),...
        raw_reg(i).coeff(1,1),raw_reg(i).coeff(2,1),raw_reg(i).coeff(3,1),...
        raw_reg(i).coeff(1,2),raw_reg(i).coeff(2,2),raw_reg(i).coeff(3,2),...
        raw_reg(i).coeff(1,3),raw_reg(i).coeff(2,3),raw_reg(i).coeff(3,3),...
        raw_reg(i).pca_extentX,raw_reg(i).theta1,90-raw_reg(i).theta1,raw_reg(i).z1,...
        raw_reg(i).pca_extentY,raw_reg(i).theta2,90-raw_reg(i).theta2,raw_reg(i).z2,...
        raw_reg(i).pca_extentZ,raw_reg(i).theta3,90-raw_reg(i).theta3,raw_reg(i).z3);
    
    fclose(fileID1);
end

