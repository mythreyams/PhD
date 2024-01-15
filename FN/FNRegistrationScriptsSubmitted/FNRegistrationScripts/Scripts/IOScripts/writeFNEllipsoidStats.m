function writeFNEllipsoidStats(scaled_FN_surface,stat_path,expname,center,evecs,radii)
% writeFNEllipsoidStats takes registered FN surface coordinates as input
% and write an csv file with the following information.
% Centroid coordinates, principal component axes, radii, PCA 
% directions (in degrees) and Volume of FN ellipsoid
%
    volume = 4*pi/3 * radii(1) * radii(2) * radii(3); 
    
    % write the stats into a spreadsheet
    [theta1,z1,r1]  = cart2sph(evecs(1,1),evecs(2,1),evecs(3,1));
    [theta2,z2,r2]  = cart2sph(evecs(1,2),evecs(2,2),evecs(3,2));
    [theta3,z3,r3]  = cart2sph(evecs(1,3),evecs(2,3),evecs(3,3));

    theta1 = abs(rad2deg(theta1));
    theta2 = abs(rad2deg(theta2));
    theta3 = abs(rad2deg(theta3));
    
    z1 = abs(rad2deg(z1));
    z2 = abs(rad2deg(z2));
    z3 = abs(rad2deg(z3));
        
    if( theta1 > 90)
        theta1 = 180 - theta1;
    end
    if( theta2 > 90)
        theta2 = 180 - theta2;
    end
    if( theta3 > 90)
        theta3 = 180 - theta3;
    end
    
    if( z1 > 90)
        z1 = 180 - z1;
    end
    if( z2 > 90)
        z2 = 180 - z2;
    end
    if( z3 > 90)
        z3 = 180 - z3;
    end
    
    fileID2 = fopen(strcat(stat_path,expname,'_', 'stats_FN_Ellipsoid.csv'),'w');
    
    fprintf(fileID2,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
                'Exp Name','Volume of FN Ellipsoid in micron cube','Center X (L-M)','Center Y (R-C)','Center Z (V-D)',...
                'Eigen Vector 1 X','Eigen Vector 1 Y ','Eigen Vector 1 Z',...
                'Eigen Vector 2 X','Eigen Vector 2 Y ','Eigen Vector 2 Z',...
                'Eigen Vector 3 X','Eigen Vector 3 Y ','Eigen Vector 3 Z',...
                'Radius 1 in Microns', 'Radius 2 in Microns', 'Radius 3 in Microns',...
                'Eigen Vector 1 angle w.r.t LM(X) axis in degrees','Eigen Vector 1 angle w.r.t RC(Y) axis in degrees','Eigen Vector 1 elevation w.r.t XY plane in degrees',...
                'Eigen Vector 2 angle w.r.t LM(X) axis in degrees','Eigen Vector 2 angle w.r.t RC(Y) axis in degrees','Eigen Vector 2 elevation w.r.t XY plane in degrees',...
                'Eigen Vector 3 angle w.r.t LM(X) axis in degrees','Eigen Vector 3 angle w.r.t RC(Y) axis in degrees','Eigen Vector 3 elevation w.r.t XY plane in degrees'...
                );

    fprintf(fileID2,'%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',...
              expname,volume,center(1),center(2),center(3),...
              evecs(1,1),evecs(2,1),evecs(3,1),...
              evecs(1,2),evecs(2,2),evecs(3,2),...
              evecs(1,3),evecs(2,3),evecs(3,3),...
              radii(1),radii(2),radii(3),...
              theta1,90-theta1,z1,...
              theta2,90-theta2,z2,...
              theta3,90-theta3,z3...
              );        
            
    fclose(fileID2);
    %writeCSV([center' evecs(:,1)' evecs(:,2)' evecs(:,3)' radii' cartesianToPolar(evecs(:,1)') cartesianToPolar(evecs(:,2)') cartesianToPolar(evecs(:,3)') volume ] ,'stats.csv')
   
end