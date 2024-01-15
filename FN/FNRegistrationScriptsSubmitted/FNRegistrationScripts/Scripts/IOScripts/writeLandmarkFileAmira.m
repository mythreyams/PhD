% Write a landmark file which can be read by amira
% Input:
% - points: cloud of points [numPtsx3]
% - outputFilename
function [] = writeLandmarkFileAmira(points,outputFilename,Coordinatedirections)
        
        points = [points(:,1)*Coordinatedirections(1) points(:,2)*Coordinatedirections(2) points(:,3)*Coordinatedirections(3)];
        
        
        fname = [outputFilename '.landmarkAscii'];
        fid = fopen(fname,'w'); 
        
        numPts = size(points,1); 
        
        if size(points,2) ~= 3
           error('Points are not 3D-coordinates!') 
        end
        
        str = sprintf(['# AmiraMesh 3D ASCII 2.0\n\n' ...
            'define Markers %d\n\n' ...
            'Parameters {\n' ...
            '\tNumSets 1,\n' ...
            '\tContentType "LandmarkSet"\n}\n\n' ...
            'Markers { float[3] Coordinates } @1\n\n' ...
            '# Data section follows\n' ...
            '@1'],numPts); 
        
        if fid ~= -1
            fprintf(fid,'%s\n',str);   % no \r
            fclose(fid);
        end
        
        dlmwrite(fname,points,'-append','delimiter',' ','newline','pc'); 
end
