function [points] = readPtsSpatialGraphAmira(filename)
    %filename = 'MG04_SpatialGraphSet.merged.am'
    % Open the text file.
    %if (~(strcmp(filename(end-3:end), '.hoc')))
    %filename = [filename '.am'];
    %end

    fileID = fopen(filename,'r');

    tline = fgetl(fileID); 

        %if ~isempty(idx)
    i = 0;

    points = [];
    idxFound = 0;
    tag = [];
    while ischar(tline)

        if (~isempty(regexp(tline,'EdgePointCoordinates')))
            tag = tline(1,size(tline,2)-1:size(tline,2));
        end
        if strcmp(tline,tag)
            idxFound = 1;
            tline = fgetl(fileID);
        end
        
        if idxFound==1

            if strcmp(tline,'')
                break;
            end

            i = i + 1;
            C = textscan(tline,'%f %f %f','delimiter',' ');
            points(i,:) = cell2mat(C);
        end

        tline = fgetl(fileID); 

    end

    fclose(fileID);
end