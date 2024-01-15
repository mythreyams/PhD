function [centered,rotated,scaled,reordered_pcas] = rigidTransformation(canonical_points, center, eig, rf_rad,rad, translate_back,scale,inverse,reorder)
% RigidTransformation performs the rigid transformation of input points.
%
% Input:
%   canonical_points: Input points to be transformed
%
%   center: Centroid of input points.
%
%   eig: Principal component directions of input points.
%
%   rf_rad: Radii of Average FN reference.
%
%   rad: Radii of the FN to be transformed
%   rf_rad/rad is used to scaling the input points to Avg Reference Frame.
%
%   translate_back: Bool flag saying whehter or not to translate the
%   transformed points back to given center.
%
%   scale: bool flag for scaling
%
%   inverse: bool flag for forward or reverese rotation
%
%   reorder: bool flag stating whether to reorder the principal components
%   in order to reorder the PCAs in order of which is closest to rf PCAs.
%   This helps align axes that are closest to reference axes with each
%   other.
%
% Output:
%   centered: Input points centered to origin (0,0,0).
%   rotated:  Input points centered and rotated as per the given axes
%   scaled: Input points centered, rotated and scaled to match the
%   reference radii.

% Check if reordering the PCAs is necessary
if(reorder)
[pcas,ind] = reorder_pca(eig);
else
pcas = eig;
ind = [1,2,3];
end

% move the points to be centered at origin
centered(:,1) = canonical_points(:,1) - center(1);
centered(:,2) = canonical_points(:,2) - center(2);
centered(:,3) = canonical_points(:,3) - center(3);

% rotate
if(inverse)
    rotated = pcas' * centered';
else
    rotated = pcas * centered';
end

rotated =  rotated';

% scale
if(scale)
    scaled(:,1) = rotated(:,1) * (rf_rad(1)/rad(ind(1)));  
    scaled(:,2) = rotated(:,2) * (rf_rad(2)/rad(ind(2))) ; 
    scaled(:,3) = rotated(:,3) * (rf_rad(3)/rad(ind(3))) ;
else
    scaled(:,1) = rotated(:,1);  
    scaled(:,2) = rotated(:,2); 
    scaled(:,3) = rotated(:,3);
end

% tranlate back to given center
if(translate_back == 1)
    scaled(:,1) =  scaled(:,1) +  center(1);
    scaled(:,2) =  scaled(:,2) +  center(2);
    scaled(:,3) =  scaled(:,3) +  center(3);    
end

reordered_pcas = pcas;

end

function [v,ind] = reorder_pca(u)
%   reorder_pca bool flag stating whether to reorder the principal components
%   in order of which is closest to reference PCAs.
%   This helps align axes that are closest to reference axes with each
%   other.
%
v = u;

[bla, orders(:,1)] = sort(abs(u(:,1)),'descend');
[bla, orders(:,2)] = sort(abs(u(:,2)),'descend');
[bla, orders(:,3)] = sort(abs(u(:,3)),'descend');

[ind1] = find(orders(:,1) == 3);
[ind2] = find(orders(:,2) == 3);
[ind3] = find(orders(:,3) == 3);

[val,indz] =  min([ind1 ind2 ind3]);

% find closest to z axis as the first guy
%[val1,ind1] = find(orders(1,:) == [3]);

if(u(3,indz) > 0)
    v(:,1) = u(:,indz);
else
    v(:,1) = -u(:,indz);
end

% find closest to x axis as the second guy skipping the first guy from
% consideration
if (indz == 1)
    [ind2] = find(orders(:,2) == 1);
    [ind3] = find(orders(:,3) == 1);
    [val,indx] =  min([10 ind2 ind3]);
elseif (indz == 2)
    [ind1] = find(orders(:,1) == 1);
    [ind3] = find(orders(:,3) == 1);
    [val,indx] =  min([ind1 10 ind3]);
else
    [ind1] = find(orders(:,1) == 1);
    [ind2] = find(orders(:,2) == 1);
    [val,indx] =  min([ind1 ind2 10]);   
end

if(u(1,indx) > 0)
    v(:,2) = u(:,indx);
else
    v(:,2) = -u(:,indx);
end

% get the third axis
for i = 1:3
    if ~((i == indz)|| (i == indx))
        indy = i;
    end
end

if(u(2,indy) > 0)
    v(:,3) = u(:,indy);
else
    v(:,3) = -u(:,indy);
end

% set the new order
ind = [indz indx indy];
end
