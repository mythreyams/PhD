% Transform coordinates from Amira Reference Frame to Canonical Reference Frame
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
function [Canonical_coordinates] = CovertAmiraToCanonicalCoordinates(Amira_coordinates)
    Canonical_coordinates = [Amira_coordinates(:,1) -Amira_coordinates(:,3) -Amira_coordinates(:,2)];
end