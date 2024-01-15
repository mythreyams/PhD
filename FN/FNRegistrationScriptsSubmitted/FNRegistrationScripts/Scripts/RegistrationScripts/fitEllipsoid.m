function [center, radii, evecs, FN_surface_points] = fitEllipsoid(FN_contours)
% Fit an ellipsoid to the FN contours using llsq method and return the
% center, radii and proncipal axes. Also get thesurface points of the fit
% FN ellipsoid.
%

%%%%%%%%%%%%%%%%%%%%%% Caution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In order for this script to produce results consistent with the paper
% the following directions need to be maintained for inputs. In cases where
% the input directions are not as expected the output can be converted
% appropriately (see code below)

% X Axis : Lateral - Medial
% Y Axis : Rostral - Caudal
% Z Axis : Ventral - Dorsal
%%%%%%%%%%%%%%%%%%%%%% Caution End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           
    % fit an ellipsoid to given points and get the center, radius (orientation and magnitude)
    [ center, radii, evecs, v, chi2 ] = ellipsoid_fit( FN_contours);

    % The desired Axes directions of Average FN are as follows
    % PCA1 : Along +ve Y axis
    % PCA2: Along +ve X axis
    % PCA3: Along +ve Z axis
    % Therefore negate if in the oppsote direction to canonical directions
    % mentioned above
    if(evecs(2,1)<0)
        evecs(:,1) = -evecs(:,1);
    end
    if(evecs(1,2)<0)
        evecs(:,2) = -evecs(:,2);
    end
    if(evecs(3,3)<0)
        evecs(:,3) = -evecs(:,3);
    end

    % create a 3D mesh points to evaluate the ellipsoid function
    mind = min( [ FN_contours(:,1) FN_contours(:,2) FN_contours(:,3) ] )-1000;
    maxd = max( [ FN_contours(:,1) FN_contours(:,2) FN_contours(:,3) ] )+1000;
    nsteps = 100;
    step = ( maxd - mind ) / nsteps;
    [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );

    % Evaluate the ellipsoid function at the grid points generated above
    Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
              2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
              2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
    p = ( isosurface( x, y, z, Ellipsoid, -v(10) ) );

    FN_surface_points = p.vertices;
    
end

function [ center, radii, evecs, v, chi2 ] = ellipsoid_fit( X )
%
% ellipsoid_fit fits an ellipsoid to the given 3D points using llsq method
% and returns the coefficients of the ellipsoid equation in quadric surface
% form along with center, radii and cosine direction vectors.
%
    x = X( :, 1 );
    y = X( :, 2 );
    z = X( :, 3 );

% fit ellipsoid in the form Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx +
% 2Hy + 2Iz + J = 0 and A + B + C = 3 constraint removing one extra
% parameter

    D = [ x .* x + y .* y - 2 * z .* z, ...
        x .* x + z .* z - 2 * y .* y, ...
        2 * x .* y, ...
        2 * x .* z, ...
        2 * y .* z, ...
        2 * x, ...
        2 * y, ...
        2 * z, ...
        1 + 0 * x ];  % ndatapoints x 9 ellipsoid parameters

% solve the normal system of equations
d2 = x .* x + y .* y + z .* z; % the RHS of the llsq problem (y's)
u = ( D' * D ) \ ( D' * d2 );  % solution to the normal equations

% find the ellipsoid parameters
% convert back to the conventional algebraic form

    v(1) = u(1) +     u(2) - 1;
    v(2) = u(1) - 2 * u(2) - 1;
    v(3) = u(2) - 2 * u(1) - 1;
    v( 4 : 10 ) = u( 3 : 9 );

v = v';

% form the algebraic form of the ellipsoid
A = [ v(1) v(4) v(5) v(7); ...
      v(4) v(2) v(6) v(8); ...
      v(5) v(6) v(3) v(9); ...
      v(7) v(8) v(9) v(10) ];
  
% find the center of the ellipsoid
center = -A( 1:3, 1:3 ) \ v( 7:9 );
% form the corresponding translation matrix
T = eye( 4 );
T( 4, 1:3 ) = center';
% translate to the center
R = T * A * T';
% solve the eigenproblem
[ evecs, evals ] = eig( R( 1:3, 1:3 ) / -R( 4, 4 ) );
radii = sqrt( 1 ./ diag( abs( evals ) ) );
sgns = sign( diag( evals ) );

% calculate difference of the fitted points from the actual data normalized by the conic radii
d = [ x - center(1), y - center(2), z - center(3) ]; % shift data to origin
d = d * evecs; % rotate to cardinal axes of the conic;
d = [ d(:,1) / radii(1), d(:,2) / radii(2), d(:,3) / radii(3) ]; % normalize to the conic radii
chi2 = sum( abs( 1 - sum( d.^2 .* repmat( sgns', size( d, 1 ), 1 ), 2 ) ) );

if abs( v(end) ) > 1e-6
    v = -v / v(end); % normalize to the more conventional form with constant term = -1
else
    v = -sign( v(end) ) * v;
end

end

