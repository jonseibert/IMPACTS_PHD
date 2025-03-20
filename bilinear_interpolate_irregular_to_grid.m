function vM = bilinear_interpolate_irregular_to_grid(y_t,x_t,lon_target,lat_target,data_src,lon_src,lat_src,nn_points)
% ------------------------------------------------------------------------
% bilinear_interpolate_irregular_to_grid.m
% 
% Purpose: Interpolate four irregularly spaced grid points using provided indices.
% 
% Author(s): Jon Seibert
% Last updated: 22 August 2022
% 
% Inputs: 
%   -x_t [Int]: X-index of latlon values in the target grid (grid being
%              interpolated to)
%   -y_t [Int]: Y-index of latlon values in target grid
%   -lon_target [Array]: 2D longitude matrix of target grid
%   -lat_target [Array]: 2D latitude matrix of target grid
%   -data_src [Array]: 2D value matrix of source grid (being interpolated from)
%   -lon_src [Array]: 2D longitude matrix of source grid
%   -lat_src [Array]: 2D latitude matrix of source grid
%   -nn_points [Array]: 8-member array containing the x and y indices of the
%                      target point's 4 Nearest Neighbor points in the source grid
% 
% Outputs: 
%   -vM [Float]: The interpolated value that would exist at the
%               target point within the source grid.
%
% TODO:
%   -
%
% NOTES:
%   - Bilinear interpolation is simply a matter of using three linear
%   interpolations to triangulate the value at the location between the
%   four cornerpoints. However, when moving from an irregular grid to a
%   regular grid, the target is a specific point, which may not be the
%   midpoint.
%   To solve this problem, use proportional distances rather than simple
%   means. It may still be computed in 3 total interpolations- one
%   between each chosen opposing edge (one of the two sets), and one
%   between the two results.
%   The math does check out on this- the x distances are proportional to
%   the y distances as long as we're only interpolating along straight
%   lines (which this code does).
% ------------------------------------------------------------------------
%%    

    % Retrieve cornerpoint coordinate indices
    xNW = nn_points(1);
    xNE = nn_points(2);
    xSW = nn_points(3);
    xSE = nn_points(4);
    yNW = nn_points(5);
    yNE = nn_points(6);
    ySW = nn_points(7);
    ySE = nn_points(8);
    
    % First, compute the location and value of a "midpoint" between each of the corner points on two
    % opposing sides (N and S) that intersects with one coordinate line of
    % the target, by using the proportional x and y distances from each
    % corner point to the target
    % Then, compute the value at the target by using the proportional Y
    % distance between the target and the two midpoints

    % Method 1: Coordinates
    
    % Theoretical
    %xN_scaling = (x_t - xNW)/(xNE - xNW);
    %lonN = x_t;
    %latN = yNW + (yNE - yNW)*xN_scaling;
    %vN = data_src(yNW,xNW) + (data_src(yNE,xNE) - data_src(yNW,xNW))*xN_scaling;

    %xS_scaling = (x_t - xSW)/(xNE - xSW);
    %lonS = x_t;
    %latS = ySW + (ySE - ySW)*xS_scaling;
    %vS = data_src(ySW,xSW) + (data_src(ySE,xSE) - data_src(ySW,xSW))*xS_scaling;
    
    %yM_scaling = (y_t - latS)/(latN - latS);
    %lonM = x_t;
    %latM = y_t;
    %vM = vS + (vS - vN)*yM_scaling;
    
    %
    % Check for NaN lat/lon values
    %ll_val = [lon_src(yNW,xNW),lon_src(yNE,xNE),lon_src(ySW,xSW),lon_src(ySE,xSE),lon_target(y_t,x_t),lat_target(y_t,x_t)];
    %if(any(isnan(ll_val)))
    %    vM = NaN;
    %    return;
    %end
    
    % Perform calculations
    
    % If the difference is 0, the N points are on the same x line.
    if(lon_src(yNE,xNE) - lon_src(yNW,xNW) == 0) % If so, take mean for lat and value.
        latN = (lat_src(yNW,xNW) + lat_src(yNE,xNE))/2;
        vN = (data_src(yNW,xNW) + data_src(yNE,xNE))/2;
    else % Otherwise, compute from proportions
        xN_scaling = (lon_target(y_t,x_t) - lon_src(yNW,xNW))/(lon_src(yNE,xNE) - lon_src(yNW,xNW));
        latN = lat_src(yNW,xNW) + (lat_src(yNE,xNE) - lat_src(yNW,xNW))*xN_scaling;
        vN = data_src(yNW,xNW) + (data_src(yNE,xNE) - data_src(yNW,xNW))*xN_scaling;
    end

    % Repeat this process for the south side.
    if(lon_src(ySE,xSE) - lon_src(ySW,xSW) == 0)
        latS = (lat_src(ySW,xSW) + lat_src(ySE,xSE))/2;
        vS = (data_src(ySW,xSW) + data_src(ySE,xSE))/2;
    else
        xS_scaling = (lon_target(y_t,x_t) - lon_src(ySW,xSW))/(lon_src(ySE,xSE) - lon_src(ySW,xSW));
        latS = lat_src(ySW,xSW) + (lat_src(ySE,xSE) - lat_src(ySW,xSW))*xS_scaling;
        vS = data_src(ySW,xSW) + (data_src(ySE,xSE) - data_src(ySW,xSW))*xS_scaling;
    end
    
    % Then interpolate the north and south points to the middle- if the N
    % and S lat values are the same, it was a perfect match, take the value
    % of any of the 4 points.
    if(latN - latS == 0)
        vM = vN;
    else
        yM_scaling = (lat_target(y_t,x_t) - latS)/(latN - latS);
        vM = vS + (vN - vS)*yM_scaling;
    end
    
end