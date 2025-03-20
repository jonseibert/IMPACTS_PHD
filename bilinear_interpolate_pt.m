function vM = bilinear_interpolate_pt(lat_target,lon_target,lat_src_grid,lon_src_grid,data_src,nn_points)
% ------------------------------------------------------------------------
% bilinear_interpolate_irregular_to_grid.m
% 
% Purpose: Interpolate four irregularly spaced grid points using provided indices.
% 
% Author(s): Jon Seibert
% Last updated: 24 Oct 2024
% 
% Inputs: 
%   -nn_points [Array]: 8-member array containing the x and y indices of the
%                      target point's 4 Nearest Neighbor points in the source grid
% 
% Outputs: 
%   -vM [Float]: The interpolated value that would exist at the
%               target point within the source grid.
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
    
    vN = linear_interpolate(lon_src_grid(yNW,xNW),lon_src_grid(yNE,xNE),lon_target,data_src(yNW,xNW),data_src(yNE,xNE));
    vS = linear_interpolate(lon_src_grid(ySW,xSW),lon_src_grid(ySE,xSE),lon_target,data_src(ySW,xSW),data_src(ySE,xSE));
    vM = linear_interpolate(lat_src_grid(yNW,xNW),lat_src_grid(ySW,xSW),lat_target,vN,vS);
    
end