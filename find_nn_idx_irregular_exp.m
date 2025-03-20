% find_nn_idx_irregular_exp.m
%
% Purpose: Local function to find 4 nearest neighbor points to the target 
%          lat/lon in an irregularly spaced grid
%
% Author(s): Jon Seibert
% Last updated: 9 June 2023
% 
% Parameters: target_lon, target_lat (point to find 4 nns of)
%             lon_grid, lat_grid (coordinate grids to search)
%             edge_buffer (radius to search from theoretical min)
%
% Outputs: List of coordinates for 4 NN points to target
% Dependencies: -
%
% Basic logic: find closest points along one horizontal and one vertical
% edge, search around intersection of the two in a radius defined by the
% edge_buffer parameter to find true closest point, determine other 3
% nearest neighbors from that
%
% Notes: 
%
%   - Returns index values, not true lat/lon
%   - Designed for a specific orientation of WRF coordinate grids

% !!!! THIS IS THE VERSION IN OPERATIONAL USE !!!!

%%
function [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE] = find_nn_idx_irregular_exp(target_lon,target_lat,lon_grid,lat_grid,edge_buffer)
    
    xNW = 0;
    xNE = 0;
    xSW = 0;
    xSE = 0;
    yNW = 0;
    yNE = 0;
    ySW = 0;
    ySE = 0;    

    if(isnan(target_lon) || isnan(target_lat)) % Immediately cancel and return 0s if inputs are NaN
        return;
    end
    
    % Get dimensions
    dimensions = size(lon_grid);
    xlen = dimensions(2);
    ylen = dimensions(1);

    % Search along edges to find the closest edge points to the target
    % point, then take the (buffer)x(buffer) square around their intersection. 
    % (If more than one closest point, pick smallest idx to base block 
    % around, ensure doesn't overflow boundaries!)
    min_x_array = [0,0];
    y_edge_array = [1 ylen];
    for y_edge_idx = 1:2
        prev_dist = 1e10;
        min_dist = 1e10;
        y = y_edge_array(y_edge_idx);
        for x = 1:xlen
            dist = sqrt((target_lon - lon_grid(y,x))^2 + (target_lat - lat_grid(y,x))^2);
            if(dist > prev_dist)
                break;
            elseif(dist < min_dist)
            %if(dist < min_dist)
                min_dist = dist;
                min_x_array(y_edge_idx) = x;
            end
            prev_dist = dist;
        end
    end
    
    min_y_array = [0,0];
    x_edge_array = [1 xlen];
    for x_edge_idx = 1:2
        prev_dist = 1e10;
        min_dist = 1e10;
        x = x_edge_array(x_edge_idx);
        for y = 1:ylen
            dist = sqrt((target_lon - lon_grid(y,x))^2 + (target_lat - lat_grid(y,x))^2);
            if(dist > prev_dist)
                break;
            elseif(dist < min_dist)
                min_dist = dist;
                min_y_array(x_edge_idx) = y;
            end
            prev_dist = dist;
        end
    end
    
    min_x = round(mean(min_x_array));
    min_y = round(mean(min_y_array));
            
    x_min = max(1, min_x-edge_buffer); x_max = min(xlen, min_x+edge_buffer);
    y_min = max(1, min_y-edge_buffer); y_max = min(ylen, min_y+edge_buffer);
    
    
    % Brute force search remaining area
    min_dist = 1e10;
    min_x = 0;
    min_y = 0;
    
    for x = x_min:x_max
        for y = y_min:y_max
            dist = sqrt((target_lon - lon_grid(y,x))^2 + (target_lat - lat_grid(y,x))^2);
            if(dist < min_dist)
                min_dist = dist;
                min_x = x;
                min_y = y;
            end
        end
    end
    
    min_dist_lat = lat_grid(min_y,min_x);
    min_dist_lon = lon_grid(min_y,min_x);
    
    %%
    % Determine other 3 closest points based on position of closest point
    % relative to target. In the event of a line or total match, duplicate
    % values are used. For WRF: 
    % latitude DECREASES significantly as Y (column) increases
    % longitude INCREASES significantly as X (row) increases
    % Therefore, Y is N -> S, X is W -> E
    % To find indices of other corners, add 1 to X to move East, vv, add 1
    % to Y to move S, vv
    if(min_dist_lat > target_lat) % Mindist lat is North of target
        if(min_dist_lon > target_lon) % Mindist lon is East of target
            % Mindist is NE of target
            xNW = min_x - 1; % Longitude increases moving East, with increasing index
            xNE = min_x;
            xSW = min_x - 1;
            xSE = min_x;
            yNW = min_y;
            yNE = min_y;
            ySW = min_y + 1; % Latitude increases moving North, with decreasing index
            ySE = min_y + 1;
        elseif(min_dist_lon < target_lon) % Mindist lon is West of target
            % Mindist is NW of target
            xNW = min_x;
            xNE = min_x + 1;
            xSW = min_x;
            xSE = min_x + 1;
            yNW = min_y;
            yNE = min_y;
            ySW = min_y + 1;
            ySE = min_y + 1;
        else % Exact lon match
            % Mindist is on lon line with target, to North
            xNW = min_x;
            xNE = min_x;
            xSW = min_x;
            xSE = min_x;
            yNW = min_y;
            yNE = min_y;
            ySW = min_y + 1;
            ySE = min_y + 1;
        end
    elseif(min_dist_lat < target_lat) % Mindist lat is South of target
        if(min_dist_lon > target_lon) % Mindist lon is East of target
            % Mindist is SE of target
            xNW = min_x - 1;
            xNE = min_x;
            xSW = min_x - 1;
            xSE = min_x;
            yNW = min_y - 1;
            yNE = min_y - 1;
            ySW = min_y;
            ySE = min_y;
        elseif(min_dist_lon < target_lon) % Mindist lon is West of target
            % Mindist is SW of target
            xNW = min_x;
            xNE = min_x + 1;
            xSW = min_x;
            xSE = min_x + 1;
            yNW = min_y - 1;
            yNE = min_y - 1;
            ySW = min_y;
            ySE = min_y;
        else % Exact lon match
            % Mindist is on lon line with target, to South
            xNW = min_x;
            xNE = min_x;
            xSW = min_x;
            xSE = min_x;
            yNW = min_y - 1;
            yNE = min_y - 1;
            ySW = min_y;
            ySE = min_y;
        end
    else % Exact lat match
        if(min_dist_lon > target_lon) % Mindist lon is East of target
            % Mindist is on lat line with target, to East
            xNW = min_x - 1;
            xNE = min_x;
            xSW = min_x - 1;
            xSE = min_x;
            yNW = min_y;
            yNE = min_y;
            ySW = min_y;
            ySE = min_y;
        elseif(min_dist_lon < target_lon) % Mindist lon is West of target
            % Mindist is on lat line with target, to West
            xNW = min_x;
            xNE = min_x + 1;
            xSW = min_x;
            xSE = min_x + 1;
            yNW = min_y;
            yNE = min_y;
            ySW = min_y;
            ySE = min_y;
        else % Exact lon match
            % Mindist is exact match with target
            xNW = min_x;
            xNE = min_x;
            xSW = min_x;
            xSE = min_x;
            yNW = min_y;
            yNE = min_y;
            ySW = min_y;
            ySE = min_y;
        end
    end
    %%
end