% neighborhood_prob_2.m
%
% Purpose: Compute neighborhood probability (NP) values for an input 2D
% binary array
% 
% Author(s): Jon Seibert
% Last updated: 25 January 2023
% 
% Inputs: MxN Binary array,  MxN lon array, MxN lat array, 
%   Radius of Influence (ROI) in meters [scalar]
% Outputs: MxN NP array
%
% Process: 
% 1) ROI checking:
% For each point no more than x gridpoints away:
%   a) Retrieve lat/lon coords. Use Haversine distance to measure distance D between centerpoint and target.
%   b) If D <= ROI, add binary value to list
% 2) Then compute mean of all values within radius
%
% Dependencies: haversine_distance.m
%%

function output = neighborhood_prob_2 (binary,lon,lat,roi)

    % Determine dimensions, set up container variables
    dims = size(binary);
    output = zeros(size(binary));
    
    % Compute gridpoint "radius" in which to search by increasing gridpoint
    % radius from 1 until haversine distance is greater than the input ROI,
    % then padding
    rad = 1;
    for idx = 1:dims(1)
        d = haversine_distance(lon(1,1),lat(1,1),lon(idx,1),lat(idx,1));
        if(d > roi)
            rad = min(idx + 3,dims(1));
            break;
        end
    end
    
    % Determine progress points for "progress bar"
    q1 = dims(1)/4;
    q2 = q1*2;
    q3 = q1*3;
    
    % For each point, compare to each other point in "radius". If D <= ROI, add value
    % to sum, increment counter, take mean at end of inner loop.
    tic;
    
    for i_c = 1:dims(1) % Centerpoint coordinates
        for j_c = 1:dims(2)
            
            sum = 0;
            n = 0;
            il_lim = max(i_c - rad,1); iu_lim = min(i_c + rad,dims(1)); % Define search box radius
            jl_lim = max(j_c - rad,1); ju_lim = min(j_c + rad,dims(2));
            
            for i_t = il_lim:iu_lim % Target point coordinates
                for j_t = jl_lim:ju_lim
                    d = haversine_distance(lon(i_c,j_c),lat(i_c,j_c),lon(i_t,j_t),lat(i_t,j_t));
                    if(d <= roi)
                        n = n+1;
                        sum = sum + binary(i_t,j_t);
                    end
                end
            end
            
            output(i_c,j_c) = sum/n; % Average neighborhood probabilities
        end

    
        % ETA output
        if(i_c == dims(1))
            disp("Done.")
        elseif(i_c > q3 && i_c <= q3 + 1)
            disp("Progress: 75%%")
        elseif(i_c > q2 && i_c <= q2 + 1)
            disp("Progress: 50%%")
        elseif(i_c > q1 && i_c <= q1 + 1)
            disp("Progress: 25%%")
        elseif(i_c == 1)
            eta = toc*dims(1)/60;
            fprintf("Estimated time to completion: %.1f minutes", eta)
        end
    
    end

end