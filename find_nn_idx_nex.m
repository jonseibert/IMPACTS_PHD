% Local function to find 4 nearest neighbor points to the "NASA" point in
% the "NEXRAD" lat/lon grid
% Basic logic: loop through each NEXRAD lat value- if target lat passed, nn
% are i and i-1. If target == current, nn are i and i. Do the same for lon.
% Returns index values, not true lat/lon
function [lon1,lon2,lat1,lat2] = find_nn_idx_nex(target_lon,target_lat,lon_values,lat_values)
    lon1 = 0;
    lon2 = 0;
    lat1 = 0;
    lat2 = 0;
    
    for i = 1:length(lon_values)
        if(lon_values(i) >= target_lon)
            if(lon_values(i) == target_lon)
                lon2 = i;
                lon1 = 1;
                break;
            else
                lon2 = i;
                lon1 = i-1;
                break;
            end
        end
    end
    
    for i = 1:length(lat_values)
        if(lat_values(i) >= target_lat)
            if(lat_values(i) == target_lat)
                lat2 = i;
                lat1 = 1;
                break;
            else
                lat2 = i;
                lat1 = i-1;
                break;
            end
        end
    end
    
end