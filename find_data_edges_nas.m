function[first_data_lat,last_data_lat,first_data_lon,last_data_lon] = find_data_edges_nas(lat,lon,refl_slice)


mid_lon = round(lon/2); % Find center of one horizontal and one vertical side, to avoid hitting an entire line of NaN
q1_lon = round(lon/4);
q3_lon = round(lon*3/4);
mid_lat = round(lat/2);
q1_lat = round(lat/4);
q3_lat = round(lat*3/4);

% Use furthest out of three points chosen to avoid accidentally hitting one
% dip in the data and cutting off too much

% Find N data border
for lat_idx = 1:lat
    if(~isnan(refl_slice(lat_idx,q1_lon)))
        lat_a = lat_idx;
        break;
    end
end
for lat_idx = 1:lat
    if(~isnan(refl_slice(lat_idx,mid_lon)))
        lat_b = lat_idx;
        break;
    end
end
for lat_idx = 1:lat
    if(~isnan(refl_slice(lat_idx,q3_lon)))
        lat_c = lat_idx;
        break;
    end
end
first_data_lat = min([lat_a,lat_b,lat_c]);

% Find S data border
for lat_idx = lat:-1:1
    if(~isnan(refl_slice(lat_idx,q1_lon)))
        lat_a = lat_idx;
        break;
    end
end
for lat_idx = lat:-1:1
    if(~isnan(refl_slice(lat_idx,mid_lon)))
        lat_b = lat_idx;
        break;
    end
end
for lat_idx = lat:-1:1
    if(~isnan(refl_slice(lat_idx,q3_lon)))
        lat_c = lat_idx;
        break;
    end
end
last_data_lat = max([lat_a,lat_b,lat_c]);

% Find W data border
for lon_idx = 1:lon
    if(~isnan(refl_slice(q1_lat,lon_idx)))
        lon_a = lon_idx;
        break;
    end
end
for lon_idx = 1:lon
    if(~isnan(refl_slice(mid_lat,lon_idx)))
        lon_b = lon_idx;
        break;
    end
end
for lon_idx = 1:lon
    if(~isnan(refl_slice(q3_lat,lon_idx)))
        lon_c = lon_idx;
        break;
    end
end
first_data_lon = min([lon_a,lon_b,lon_c]);

% Find E data border
for lon_idx = lon:-1:1
    if(~isnan(refl_slice(q1_lat,lon_idx)))
        lon_a = lon_idx;
        break;
    end
end
for lon_idx = lon:-1:1
    if(~isnan(refl_slice(mid_lat,lon_idx)))
        lon_b = lon_idx;
        break;
    end
end
for lon_idx = lon:-1:1
    if(~isnan(refl_slice(q3_lat,lon_idx)))
        lon_c = lon_idx;
        break;
    end
end
last_data_lon = max([lon_a,lon_b,lon_c]);


end