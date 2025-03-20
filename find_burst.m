
max_q = 0;
max_q_lon = 0;
max_q_lat = 0;
max_q_x = 0;
max_q_y = 0;
for y = 1:ylen_wrf
    for x = 1:xlen_wrf
        if(lon_wrf(y,x) > -75 && lon_wrf(y,x) < -72 && lat_wrf(y,x) > 42 && lat_wrf(y,x) < 44)
            if(data_a(y,x) > max_q)
                max_q = data_a(y,x);
                max_q_lon = lon_wrf(y,x);
                max_q_lat = lat_wrf(y,x);
                max_q_x = x;
                max_q_y = y;
            end
        end
    end
end