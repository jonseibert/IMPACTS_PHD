function xm = value_interpolate(x1,x2,v1,v2,vm)
    % Linear interpolation uses the distance between two points to compute
    % the position of a desired value that lies between the two.
    v_dif = v2 - v1;
    x_dif = x2 - x1;
    v_proportion = (vm - v1)/v_dif;
    %x_proportion = (xm - x1)/x_dif;
    %vm = v1 + x_proportion*v_dif;
    xm = x1 + v_proportion*x_dif;
end