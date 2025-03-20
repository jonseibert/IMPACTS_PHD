function vm = linear_interpolate(x1,x2,xm,v1,v2)
    % Linear interpolation uses the distance between two points to compute
    % the value at the chosen point between them.
    v_dif = v2 - v1;
    x_dif = x2 - x1;
    x_proportion = (xm - x1)/x_dif;
    vm = v1 + x_proportion*v_dif;
end