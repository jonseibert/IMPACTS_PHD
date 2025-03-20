function c = tmap()
% reflmap.m
% Purpose: Custom radar reflectivity colormap
% Author(s): Jon Seibert
% Last updated: 18 Nov 2021
% Parameters: none
% Output: c (3x21 colormap)

    % Designed to fit on a 5-degree interval for [-20 20] Celsius
    c =    [10 75 255; %
            60 120 255;
            0 255 255;
            0 200 150;
            255 0 0;	  % 50 <= x < 55 dbz (Scarlet)
            100 200 50;
            150 255 50;
            200 255 50;
            255 255 50];
    c = c/255;

end