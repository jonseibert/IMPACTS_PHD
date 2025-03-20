function c = tmap()
% tmap.m
% Purpose: Custom colormap for temperature contour plots. Highlights
% freezing line.
% Author(s): Jon Seibert
% Last updated: 13 Jan 2023
% Parameters: none
% Output: c (3x9 colormap)

    % Designed to fit on a 5-degree interval for [-20 20] Celsius
    c =    [150 0 150;
            0 40 255; % -20: Darkish blue
            60 120 255;
            0 255 255;
            0 200 150;
            255 0 0;	  % 0 C (Scarlet)
            0 200 0;
            100 210 50;
            150 220 50;
            175 230 50;
            200 200 25]; % 20:
    c = c/255;

end