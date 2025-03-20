function c = reflmap()
% reflmap.m
% Purpose: Custom radar reflectivity colormap
% Author(s): Jon Seibert
% Last updated: 18 Nov 2021
% Parameters: none
% Output: c (3x21 colormap)

    % Cutting out one of the greys makes the colors line up to the numbers properly on caxis([-30 75])
    c = [255 255 255; % -35 (and lower) < x < -30 (White)
            %230 230 230;  % -30:-25 (Light Grey)  
            230 230 230;  % -25:-20 (Light Grey)
            230 230 230;  % -20:-15 (Light Grey)
            230 230 230;  % -15:-10 (Light Grey)
            230 230 230;  % -10:-5 (Light Grey)
            230 230 230;  % -5 <= x < 0 dbz (Light Grey)
            0 255 255;	  % 0 <= x < 5 dbz (Cyan)
            75 200 255;   % 5 <= x < 10 dbz (Light Blue)
            30 144 255;   % 10 <= x < 15 dbz (Blue)
            0 0 255;	  % 15 <= x < 20 dbz (Indigo)
            0 255 0; 	  % 20 <= x < 25 dbz (Lime)
            0 205 0; 	  % 25 <= x < 30 dbz (Green)
            0 139 0; 	  % 30 <= x < 35 dbz (Dark Green)
            255 255 0; 	  % 35 <= x < 40 dbz (Yellow)
            255 215 0;	  % 40 <= x < 45 dbz (Gold)
            255 127 0; 	  % 45 <= x < 50 dbz (Orange)
            255 0 0;	  % 50 <= x < 55 dbz (Scarlet)
            205 0 0; 	  % 55 <= x < 60 dbz (Red)
            139 0 0; 	  % 60 <= x < 65 dbz (Dark Red)
            255 0 255;	  % 65 <= x < 70 dbz (Fuschia)
            145 44 238];  % 70 <= x < 75 dbz (Purple)

    c = c/255;

end