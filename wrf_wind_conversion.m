% wrf_wind_conversion.m
% Purpose: Calculate partial-3D radial velocity from 3D wind fields relative to a single radar point
%   - For now: KBUF
%   - Elevation accommodated for by adding the elevation of KBUF to all
%     height values
% 
% Author(s): Jon Seibert
% Last updated: 18 January 2022

% TODO:
% - HOW TO PLOT NANs IN A DIFFERENT COLOR (GREY)
% - Generalize to multiple/other choice of radar
%   - Next step: use filename formats for saves/loads and plot saves

%% 0a. General Options

% Specify paths
base_path = pwd;
input_path = "wrf_data/m1";
intermediate_path = "intermediate";
output_path = "output/radial_velocity/wrf";

% File lists
mean_files = string(ls("wrf_data/ncea/wrf*"));
mem_files = string(ls("wrf_data/m1/wrf*"));
ex_files = string(ls("nasa_data/exrad_impacts*"));
radar_data_filename = "NEXRAD_data.mat";

% Standardize output filepaths and names
filepath_format_vr = '%s/Vr_WRF_%s_%s_%d.mat'; % intermediate_path,fit/full,type_label,wrf_file_choice
filepath_format_plot = '%s/Vr_WRF_%s_%s_%s_%s_e%d.png'; % output_path,fit/full/dif,comparison_data,type_label,hour,angle_idx
    % ex. Vr_WRF_full_wrf_001_14_e5.png, Vr_WRF_dif_ex_mean_14_e1.png
%title_format = '%s';

% Specify the list of radar elevation angles
% (0.5-19.5 is WSR-88D range, in increments of 0.5*)
elevation_angles = single(0.5:0.5:19.5); 
radar_pad = 0.5; % Amount of +- padding on the radar elevation angle to determine beam coverage

% Specify maximum range of radar
radar_max_dist = 200*1.852*1000; %200 nautical miles * (1.852 km/nm) * (1000 m/km) [FOR WIND- REFL MAX DIST IS 150KM, FROM GEM]

% Paired list of radar sites and their elevations
% KAPX, KBGM, KTYX, KBUF, KCLE, KDTX, KGRR, KENX, KCCX, KCXX, KPBZ
radar_elevations = [0,0,0,708.7,0,0,0,0,0,0,0];

% Plot dimensions
plot_width = 1600;
plot_height = 900;
cmap = 1; % Which colormap to plot Vr values with (1: jet, 2: redblue)

% Plot font sizes
title_font_size = 18;
label_font_size = 14;
axes_font_size = 12;

% Plot angle cutoff (only plot the first X elevation angles)
angle_cutoff = 16;

%% 0b. Data options

radar_idx = 4; % Radar choice index [4 = KBUF]
wrf_file_choice = 2; % Which WRF file to sample for development
type_choice = 1; % Use member (1) or mean (2)
exrad_file_choice = 1; % Index of EXRAD file that goes with chosen WRF file, timewise
exrad_comparison_filename = 'EXRAD_Vr_angled_1.mat'; % Saved EXRAD Vr file. Assumes variable is named "Vr_angled".

%% 1. Extract WRF 3D wind data

% Determine if member or mean is being used
if(type_choice == 1)
    filename = sprintf("%s/%s",input_path,mem_files(wrf_file_choice));
    type_label = '001';
else
    filename = sprintf("%s/%s",input_path,mean_files(wrf_file_choice));
    type_label = 'mean';
end

% Extract data
longitude = ncread(filename,"XLONG"); % WRF lon mostly increases with x, lat increases with y [(x,y)]
latitude = ncread(filename,"XLAT");
reflectivity = ncread(filename,"REFL_10CM");
u_offset = ncread(filename,"U");
v_offset = ncread(filename,"V");
w_offset = ncread(filename,"W");

% Need these two to compute z: z = (PH+PHB)/g  where g is 9.81 m/s2 
PH = ncread(filename,"PH"); % 'perturbation geopotential'
PHB = ncread(filename,"PHB"); % 'base-state geopotential'
g = 9.81; % m/s^2

% Get string of hour in timestamp
timestamp = convertCharsToStrings(ncread(filename,"Times"));
hour = extractBetween(timestamp,12,13);
hour = hour{1};

% Remake timestamp without "_" in the string, to avoid funkiness
ts_1 = extractBetween(timestamp,1,10);
ts_2 = extractBetween(timestamp,12,19);
timestamp = sprintf("%s-%s",ts_1{1},ts_2{1});

% Specify dimensions
dims = size(reflectivity); 
xlen = dims(1);
ylen = dims(2);
zlen = dims(3);

% Define number of elevation angles on radar
num_angles = length(elevation_angles);

% Elevation above sea level of chosen radar
radar_height = radar_elevations(radar_idx); 

%% 2. Convert offset wind grids to singular central wind grid

% lat,lon,time / x,y,t / row,col,-
% latlon_u has 1 extra entry in x/latitude for grid edge winds
% latlon_v has 1 extra in y/lon

% Set up holding containers
u_full = zeros(size(reflectivity));
v_full = zeros(size(reflectivity));
w_full = zeros(size(reflectivity));
PH_centered = zeros(size(reflectivity));
PHB_centered = zeros(size(reflectivity));

dims_centered = size(reflectivity);

for x = 1:dims_centered(1)
    for y = 1:dims_centered(2)
        u_full(x,y,:) = (u_offset(x,y,:) + u_offset(x+1,y,:))/2;
        v_full(x,y,:) = (v_offset(x,y,:) + v_offset(x,y+1,:))/2;
        for z = 1:dims_centered(3)
            PH_centered(x,y,z) = (PH(x,y,z) + PH(x,y,z+1)/2);
            PHB_centered(x,y,z) = (PHB(x,y,z) + PHB(x,y,z+1)/2);
            w_full(x,y,z) = (w_offset(x,y,z) + w_offset(x,y,z+1))/2;
        end
    end
end

%% 3. Calculate Vr (radial velocity) per elevation angle
% Instead of computing and interpolating Vr, interpolate 3D wind and use it
% to directly compute angle-relative Vr

% Radial velocity: Vr
% Elevation angle: theta-e
% Azimuthal angle: phi
% Ground-relative range from the radar: r

% Vr = u*sin(phi)*cos(theta_e) + v*cos(phi)*cos(theta_e) + w*sin(theta_e)

% Radar azimuth angle defined as clockwise from N
% Math angle defined as counter-clockwise from E
% To convert to math form in radians from radar form in degrees:
% phi_math = -1*(phi - 90)*(pi/180)

% HOW TO COMPUTE WHETHER A POINT IS CLOSE ENOUGH TO THE RADAR BEAM TO BE
% INCLUDED AT A GIVEN ELEVATION ANGLE:
%   - Beam width can be considered an extra 0.5* on either side of the
%   central angle: elevation 0.5 includes 0*-1*
%   - therefore, to determine if a point is within the beam at a specific
%   elevation angle, calculate the height of the top and bottom of the beam
%   and compare to height of point

% z = (PH+PHB)/g
heights_full = (PH_centered+PHB_centered)/g; % m above sea level

% Load radar data
load(sprintf("%s/%s",intermediate_path,radar_data_filename),'lons','lats','ids');

id_radar = ids(radar_idx);
lat_radar = lats(radar_idx); %lat_radar = 42.9489 for KBUF
lon_radar = lons(radar_idx); %lon_radar = 78.7367 for KBUF

Vr_angled = single(zeros(xlen,ylen,num_angles)); % Empty holding variable for radial velocity

% Loop over all points in the 3D wind data for each angle
% For each point: is it in range of the radar? If true, interpolate the 3D
% winds to the radar beam centerpoint and compute Vr if possible (has above
% and below to interpolate from)
% otherwise, NaN
for x = 1:xlen
    for y = 1:ylen
        % Check if point is within horizontal range (ignores height)
        % Note that lat/lon does not vary with height for WRF data
        dist = haversine_distance(longitude(x,y,1),latitude(x,y,1),lon_radar,lat_radar);
        if(dist <= radar_max_dist)

            % Compute azimuth angle
            lat_target = latitude(x,y,1);
            lon_target = longitude(x,y,1);
            az = azimuth('gc',lat_radar,lon_radar,lat_target,lon_target); % Great circle azimuth
            az_math = -1*(az - 90)*(pi/180); % convert azimuth to counterclockwise from East, in radians [math form]

            for angle_idx = 1:length(elevation_angles)
                angle = elevation_angles(angle_idx);

                % Compute height of beam center
                beam_height = radar_height + dist*tand(angle);

                % If beam height falls within two interpolate-able heights,
                % compute Vr: otherwise NaN
                Vr_angled(x,y,angle_idx) = NaN;
                for z = 1:zlen
                    h = heights_full(x,y,z);
                    
                    if(h > beam_height && z > 1) % If we passed the beam height moving upward
                        lower_bound = heights_full(x,y,z-1);
                        upper_bound = h;
                        
                        % Compute proportional distance between beam center
                        % and neighbor heights
                        proportion = (beam_height - lower_bound)/(upper_bound - lower_bound);
                        
                        % Extract wind values at the chosen point
                        u = u_full(x,y,z-1)*(1-proportion) + u_full(x,y,z)*(proportion);
                        v = v_full(x,y,z-1)*(1-proportion) + v_full(x,y,z)*(proportion);
                        w = w_full(x,y,z-1)*(1-proportion) + w_full(x,y,z)*(proportion);
                        
                        % Compute Vr
                        Vr_angled(x,y,angle_idx) = u*sin(az_math)*cosd(angle) + v*cos(az_math)*cosd(angle) + w*sind(angle);
                        break; % If found the correct heights and computed, no need to keep looping for this [y x angle] instance
                    end
                end
            end
        else
            Vr_angled(x,y,:) = NaN;
        end
    end
end

% Save for easy re-access
save(sprintf(filepath_format_vr,intermediate_path,'full',type_label,wrf_file_choice),'Vr_angled');

%% 3b. Load prior saved results [OPTIONAL]

load(sprintf(filepath_format_vr,intermediate_path,'full',type_label,wrf_file_choice),'Vr_angled');

%% 3c. Plot along elevation angles

for angle_idx = 1:angle_cutoff
    angle = elevation_angles(angle_idx);

    % Plot
    figure('Position',[0 200 1240 620]); % Set dimensions of figure
    h = pcolor(longitude(:,:,1),latitude(:,:,1),Vr_angled(:,:,angle_idx)); % Plot the data
    hold on; % Allow further plotting on same figure
    scatter(lon_radar,lat_radar,80,'p','black','MarkerFaceColor','white'); % Plot star at radar station on the map
    set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
    shading interp; % Smooth out plot from grid-boxes
    c = colorbar('FontSize',axes_font_size); % Make colorbar
    if(cmap == 1)
        cmap_name = "jet";
        colormap(jet);
    else
        cmap_name = "redblue";
        colormap(redblue)
    end
    caxis([-60 60]); % Limit colorbar values for readability/consistency

    ax = gca();
    ax.FontSize = axes_font_size; % Set axis font size

    borders('continental us','black','linewidth',1); % Plot state borders
    hold off;

    xlabel('Longitude','FontSize',label_font_size); % Set other labels
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf("Radial Velocity: WRF | %s | %g* elevation | KBUF",timestamp,angle),'FontSize',title_font_size);  
    
    saveas(h,sprintf(filepath_format_plot,output_path,'full','wrf',type_label,hour,angle_idx)); % Save as .png
    close all; % Close all figure popups

end

%% 4. Regrid radial velocity onto EXRAD grid

% Read in EXRAD data
filename = ex_files(exrad_file_choice);
exrad_timestamp = extractBetween(filename,15,36);
filepath = sprintf("%s/%s","nasa_data",filename);

% Dims: lat x lon x 57 (across x along x height)
lat_ex = ncread(filepath, "latitude"); % 3D [degrees]
lon_ex = ncread(filepath, "longitude"); % 3D [degrees]
refl_ex = ncread(filepath, "Ku_band_reflectivity"); % 3D [dBZ]

% Define dimensions
lon_ex = lon_ex - 360; % Convert to negative values
dims_ex = size(lon_ex);

% Reference slices
lon_slice = lon_ex(:,:,1);
lat_slice = lat_ex(:,:,1);
refl_slice = refl_ex(:,:,1);
    
% Calculate the bounds of the real data for proper comparison with
% wrf / count how many NaN rows and columns there are
[first_data_lat,last_data_lat,first_data_lon,last_data_lon] = find_data_edges_nas(dims_ex(1),dims_ex(2),refl_slice);

% 2. Brute force bilinear interp to make the WRF data match the shape of the EXRAD data
% NOTE: Assumes that the WRF data fully encompasses the EXRAD data, so that
% edge cases can work as normal!

% Define empty grid to fit wrf data into: using radar angle instead of a
% true vertical coordinate
lon_wrf = longitude;
lat_wrf = latitude;
vr_wrf = Vr_angled;
vr_wrf_fitted = zeros(dims_ex(1),dims_ex(2),num_angles); % Radial velocity needs vertical grid too

% Loop through each EXRAD point and find corresponding WRF value at each
% radar angle therein
for x = 1:dims_ex(1)
    for y = 1:dims_ex(2)
        [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE] = find_nn_idx_irregular(lon_slice(x,y),lat_slice(x,y),lon_wrf,lat_wrf); % Get nn index values
        nn_points = [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE];
        if(any(nn_points == 0)) % If indices could not be found / reference points are NaN
            vr_wrf_fitted(x,y,:) = -30; % Lowest dBZ value used on color scale
        else
            for angle_idx = 1:num_angles
                vr_wrf_fitted(x,y,angle_idx) = bilinear_interpolate_irregular_to_grid(x,y,lon_slice,lat_slice,vr_wrf(:,:,angle_idx),lon_wrf,lat_wrf,nn_points);
            end
        end
    end
end

% Fill in NaNs on the border where they are in EXRAD
if(first_data_lat > 1)
    vr_wrf_fitted(1:first_data_lat,:,:) = NaN;
end
if(last_data_lat < dims_ex(1))
    vr_wrf_fitted(last_data_lat:dims_ex(1),:,:) = NaN;
end
if(first_data_lon > 1)
    vr_wrf_fitted(:,1:first_data_lon,:) = NaN;
end
if(last_data_lon < dims_ex(2))
    vr_wrf_fitted(:,last_data_lon:dims_ex(2),:) = NaN;
end

% Save for exterior use
save(sprintf(filepath_format_vr,intermediate_path,'fit',type_label,wrf_file_choice),'vr_wrf_fitted');

%% 4b. Load prior saved results [OPTIONAL]

load(sprintf(filepath_format_vr,intermediate_path,'fit',type_label,wrf_file_choice),'vr_wrf_fitted');

%% 4c. Plot regridded WRF Vr

% Plot
for angle_idx = 1:num_angles

    figure('Position',[0 200 2100 620]);
    h = pcolor(lon_slice,lat_slice,vr_wrf_fitted(:,:,angle_idx));
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    if(cmap == 1)
        cmap_name = "jet";
        colormap(jet);
    else
        cmap_name = "redblue";
        colormap(redblue)
    end
    caxis([-60 60]);

    ax = gca();
    ax.FontSize = axes_font_size;

    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf("WRF Vr (%s) on EXRAD grid | %g* | %s",type_label,elevation_angles(angle_idx),exrad_timestamp),'FontSize',title_font_size);

    saveas(h,sprintf(filepath_format_plot,output_path,'fit','ex',type_label,hour,angle_idx)); % Save as .png
end

close all;

%% 5. Calculate difs and RMSE between EXRAD and regridded others, plot difs

% WRF-EXRAD

load(sprintf('%s/%s',intermediate_path,exrad_comparison_filename),'Vr_angled');

vr_ex = Vr_angled;

for angle_idx = 1:num_angles
    
    sprintf('%d',angle_idx)

    rmse_wrf_ex = RMSE(vr_wrf_fitted(:,:,angle_idx),vr_ex(:,:,angle_idx));
    dif_wrf_ex = vr_wrf_fitted(:,:,angle_idx) - vr_ex(:,:,angle_idx);

    sprintf('RMSE for WRF[%s](c)-EXRAD(c) Vr dif case %d = %d',type_label,exrad_file_choice,rmse_wrf_ex)

    % Plot dif
    figure('Position',[0 200 2100 620]);
    h = pcolor(lon_slice,lat_slice,dif_wrf_ex);
    hold on;
    set(h, 'EdgeColor', 'none');
    shading interp;
    c = colorbar('FontSize',axes_font_size);
    colormap(redblue);
    caxis([-30 30]);
    ax = gca();
    ax.FontSize = axes_font_size;
    hold off;

    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    title(sprintf('WRF(%s)-EXRAD Vr Dif | %g* | %s ',type_label,elevation_angles(angle_idx),exrad_timestamp),'FontSize',title_font_size);

    saveas(h,sprintf(filepath_format_plot,output_path,'dif','ex',type_label,hour,angle_idx)); % Save as .png
    
    close all;

end

% WRF-NEXRAD: RMSE

% WRF-NEXRAD: Dif

% EXRAD-NEXRAD: RMSE

% EXRAD-NEXRAD: Dif



% END