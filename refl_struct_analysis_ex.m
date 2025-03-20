% refl_struct_analysis_ex (v2.1.2)
%
% Purpose: Generate 2D map of the vertical heights / pressure levels of the
% maximum radar reflectivity per column & vertical slices, for EXRAD
% aircraft data
% 
% Author(s): Jon Seibert
% Last updated: 25 Sept 2024
% 
% Inputs: [WRF Outputs].nc
% Outputs: [Vertical slice].png, [Height map].png
% Dependencies: reflmap.m, borders.m, linear_interpolate.m, 
%       latlon_[]_trim.mat (predefined domain boundaries), [WRF/GIS mask].mat
%
% NOTES:
%  - TIME_COUNT is not robust when crossing month boundaries.
%  - Not robust for use in southern hemisphere (lat < 0)
%  - idx_nn removed to compute more accurately
%
% TODO:
%  - 
script_name = 'refl_struct_analysis_ex.m';

local_wd = "C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code";
server_wd = "/storage/home/jjs5895/projects/IMPACTS/code";

% Check for correct execution mode
working_dir = strrep(pwd,"\","/");
if(working_dir == local_wd)
    local = 1;
elseif(working_dir == server_wd)
    local = 0;
else
    error('Unknown environment detected.')
end
   
%% 0A. Script Controls

run_name = 'vslice_9_ex';

% Local or server execution mode
%local = 1; % true = running on local device, false = running on server

data_type = "REFL"; % Working name of variable being analyzed (REFL, T, MSLP, GPH, VORT, WIND, Q, OMEGA)

% Date & time
% Range: 2020-02-07-1400 : 2020-02-07-1800; 14-00, 18-00
% 1400 is same across all experiments
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 7;
start_hour = 15;
end_hour = 15;

is_3d = 1;

% Note whether transposition or compositing is necessary
transpose_data = true;
flip_data = true;
z_composite_data = false; % Unused if use_base_ref is true

limit_raw_colorbar = false;
clim_lower = -30; % Colorbar limits
clim_upper = 75;
clim_lower_dif = -10; %  Overridden by some choices of data variables
clim_upper_dif = 10;
pressure_min = 300;
pressure_max = 1000;

% Plot selections
remake_plots = 1; % If true, overrides existing outputs. If false, only plots 'new' files.
show_plots = 1;
show_progress = 1; % If true, prints out progress bar/debug messages
announce_plots = 0;
use_auto_run_name = 0;

% Vertical slice boundaries
SW_lon = -75.7; 
SW_lat = 42; 
NE_lon = -73.4;
NW_lat = 43.8;
min_z = 1;
max_z = 50;

plot_s_slice_edge = 39;
plot_n_slice_edge = 46;
plot_w_slice_edge = -77;
plot_e_slice_edge = -71;

plot_height_maps = 0;
plot_beam_line = 0;
plot_all_beams = 1;
%use_advanced_br = 1; % This is now always enabled

radar_heights_filename = 'NEXRAD_radar_heights_wrf_grid_trim.mat';
wrf_heights_filename = 'wrf_base_refl_heights.mat';

max_dist = 230000; %(m): radar radius

z_limit = 7800; % m in altitude
print_diagnostics = 0;

xlen_max = 512;
ylen_max = 41;
zlen_max = 57;
% 512, 41, 57

make_st_locs = 0;

%% 0B. General Settings

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020';
    input_path_gis = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    input_path_nas = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/exrad_data/';
    intermediate_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/testing';
    path_to_code = "C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code";
    path_to_extra_code = './downloaded_code';
    input_path_ecmwf = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Data/ECMWF';
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/';
    input_path_gis = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/GIS';
    intermediate_path = '/storage/home/jjs5895/projects/IMPACTS/intermediate';
    output_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/ptp';
    path_to_code = "/storage/home/jjs5895/projects/IMPACTS/code";
    path_to_extra_code = '/storage/home/jjs5895/projects/IMPACTS/code/downloaded_code'; % Specify path to borders.m
    input_path_ecmwf = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/ECMWF';
end

addpath(path_to_extra_code);

hour_step = 1; % Hour increment size

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 925;
fig_height = 900;
fig_width_small = 350;
fig_height_small = 350;

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width_slice = 1225;
fig_height_slice = 500;
fig_width_slice_small = 600;
fig_height_slice_small = 300;

fig_x_d01 = 100;
fig_y_d01 = 100;
fig_width_d01 = 925;
fig_height_d01 = 675;
fig_width_d01_small = 410;
fig_height_d01_small = 300;

% Spatial limits of analysis and plot (degrees lat/lon)
w_lim = -79;
e_lim = -69.75;
s_lim = 36;
n_lim = 46;
w_lim_d01 = -97;
e_lim_d01 = -67;
s_lim_d01 = 29.5;
n_lim_d01 = 48;
limit_borders = true; % Whether to apply spatial limits to plot
limit_slice = true;
trim_e = false; % Whether to trim the eastern edge

if(local)
    % Figure font sizes
    title_font_size = 18;
    label_font_size = 18;
    axes_font_size = 16;
    contour_font_size = 14;
else
    % Figure font sizes
    title_font_size = 24;
    label_font_size = 24;
    axes_font_size = 22;
    contour_font_size = 20;
end

%dbz_nodata = -30;
%dbz_nodata_gis = -66;
nodata_val = NaN;
nodata_val_gis = NaN;

Re = 6.3781e6; % Radius pf Earth in m
g = 9.80665; % m/s^2 standard gravity at sea level
R = 287.0600676; % J/kg/K, specific gas constant

ecmwf_filename = 'ecmwf_mslp_d01_20200207.nc';
       
beam_width = 0.9; % Degrees, DIAMETER
beam_angle = 0.5;

%% NASA IMPACTS: EXRAD comparison data

% NASA Timegroups: each pair describes the relevant timestamps contained
% within (1&2, 3&4, etc)
timegroups = [135456,142903,143307,145251,145658,150453,150806,152406,152952,155034,155544,161317,161829,170931];
dategroups = [20200207135456,20200207142903,20200207143307,20200207145251,20200207145658,20200207150453,20200207150806,20200207152406,20200207152952,20200207155034,20200207155544,20200207161317,20200207161829,20200207170931];

date_strings = string(dategroups);
dt = datetime(date_strings,'InputFormat','yyyyMMddHHmmss'); % DateTimes

midpoints = [141159, 144259, 150055, 151606, 154013, 160430, 164400]; %(All with 20200207 in front)
midpoint_strings = string(midpoints);

years = ones(7,1)*2020;
months = ones(7,1)*2;
days = ones(7,1)*7;
hours = [14 14 15 15 15 16 16];
minutes = [11 42 0 16 40 04 44];
seconds = [59 59 55 06 13 30 00];

closest_hours = [14 14 15 15 16 16 17];

nasa_files = string(ls("input/exrad_data/exrad_impacts*"));

filename_format_free = 'wrfout_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00'; % domain, year, month, day, hour

% Dims: lat x lon x 57 (across x along x height)
%across_track_grid = ncread(filepath, "across_track_grid"); % 41x1 [km]
%along_track_grid = ncread(filepath, "along_track_grid"); % 517x1 [km]
%vertical_grid = ncread(filepath, "vertical_grid"); % 57x1 [km]

data_name_u = "zonal_wind"; % Wind 
data_name_v = "meridional_wind"; % Wind
data_name_w = "vertical_wind"; % Wind
data_name_refl = "Ku_band_reflectivity";
data_name_lat = "latitude";
data_name_lon = "longitude";
data_name_vertical_grid = "vertical_grid";

data_type = 'REFL';
data_src = 'nas';
exp_name = 'nas';


%% 2. Preliminary Processing & Setup

if(~show_plots)
    set(groot,'DefaultFigureVisible','off') % Turn off figure popups for local
end

% Announce operating mode
fprintf('Starting %s in %s mode.\n',script_name,upper(mode));

tic; % Start script timer

% Check for correct execution mode
current_dir = pwd;
current_dir = strrep(current_dir,"\","/");
if(current_dir ~= path_to_code)
    error("Operating mode does not match execution directory.")
end

% Lay out title and filename formats
title_format =   '[%s-%s|d0%d|%s] %s: %s(%s)|%s'; % exp name 1, exp name 2, domain, member #/mean, background/analysis/obs (bao), bao, plot_type, data type, units, timestamp
output_format =         '%s/%s_%s_d0%d_sbd%d_%s_%s_%s_%s_%s.png'; % output path, exp name 1, exp name 2, domain, subdomain, mem/mean, bao_short, plot_type, data type, timestamp
datetime_format_file =  '%04.f%02.f%02.f%02.f%02.f'; % year, month, day, hour, min
datetime_format_title = '%04.f-%02.f-%02.f-%02.f%02.f'; % year, month, day, hour, min

fprintf('Run designation: %s.\n',run_name);
if(show_progress)
    toc
    if(remake_plots)
        fprintf('WARNING: Overriding existing designated plots.\n');
    end
    fprintf('Performing initial setup.\n');
end

year = years(1);
month = months(1);
day = days(1);
hour = hours(1);
minute = minutes(1);

% Create output folders
% Directory structure:
% output/ptp/(run_name)/(variable_name)/(plot_size)
output_path_large = sprintf('%s/%s/%s/%s/%s',output_path_base,run_name,data_type,'large');
output_path_small = sprintf('%s/%s/%s/%s/%s',output_path_base,run_name,data_type,'small');

% If output folders do not exist, create them
if(~isfolder(output_path_large) || ~isfolder(output_path_small))
    mkdir(output_path_large);
    mkdir(output_path_small);
end

% REFL_10CM/REFL (dBZ), T/T (deg C), P/MSLP (mb), GPH/GPH(m),
% VORT/VORT(10^-5 s^-1), U-V-W/WIND(kt)
% Define units
switch data_type
    case "REFL"
        units = "dBZ";
    case "WIND"
        units = "kt";
        data_name = "";
        contour_steps = 6; % P: 4, T: 2, GPH: 6 (plots with GPH)
        use_contours = 0;
end

%% MAIN

num_timestamps = 4;

lon_list_box = zeros(num_timestamps,xlen_max);
lat_list_box = zeros(num_timestamps,ylen_max);
hgt_list_box = zeros(num_timestamps,1);
dims_list_box = zeros(num_timestamps,3);
slice_line_coords = zeros(num_timestamps,4);
timestamps = ["0" "0" "0" "0"];

%% 3B. For each timestamp being analyzed:
for time_idx = [1 2 3 4]

    year = years(time_idx);
    month = months(time_idx);
    day = days(time_idx);
    hour = hours(time_idx);
    minute = minutes(time_idx);
   
    % Set timestamps and experiment names
    timestamp_title = sprintf(datetime_format_title,year,month,day,hour,minute);
    timestamp_file = sprintf(datetime_format_file,year,month,day,hour,minute);
    timestamps(time_idx) = timestamp_file;
    %timestamp_file_background = sprintf(datetime_format_file,year,month,day,hour-1);
    
        
    member_string = sprintf('%03.f',0);


    filename = nasa_files(time_idx);
    filepath = sprintf("%s/%s",input_path_nas,filename);

    u = ncread(filepath, data_name_u); % 3D [m/s]
    v = ncread(filepath, data_name_v); % 3D [m/s]
    w = ncread(filepath, data_name_w); % 3D [m/s]
    refl = ncread(filepath, data_name_refl); % 3D [dBZ]
    lat_nas_full = ncread(filepath, data_name_lat); % 3D [degrees]
    lon_nas_full = ncread(filepath, data_name_lon); % 3D [degrees]
    model_height_list = ncread(filepath, data_name_vertical_grid); %km
    model_height_list = model_height_list*1000; % Convert from km to m

    lon_nas_full = lon_nas_full - 360;

    lon_nas = lon_nas_full(:,:,28);
    lat_nas = lat_nas_full(:,:,28);

    ylen_nas = size(refl,1);
    xlen_nas = size(refl,2);
    zlen_nas = size(refl,3);


    %% Calculate vertical slices: x and y-aligned only
    
    % Existing model height field name is "model_height" (above ground)
    % Load NEXRAD beam heights and WRF Base Refl height map
    %load(sprintf('%s/%s',intermediate_path,radar_heights_filename),'radar_height_grid');
    %radar_height_grid_full = radar_height_grid;
    load(sprintf('%s/%s',intermediate_path,'nexrad_stations.mat')); % wban, station_ids,_station_names,lon_deg,lat_deg,elevations,tower_heights

    %for slice_idx = 1:1
    slice_idx = time_idx;

        %clean_data = refl(~isnan(refl));

        % Slice 3: X-Z, along longitude:
        y_axis = 0;
        slice_lat = mean(lat_nas(~isnan(refl(:,:,28))),"omitnan");
        
        

        % Calculate idx_a and idx_b (indices of true grid slices
        % bracketing the desired slice lon/lat
        if(y_axis) % Slice is along y-axis
            slice_lon_list = ones(xlen_nas,1).*slice_lon;
            slice_lat_list = zeros(ylen_nas,1);
            idx_a = 0; idx_b = 0;
            y_mid = floor(ylen_nas/2);
            x_dir = 1;
            if(lon_nas(y_mid,2) < lon_nas(y_mid,1))
                x_dir = -1;
            end
            for x = 1:xlen_nas
                if(lon_nas(y_mid,x) == slice_lon)
                    idx_a = x; idx_b = x;
                    break;
                elseif(xdir == 1 && lon_nas(y_mid,x) > slice_lon)
                    idx_a = max(1,x-1); idx_b = x;
                    break;
                elseif(xdir == -1 && lon_nas(y_mid,x) < slice_lon)
                    idx_a = max(1,x-1); idx_b = x;
                    break;
                end
            end

            ab_dif = abs(lon_nas(y_mid,idx_a) - lon_nas(y_mid,idx_b));
            if(ab_dif == 0)
                wa = 0.5; wb = 0.5;
            else
                wa = abs(lon_nas(y_mid,idx_a) - slice_lon)/ab_dif;
                wb = abs(lon_nas(y_mid,idx_b) - slice_lon)/ab_dif;
            end
            for y = 1:ylen_nas
                slice_lat_list(y) = lat_nas(y,idx_a).*wb +lat_nas(y,idx_b).*wa;
            end

        else % Slice is along x-axis
            slice_lon_list = zeros(xlen_nas,1);
            slice_lat_list = ones(ylen_nas,1).*slice_lat;
            idx_a = 0; idx_b = 0;
            x_mid = floor(xlen_nas/2);
            y_dir = 1;
            if(lat_nas(2,x_mid) < lat_nas(1,x_mid))
                y_dir = -1;
            end
            for y = 1:ylen_nas
                if(lat_nas(y,x_mid) == slice_lat)
                    idx_a = y; idx_b = y; 
                    break;
                elseif(y_dir == 1 && lat_nas(y,x_mid) > slice_lat)
                    idx_a = max(1,y-1); idx_b = y;
                    break;
                elseif(y_dir == -1 && lat_nas(y,x_mid) < slice_lat)
                    idx_a = max(1,y-1); idx_b = y;
                    break;
                end
            end

            ab_dif = abs(lat_nas(idx_a,x_mid) - lat_nas(idx_b,x_mid));
            if(ab_dif == 0)
                wa = 0.5; wb = 0.5;
            else
                wa = abs(lat_nas(idx_a,x_mid) - slice_lat)/ab_dif;
                wb = abs(lat_nas(idx_b,x_mid) - slice_lat)/ab_dif;
            end
            
            for x = 1:xlen_nas
                slice_lon_list(x) = lon_nas(idx_a,x)*wb + lon_nas(idx_b,x)*wa;
            end
            
        end

        % Remove NaNs from latlon lists
        clean_lon_idx = find(~isnan(slice_lon_list));
        slice_lon_list_full = slice_lon_list;
        slice_lon_list = slice_lon_list(clean_lon_idx);

        clean_lat_idx = find(~isnan(slice_lat_list));
        slice_lat_list_full = slice_lat_list;
        slice_lat_list = slice_lat_list(clean_lat_idx);

        xlen_nas_clean = length(slice_lon_list);
        ylen_nas_clean = length(slice_lat_list);
        zlen_nas_clean = zlen_nas;

        slice_lat_grid = ones(ylen_nas_clean,zlen_nas_clean);
        slice_lon_grid = ones(xlen_nas_clean,zlen_nas_clean);

        for z = 1:zlen_nas_clean
            slice_lat_grid(:,z) = slice_lat_list;
            slice_lon_grid(:,z) = slice_lon_list;
        end

        data = refl(clean_lat_idx,clean_lon_idx,:);
        lat_nas_clean = lat_nas(clean_lat_idx,clean_lon_idx,:);
        lon_nas_clean = lon_nas(clean_lat_idx,clean_lon_idx,:);

        if(y_axis)
            slice_values = zeros(ylen_nas_clean,zlen_nas_clean);
            slice_line_coords(time_idx,1) = slice_lon;
            slice_line_coords(time_idx,2) = slice_lat_list(1);
            slice_line_coords(time_idx,3) = slice_lon;
            slice_line_coords(time_idx,4) = slice_lat_list(ylen_nas_clean);
        else
            slice_values = zeros(xlen_nas_clean,zlen_nas_clean);
            slice_line_coords(time_idx,1) = slice_lon_list(1);
            slice_line_coords(time_idx,2) = slice_lat;
            slice_line_coords(time_idx,3) = slice_lon_list(xlen_nas_clean);
            slice_line_coords(time_idx,4) = slice_lat;
        end

        % Set up proper grid for "model" height (altitude)
        model_height = ones(ylen_nas_clean,xlen_nas_clean,zlen_nas_clean);
        for z = 1:zlen_nas_clean
            model_height(:,:,z) = model_height(:,:,z).*model_height_list(z);
        end

        % Recalculate idx_a and idx_b (indices of true grid slices
        % bracketing the desired slice lon/lat
        if(y_axis) % Slice is along y-axis
            idx_a = 0; idx_b = 0;
            y_mid = floor(ylen_nas_clean/2);
            y_dir = 1;
            if(lon_nas_clean(y_mid,2) < lon_nas_clean(y_mid,1))
                y_dir = -1;
            end
            for x = 1:xlen_nas_clean
                if(lon_nas_clean(y_mid,x) == slice_lon)
                    idx_a = x; idx_b = x;
                    break;
                elseif(ydir == 1 && lon_nas_clean(y_mid,x) > slice_lon)
                    idx_a = max(1,x-1); idx_b = x;
                    break;
                elseif(ydir == -1 && lon_nas_clean(y_mid,x) < slice_lon)
                    idx_a = max(1,x-1); idx_b = x;
                    break;
                end
            end
        else % Slice is along x-axis
            idx_a = 0; idx_b = 0;
            x_mid = floor(xlen_nas_clean/2);
            x_dir = 1;
            if(lat_nas(2,x_mid) < lat_nas(1,x_mid))
                x_dir = -1;
            end
            for y = 1:ylen_nas
                if(lat_nas(y,x_mid) == slice_lat)
                    idx_a = y; idx_b = y; 
                    break;
                elseif(x_dir == 1 && lat_nas(y,x_mid) > slice_lat)
                    idx_a = max(1,y-1); idx_b = y;
                    break;
                elseif(x_dir == -1 && lat_nas(y,x_mid) < slice_lat)
                    idx_a = max(1,y-1); idx_b = y;
                    break;
                end
            end
        end

        %% Compute radar beam intersections

        % Existing model height field name is "model_height" (above ground)
        load(sprintf('%s/%s',intermediate_path,radar_heights_filename),'radar_height_grid');
        %load(sprintf('%s/%s',intermediate_path,wrf_heights_filename),'wrf_height_map');
        radar_height_grid_full = radar_height_grid;

        %D = haversine_distance(rad_lon,rad_lat,wrf_lon_val,wrf_lat_val);
        %R = D/(cosd(beam_angle));
        %H = sqrt(R^2 + ((4/3)*Re)^2 + 2*R*(4/3)*Re*sind(beam_angle)) - (4/3)*Re + rad_elev + rad_th; % Keenan version adjusted for m

        % Find the radars that intersect with the slice plane
        num_radars_all = size(radar_height_grid_full,1);
        num_radars = 0;
        radar_checklist = zeros(num_radars_all,1);
        for rad_idx = 1:num_radars_all
            if(y_axis)
                for lat_idx = 1:length(slice_lat_list)
                    slice_lat = slice_lat_list(lat_idx);
                    D = haversine_distance(slice_lon,slice_lat,lon_deg(rad_idx),lat_deg(rad_idx));
                    if(D < max_dist)
                        % Mark as intersecting somewhere
                        num_radars = num_radars+1;
                        radar_checklist(rad_idx) = 1;
                        break;
                    end
                end
            else
                for lon_idx = 1:length(slice_lon_list)
                    slice_lon = slice_lon_list(lon_idx);
                    D = haversine_distance(slice_lon,slice_lat,lon_deg(rad_idx),lat_deg(rad_idx));
                    if(D < max_dist)
                        % Mark as intersecting somewhere
                        num_radars = num_radars+1;
                        radar_checklist(rad_idx) = 1;
                        break;
                    end
                end
            end
        end

        radars_used = find(radar_checklist == 1);

        %% Compute aircraft-localized radar height grid
    
        radar_height_grid = zeros(num_radars,ylen_nas_clean,xlen_nas_clean);
        radar_height_grid(:,:,:) = NaN;
        for station_idx = 1:num_radars
            % Retrieve data
            rad_lon = lon_deg(radars_used(station_idx));
            rad_lat = lat_deg(radars_used(station_idx));
            rad_id = station_ids(radars_used(station_idx));
            rad_name = station_names(radars_used(station_idx));
            rad_elev = elevations(radars_used(station_idx));
            rad_th = tower_heights(radars_used(station_idx));
    
            % Compute wrf grid indices and beam height
            for y = 1:ylen_nas_clean
                for x = 1:xlen_nas_clean
                    D = haversine_distance(rad_lon,rad_lat,lon_nas_clean(y,x),lat_nas_clean(y,x));
                    if(D <= max_dist)
                        % Compute beam height H from slant distance R and ground distance D
                        R = D/(cosd(beam_angle));
                        %H = R*sind(beam_angle) + (R^2)/(2*(4/3)*Re) + rad_elev + rad_th; % Jon version
                        H = sqrt(R^2 + ((4/3)*Re)^2 + 2*R*(4/3)*Re*sind(beam_angle)) - (4/3)*Re + rad_elev + rad_th; % Keenan version adjusted for m
                        radar_height_grid(station_idx,y,x) = H;
                    end
                end
            end
        end

        %% 
        % Find each intersection
        if(y_axis)
            intersection_height_mult = zeros(num_radars,ylen_nas_clean);
            refl_value_seen = zeros(size(intersection_height_mult));
            model_height_slice = zeros(size(slice_values));
            for rad_idx = 1:num_radars
                for y = 1:ylen_nas_clean
                    rad_height = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,radar_height_grid(rad_idx,y,idx_a),radar_height_grid(rad_idx,y,idx_b));
                    intersection_height_mult(rad_idx,y) = rad_height;    
                    for z = 1:zlen_nas_clean
                        model_height_slice(y,z) = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,model_height(y,idx_a,z),model_height(y,idx_b,z));
                        slice_values(y,z) = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,data(y,idx_a,z),data(y,idx_b,z));
                    end
                end  
            end              
        else
            intersection_height_mult = zeros(num_radars,xlen_nas_clean);
            refl_value_seen = zeros(size(intersection_height_mult));
            model_height_slice = zeros(size(slice_values));
            for rad_idx = 1:num_radars
                for x = 1:xlen_nas_clean
                    rad_height = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,radar_height_grid(rad_idx,idx_a,x),radar_height_grid(rad_idx,idx_b,x));
                    intersection_height_mult(rad_idx,x) = rad_height;
                    for z = 1:zlen_nas_clean
                        model_height_slice(x,z) = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,model_height(idx_a,x,z),model_height(idx_b,x,z));
                        slice_values(x,z) = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,data(idx_a,x,z),data(idx_b,x,z));
                    end
                end  
            end    
        end

        %% Compute real base reflectivity heights
        
        %if(use_advanced_br)
        
            % Determine which radar refl value is greatest, take that
            % beam's intersection height and compute pressure
            if(y_axis) % Slice is along y-axis
                intersection_height = zeros(ylen_nas_clean,1);
                intersection_height(:) = NaN;
                for y = 1:ylen_nas_clean
                    max_dbz_idx = 0;
                    max_dbz = -999;
                    for rad_idx = 1:num_radars
                        rad_height = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,radar_height_grid(rad_idx,y,idx_a),radar_height_grid(rad_idx,y,idx_b));
                        for z = 1:zlen_nas_clean
                            if(model_height_slice(y,z) > rad_height)
                                rad_dbz = linear_interpolate(model_height_slice(y,max(1,z-1)),model_height_slice(y,z),rad_height,slice_values(y,max(1,z-1)),slice_values(y,z));
                                refl_value_seen(rad_idx,y) = rad_dbz;
                                break;
                            end
                        end
                        if(rad_dbz > max_dbz)
                            max_dbz = rad_dbz;
                            max_dbz_idx = rad_idx;
                        end
                    end
                    if(max_dbz_idx > 0)
                        intersection_height(y) = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,radar_height_grid(max_dbz_idx,y,idx_a),radar_height_grid(max_dbz_idx,y,idx_b));
                    end
                end
            else % Slice is along x-axis
                intersection_height = zeros(xlen_nas_clean,1);
                intersection_height(:) = NaN;
                for x = 1:xlen_nas_clean
                    max_dbz_idx = 0;
                    max_dbz = -999;
                    for rad_idx = 1:num_radars
                        rad_height = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,radar_height_grid(rad_idx,idx_a,x),radar_height_grid(rad_idx,idx_b,x));
                        for z = 1:zlen_nas_clean
                            if(model_height_slice(x,z) > rad_height)
                                rad_dbz = linear_interpolate(model_height_slice(x,max(1,z-1)),model_height_slice(x,z),rad_height,slice_values(x,max(1,z-1)),slice_values(x,z));
                                refl_value_seen(rad_idx,x) = rad_dbz;
                                break;
                            end
                        end
                        if(rad_dbz > max_dbz)
                            max_dbz = rad_dbz;
                            max_dbz_idx = rad_idx;
                        end
                    end
                    if(max_dbz_idx > 0)
                        intersection_height(x) = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,radar_height_grid(max_dbz_idx,idx_a,x),radar_height_grid(max_dbz_idx,idx_b,x));
                    end
                end
            end
        %end

        % Test: determine max height needed for plots
        if(print_diagnostics)
            highest_data = 0;
            for z = 1:zlen_nas_clean
                if(any(slice_values(:,z)))
                    highest_data = z;
                else
                    break;
                end
            end
            fprintf('Timestamp: %s\n',timestamp_title)
            fprintf('Maximum data altitude: z=%d, %f m\n',z,model_height_slice(1,z));

            % ---------------

            fprintf('Dimensions: X %d, Y %d, Z %d\n',xlen_nas_clean,ylen_nas_clean,zlen_nas_clean);
        end

        lon_list_box(time_idx,(1:xlen_nas_clean)) = slice_lon_list;
        lat_list_box(time_idx,(1:ylen_nas_clean)) = slice_lat_list;
        hgt_list_box(time_idx,1) = z_limit;
        dims_list_box(time_idx,:) = [ylen_nas_clean,xlen_nas_clean,zlen_nas_clean];

        %% Plot vertical slices

        %beam_color = "#757575";
        %beam_color = "#0076ff";
        %beam_color = "#ff00f3";
        beam_color = "#000000";
        multibeam_color = "#ff00f3";

        exp_name_plot = 'EXRAD';
        bao_short = 'oo';

        plot_type = 'vsl';

        if(plot_all_beams)
            plot_type = strcat(plot_type,'-all');
        end
        if(plot_beam_line)
            plot_type = strcat(plot_type,'-base');
        end

        data_to_plot = slice_values;
        cmap = 'reflmap';
        output_format_slice = '%s/%s_%s_%s_%s_%s_s%d_%s.png'; % path, exp name a, exp name b, domain, subdomain, member string, bao, plot type, data type, slice idx, timestamp
        plot_filename = sprintf(output_format_slice,output_path_large,lower(exp_name_plot),lower(exp_name_plot),bao_short,plot_type,data_type,slice_idx,timestamp_file);
        plot_filename_small = sprintf(output_format_slice,output_path_small,lower(exp_name_plot),lower(exp_name_plot),bao_short,plot_type,data_type,slice_idx,timestamp_file);

        if(isfile(plot_filename) && ~remake_plots) % If shouldn't override existing plots, don't
            continue;
        end

        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width_slice fig_height_slice]); % Create initial blank figure

        if(show_progress && announce_plots)
            toc
            fprintf('Plotting %s...\n',plot_type); % debug statement
        end


        if(y_axis)
            h = pcolor(slice_lat_grid,model_height_slice,data_to_plot); % Plot the data
        else
            h = pcolor(slice_lon_grid,model_height_slice,data_to_plot); % Plot the data
        end

        set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
        shading interp; % Smooth out plot from grid-boxes
        c = colorbar('FontSize',axes_font_size); % Make colorbar
        colormap(cmap); % Set colors
        ylabel(c,'Radar Reflectivity (dBZ)');
        caxis([clim_lower clim_upper]);
        %set (gca,'YDir','reverse')

        hold on;

        if(plot_all_beams)
            if(y_axis)
                for rad_idx = 1:num_radars
                    h = plot(slice_lat_list,intersection_height_mult(rad_idx,:,:),'Color',multibeam_color,'LineWidth',2,'LineStyle','--');
                    %h = scatter(slice_lat_list,intersection_pressure_mult(rad_idx,:,:),5,'MarkerEdgeColor',multibeam_color,'MarkerFaceColor',multibeam_color);
                end
            else
                for rad_idx = 1:num_radars
                    h = plot(slice_lon_list,intersection_height_mult(rad_idx,:,:),'Color',multibeam_color,'LineWidth',2,'LineStyle','--');
                end
            end
        end
        if(plot_beam_line)
            if(y_axis)
                h = plot(slice_lat_list,intersection_height,'Color',beam_color,'LineWidth',2,'LineStyle','-');
            else
                h = plot(slice_lon_list,intersection_height,'Color',beam_color,'LineWidth',2,'LineStyle','-');
            end
        end
        if(limit_slice)
                ylim([0 z_limit]);
        end

        % Apply labels

        if(y_axis)
            title(sprintf('[%s] Refl VSlice #%d (Lon %02.02f) | %s',exp_name_plot,slice_idx,slice_lon,timestamp_title),'FontSize',title_font_size); 
            xlabel('Latitude (deg)','FontSize',label_font_size);
            ylabel('Altitude (m)','FontSize',label_font_size);
        else
            title(sprintf('[%s] Refl VSlice #%d (Lat %02.02f) | %s',exp_name_plot,slice_idx,slice_lat,timestamp_title),'FontSize',title_font_size); 
            xlabel('Longitude (deg)','FontSize',label_font_size);
            ylabel('Altitude (m)','FontSize',label_font_size);
        end
        set(gca,'Fontsize',axes_font_size);

        if(show_progress && announce_plots)
            toc
            fprintf('Saving...\n'); % debug statement
        end

        saveas(gcf,plot_filename); % Save as .png

        if(show_progress && announce_plots)
            toc
            fprintf('Making small...\n'); % debug statement
        end

        % Plot SMALL
        f.Position = [fig_x fig_y fig_width_slice_small fig_height_slice_small]; % Shrink figure
        title(sprintf('[%s] Refl VSL #%d',exp_name_plot,slice_idx),'FontSize',title_font_size); 

        saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
        close('all');

        clearvars f h data_to_plot;
end

%% Plot slice locations

plot_type = "slloc";
%output_format_slice = '%s/%s_%s_%s_%s_%s_s%d_%s.png'; % path, exp name a, exp name b, domain, subdomain, member string, bao, plot type, data type, slice idx, timestamp
plot_filename = sprintf(output_format_slice,output_path_large,lower(exp_name_plot),lower(exp_name_plot),bao_short,plot_type,data_type,slice_idx,timestamp_file);
plot_filename_small = sprintf(output_format_slice,output_path_small,lower(exp_name_plot),lower(exp_name_plot),bao_short,plot_type,data_type,slice_idx,timestamp_file);

colors = ["red" "#438d02" "blue" "#8c01d6"];

% Plot LARGE
f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure

hold on;
borders('continental us','black','linewidth',1); 

for plot_idx = 1:num_timestamps
    plot(slice_line_coords(plot_idx,1:2:3),slice_line_coords(plot_idx,2:2:4),'Color',colors(plot_idx),'LineWidth',3);
end

save("intermediate/exrad_slice_line_coords.mat","slice_line_coords");

limit_borders = 1;
if(limit_borders)
    xlim([w_lim e_lim]);
    ylim([s_lim n_lim]);
end


% TEMP PLOT
if(make_st_locs)
    load('st_latlon.mat','st_lon_list','st_lat_list','st_name_list');
    h = scatter(st_lon_list,st_lat_list,200,"yellow","filled","pentagram");
    h.MarkerEdgeColor = "#000000";
    h.MarkerFaceColor = "#f4dc00";
    hold off;
    
    st_lon_list_labels = st_lon_list;
    st_lat_list_labels = st_lat_list;
    st_lon_list_labels(1) = st_lon_list(1) - 0.1;
    st_lat_list_labels(1) = st_lat_list(1) - 0.3;
    st_lon_list_labels(2) = st_lon_list(2) + 0.6;
    st_lat_list_labels(2) = st_lat_list(2) - 0.35;
    st_lon_list_labels(3) = st_lon_list(3) - 0.1;
    st_lat_list_labels(3) = st_lat_list(3) - 0.3;
    st_lon_list_labels(4) = st_lon_list(4) + 0.9;
    st_lat_list_labels(4) = st_lat_list(4) + 0.12;
    
    title(sprintf('EXRAD Vertical Slice Locations'),'FontSize',title_font_size); 
    xlabel('Longitude (deg)','FontSize',label_font_size);
    ylabel('Latitude (deg)','FontSize',label_font_size);
    legend(["" "1411" "1442" "1500" "1516"]);
    set(gca,'Fontsize',axes_font_size);
    labelpoints(st_lon_list_labels,st_lat_list_labels,st_name_list,'FontSize',axes_font_size-5,'FontWeight','bold');
    
    saveas(gcf,plot_filename); % Save as .png
end

%%

save('intermediate/nas_slice_dims.mat','lon_list_box','lat_list_box','hgt_list_box','dims_list_box');

% Get total runtime and print to stdout
runtime = toc;
hours = floor(runtime/3600);
mins = floor((runtime/60) - (hours*60));
secs = toc - (hours*3600) - (mins*60);
fprintf('Done. Total script runtime = %02.f:%02.f:%02.f\n',hours,mins,secs)

% END      