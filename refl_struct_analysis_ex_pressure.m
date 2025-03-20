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

run_name = 'vslice_2';

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
start_hour = 14;
end_hour = 18;

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
show_plots = 0;
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
plot_beam_line = 1;
plot_all_beams = 1;
%use_advanced_br = 1; % This is now always enabled

radar_heights_filename = 'NEXRAD_radar_heights_wrf_grid_trim.mat';
wrf_heights_filename = 'wrf_base_refl_heights.mat';

max_dist = 230000; %(m): radar radius

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
datetime_format_file =  '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

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

%% 3B. For each timestamp being analyzed:
for time_idx = [1 2 3 4]

    year = years(time_idx);
    month = months(time_idx);
    day = days(time_idx);
    hour = hours(time_idx);
    minute = minutes(time_idx);
   
    % Set timestamps and experiment names
    timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
    timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
    timestamp_file_background = sprintf(datetime_format_file,year,month,day,hour-1);
    
        
    member_string = sprintf('%03.f',0);


    filename = nasa_files(time_idx);
    filepath = sprintf("%s/%s",input_path_nas,filename);

    u = ncread(filepath, data_name_u); % 3D [m/s]
    v = ncread(filepath, data_name_v); % 3D [m/s]
    w = ncread(filepath, data_name_w); % 3D [m/s]
    refl = ncread(filepath, data_name_refl); % 3D [dBZ]
    lat_nas_full = ncread(filepath, data_name_lat); % 3D [degrees]
    lon_nas_full = ncread(filepath, data_name_lon); % 3D [degrees]
    model_height = ncread(filepath, data_name_vertical_grid); %km
    model_height = model_height*1000; % Convert from km to m

    lon_nas_full = lon_nas_full - 360;

    lon_nas = lon_nas_full(:,:,28);
    lat_nas = lat_nas_full(:,:,28);

    load(sprintf('%s/%s',intermediate_path,sprintf('mean_storage_CONV_20200207%d00.mat',closest_hours(time_idx))),'data_p','model_height');
    load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'lon_wrf_trim','lat_wrf_trim');
    height_grid_wrf = model_height;
    pressure_wrf = data_p;
    lon_wrf = lon_wrf_trim;
    lat_wrf = lat_wrf_trim;

    ylen_nas = size(refl,1);
    xlen_nas = size(refl,2);
    zlen_nas = size(refl,3);
    ylen_wrf = size(pressure_wrf,1);
    xlen_wrf = size(pressure_wrf,2);
    zlen_wrf = size(pressure_wrf,3);


    %% Plot vertical slices: x and y-aligned only
    
    % Existing model height field name is "model_height" (above ground)
    % Load NEXRAD beam heights and WRF Base Refl height map
    %load(sprintf('%s/%s',intermediate_path,radar_heights_filename),'radar_height_grid');
    %radar_height_grid_full = radar_height_grid;
    load(sprintf('%s/%s',intermediate_path,'nexrad_stations.mat')); % wban, station_ids,_station_names,lon_deg,lat_deg,elevations,tower_heights



    %%

    for slice_idx = 1:1

        clean_data = refl(~isnan(refl));

        % Slice 3: X-Z, along longitude:
        y_axis = 0;
        slice_lat = mean(lat_nas(~isnan(refl(:,:,28))),"omitnan");
        

        if(y_axis) % Slice is along y-axis
            slice_lon_list = ones(xlen_nas,1).*slice_lon;
            slice_lat_list = zeros(ylen_nas,1);
            idx_a = 0; idx_b = 0;
            y_mid = floor(ylen_nas/2);
            for x = 1:xlen_nas
                if(lon_nas(y_mid,x) == slice_lon)
                    idx_a = x; idx_b = x;
                    break;
                elseif(lon_nas(y_mid,x) > slice_lon)
                    idx_a = max(1,x-1); idx_b = x;
                    break;
                end
            end

            idx_a_wrf = 0; idx_b_wrf = 0;
            y_mid_wrf = floor(ylen_wrf/2);
            for x = 1:xlen_wrf
                if(lon_wrf(y_mid_wrf,x) == slice_lon)
                    idx_a_wrf = x; idx_b_wrf = x;
                    break;
                elseif(lon_wrf(y_mid_wrf,x) > slice_lon)
                    idx_a_wrf = max(1,x-1); idx_b_wrf = x;
                    break;
                end
            end

            pressure_slice_a = squeeze(pressure_wrf(:,idx_a_wrf,:));
            pressure_slice_b = squeeze(pressure_wrf(:,idx_b_wrf,:));
            height_slice_a = squeeze(height_wrf(:,idx_a_wrf,:));
            height_slice_b = squeeze(height_wrf(:,idx_b_wrf,:));
            wa = abs(lon_wrf(y_mid_wrf,idx_a_wrf) - slice_lon)/abs(lon_wrf(y_mid_wrf,idx_a_wrf) - lon_wrf(y_mid_wrf,idx_b_wrf));
            wb = abs(lon_wrf(y_mid_wrf,idx_b_wrf) - slice_lon)/abs(lon_wrf(y_mid_wrf,idx_a_wrf) - lon_wrf(y_mid_wrf,idx_b_wrf));
            pressure_slice_wrf = pressure_slice_a.*wb + pressure_slice_b.*wa;
            height_slice_wrf = height_slice_a.*wb + height_slice_b.*wa;


            pressure_max = max(pressure_slice_wrf,[],"all");

            y_mid = floor(ylen_nas/2);
            for y = 1:ylen_nas
                wa = abs(lon_nas(y_mid,idx_a) - slice_lon)/abs(lon_nas(y_mid,idx_a) - lon_nas(y_mid,idx_b));
                wb = abs(lon_nas(y_mid,idx_b) - slice_lon)/abs(lon_nas(y_mid,idx_a) - lon_nas(y_mid,idx_b));
                slice_lat_list(y) = lat_nas(y,idx_a).*wb +lat_nas(y,idx_b).*wa;
            end

        else % Slice is along x-axis
            slice_lon_list = zeros(xlen_nas,1);
            slice_lat_list = ones(ylen_nas,1).*slice_lat;
            idx_a = 0; idx_b = 0;
            x_mid = floor(xlen_nas/2);
            for y = 1:ylen_nas
                if(lat_nas(y,x_mid) == slice_lat)
                    idx_a = y; idx_b = y; 
                    break;
                elseif(lat_nas(y,x_mid) > slice_lat)
                    idx_a = max(1,y-1); idx_b = y;
                    break;
                end
            end

            idx_a_wrf = 0; idx_b_wrf = 0;
            x_mid_wrf = floor(xlen_wrf/2);
            y_mid_wrf = floor(ylen_wrf/2);
            for y = 1:ylen_wrf
                if(lat_wrf(y,x_mid) <= slice_lat)
                    y_mid_wrf = y;
                    break;
                end
            end
            for x = 1:xlen_wrf
                if(lon_wrf(y_mid_wrf,x) >= mean(lon_nas,"all","omitnan"))
                    x_mid_wrf = x;
                    break;
                end
            end
            for y = ylen_wrf:-1:1
                if(lat_wrf(y,x_mid_wrf) == slice_lat)
                    idx_a_wrf = y; idx_b_wrf = y; 
                    break;
                elseif(lat_wrf(y,x_mid_wrf) > slice_lat)
                    idx_a_wrf = max(1,y-1); idx_b_wrf = y;
                    break;
                end
            end

            pressure_slice_a = squeeze(pressure_wrf(idx_a_wrf,:,:));
            pressure_slice_b = squeeze(pressure_wrf(idx_b_wrf,:,:));
            height_slice_a = squeeze(height_grid_wrf(idx_a_wrf,:,:));
            height_slice_b = squeeze(height_grid_wrf(idx_b_wrf,:,:));
            wa_wrf = abs(lat_wrf(idx_a_wrf,x_mid_wrf) - slice_lat)/abs(lat_wrf(idx_a_wrf,x_mid_wrf) - lat_wrf(idx_b_wrf,x_mid_wrf));
            wb_wrf = abs(lat_wrf(idx_b_wrf,x_mid_wrf) - slice_lat)/abs(lat_wrf(idx_a_wrf,x_mid_wrf) - lat_wrf(idx_b_wrf,x_mid_wrf));
            pressure_slice_wrf = pressure_slice_a.*wb_wrf + pressure_slice_b.*wa_wrf;
            height_slice_wrf = height_slice_a.*wb_wrf + height_slice_b.*wa_wrf;

            pressure_max = max(pressure_slice_wrf,[],"all");

            x_mid = floor(xlen_nas/2);
            wa = abs(lat_nas(idx_a,x_mid) - slice_lat)/abs(lat_nas(idx_a,x_mid) - lat_nas(idx_b,x_mid));
            wb = abs(lat_nas(idx_b,x_mid) - slice_lat)/abs(lat_nas(idx_a,x_mid) - lat_nas(idx_b,x_mid));
            for x = 1:xlen_nas
                slice_lon_list(x) = lon_nas(idx_a,x)*wb + lon_nas(idx_b,x)*wa;
            end
        end

        clean_lon_idx = find(~isnan(slice_lon_list));
        slice_lon_list_full = slice_lon_list;
        slice_lon_list = slice_lon_list(clean_lon_idx);

        clean_lat_idx = find(~isnan(slice_lat_list));
        slice_lat_list_full = slice_lat_list;
        slice_lat_list = slice_lat_list(clean_lat_idx);

        xlen_nas_clean = length(slice_lon_list);
        ylen_nas_clean = length(slice_lat_list);
        zlen_nas_clean = zlen_nas;

        slice_lat_grid = ones(zlen_nas_clean,ylen_nas_clean);
        slice_lon_grid = ones(zlen_nas_clean,xlen_nas_clean);

        for z = 1:zlen_nas_clean
            slice_lat_grid(z,:) = slice_lat_list;
            slice_lon_grid(z,:) = slice_lon_list;
        end

        data = refl(clean_lat_idx,clean_lon_idx,:);
        lat_nas_clean = lat_nas(clean_lat_idx,clean_lon_idx,:);
        lon_nas_clean = lon_nas(clean_lat_idx,clean_lon_idx,:);

        if(y_axis)
            slice_values = zeros(ylen_nas_clean,zlen_nas_clean);
        else
            slice_values = zeros(xlen_nas_clean,zlen_nas_clean);
        end

        height_grid_nas = ones(ylen_nas_clean,xlen_nas_clean,zlen_nas_clean);
        for z = 1:zlen_nas_clean
            height_grid_nas(:,:,z) = height_grid_nas(:,:,z).*model_height(z);
        end

        % Use bilinear interpolation to create a nas-grid pressure and
        % height field for the slice
        pressure_slice = zeros(size(slice_values));
        height_slice = zeros(size(slice_values));

        if(~y_axis)
            for x = 1:xlen_nas_clean
                cur_lat = slice_lat;
                cur_lon = slice_lon_list(x);
                for z = 1:zlen_nas_clean
                    cur_hgt = height_grid_nas(idx_a,x,z)*wb + height_grid_nas(idx_b,x,z)*wa;
                    height_slice(x,z) = cur_hgt;
                    x_idx_wrf = 0;
                    z_idx_wrf_a = 0;
                    z_idx_wrf_b = 0;
                    for x_wrf = 1:xlen_wrf
                        % Find approximate longitude of WRF at this x from
                        % the idx_a_Wrf and idx_b_wrf???
                        cur_lon_wrf = linear_interpolate();
                        if(cur_lon_wrf >= cur_lon)
                            %;
                            % Find neighbor indices for z, to interpolate
                            % below
                            for z_wrf = 1:zlen_wrf
    
                            end
                        end
                        
                    end
                    pressure_slice(x,z) = linear_interpolate(height_slice_wrf(x_idx_wrf,z_idx_wrf_a),height_slice_wrf(x_idx_wrf,z_idx_wrf_b),cur_hgt,pressure_slice_wrf(x_idx_wrf,z_idx_wrf_a),pressure_slice_wrf(x_idx_wrf,z_idx_wrf_b));
                    %height_slice_wrf(x,z) = 0;
                    %pressure_slice_wrf(x,z) = 0;
                end
            end
        
        end

        %% Compute radar beam intersections

        % Existing model height field name is "model_height" (above ground)
        %load(sprintf('%s/%s',intermediate_path,radar_heights_filename),'radar_height_grid');
        %load(sprintf('%s/%s',intermediate_path,wrf_heights_filename),'wrf_height_map');

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

        % Find each intersection
        if(y_axis)
            intersection_pressure_mult = zeros(num_radars,ylen_nas_clean);
            refl_value_seen = zeros(size(intersection_pressure_mult));
            model_height_interp = zeros(size(slice_values));
            data_interp = zeros(size(slice_values));
            for rad_idx = 1:num_radars
                for y = 1:ylen_nas_clean
                    rad_height = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,radar_height_grid(rad_idx,y,idx_a),radar_height_grid(rad_idx,y,idx_b));
                    for z = 1:zlen_nas_clean
                        model_height_interp(y,z) = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,height_grid_nas(y,idx_a,z),height_grid_nas(y,idx_b,z));
                        data_interp(y,z) = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,data(y,idx_a,z),data(y,idx_b,z));
                        if(model_height_interp(y,z) > rad_height)
                            intersection_pressure_mult(rad_idx,y) = linear_interpolate(model_height_interp(y,max(z-1,1)),model_height_interp(y,z),rad_height,pressure_slice_wrf(y,max(z-1,1)),pressure_slice_wrf(y,z));
                            break;
                        end
                    end
                end  
            end              
        else
            intersection_pressure_mult = zeros(num_radars,xlen_nas_clean);
            refl_value_seen = zeros(size(intersection_pressure_mult));
            model_height_interp = zeros(size(slice_values));
            data_interp = zeros(size(slice_values));
            for rad_idx = 1:num_radars
                for x = 1:xlen_nas_clean
                    p_idx_a = 0;
                    p_idx_b = 0;
                    for x_wrf = 1:xlen_wrf
                        p_idx_a = x_wrf;
                        p_idx_b = max(1,x_wrf-1);
                        if(lon_wrf(idx_a_wrf,x_wrf) >= lon_nas_clean(idx_a,x,1))
                            break;
                        end
                    end
                    
                    rad_height = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,radar_height_grid(rad_idx,idx_a,x),radar_height_grid(rad_idx,idx_b,x));
                    for z = 1:zlen_nas_clean
                        model_height_interp(x,z) = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,height_grid_nas(idx_a,x,z),height_grid_nas(idx_b,x,z));
                        data_interp(x,z) = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,data(idx_a,x,z),data(idx_b,x,z));
                        if(model_height_interp(x,z) > rad_height)

                            %% THIS DOESNT WORK- USING NAS Z AXIS ON WRF GRIDS

                            pressure_a = linear_interpolate(lon_wrf(idx_a_wrf,p_idx_a),lon_wrf(idx_a_wrf,p_idx_b),lon_nas_clean(idx_a,x,max(z-1,1)),pressure_slice_wrf(p_idx_a,max(z-1,1)),pressure_slice_wrf(p_idx_b,max(z-1,1)));          
                            pressure_b = linear_interpolate(lon_wrf(idx_a_wrf,p_idx_a),lon_wrf(idx_a_wrf,p_idx_b),lon_nas_clean(idx_a,x,z),pressure_slice_wrf(p_idx_a,z),pressure_slice_wrf(p_idx_b,z));          
                            intersection_pressure_mult(rad_idx,x) = linear_interpolate(model_height_interp(x,max(z-1,1)),model_height_interp(x,z),rad_height,pressure_a,pressure_b);
                            break;
                        end
                    end
                end  
            end    
        end

        %% Compute real base reflectivity heights
        
        %if(use_advanced_br)
        
            % Determine which radar refl value is greatest, take that
            % beam's intersection height and compute pressure
            if(y_axis) % Slice is along y-axis
                height_values = zeros(ylen_nas_clean,1);
                intersection_pressure = zeros(ylen_nas_clean,1);
                intersection_pressure(:) = NaN;
                for y = 1:ylen_nas_clean
                    max_dbz_idx = 0;
                    max_dbz = -999;
                    for rad_idx = 1:num_radars
                        rad_height = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,radar_height_grid(rad_idx,y,idx_a),radar_height_grid(rad_idx,y,idx_b));
                        for z = 1:zlen_nas_clean
                            if(model_height_interp(y,z) > rad_height)
                                rad_dbz = linear_interpolate(model_height_interp(y,max(1,z-1)),model_height_interp(y,z),rad_height,data_interp(y,max(1,z-1)),data_interp(y,z));
                                refl_value_seen(rad_idx,y) = rad_dbz;
                                break;
                            end
                        end
                        if(rad_dbz > max_dbz)
                            max_dbz = rad_dbz;
                            max_dbz_idx = rad_idx;
                        end
                    end
                    height_values(y) = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,radar_height_grid(max_dbz_idx,y,idx_a),radar_height_grid(max_dbz_idx,y,idx_b));
                    pressure_found = 0;
                    for z = 1:zlen_nas_clean
                        slice_values(y,z) = linear_interpolate(lon_nas_clean(y,idx_a),lon_nas_clean(y,idx_b),slice_lon,data(y,idx_a,z),data(y,idx_b,z));
                        if(model_height_interp(y,z) >= height_values(y) && ~pressure_found)
                            intersection_pressure(y) = linear_interpolate(model_height_interp(y,max(z-1,1)),model_height_interp(y,z),height_values(y),pressure_slice_wrf(y,max(z-1,1)),pressure_slice_wrf(y,z));
                            pressure_found = 1;
                        end
                    end
                end
            else % Slice is along x-axis
                height_values = zeros(xlen_nas_clean,1);
                intersection_pressure = zeros(xlen_nas_clean,1);
                intersection_pressure(:) = NaN;
                for x = 1:xlen_nas_clean
                    max_dbz_idx = 0;
                    max_dbz = -999;
                    for rad_idx = 1:num_radars
                        rad_height = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,radar_height_grid(rad_idx,idx_a,x),radar_height_grid(rad_idx,idx_b,x));
                        for z = 1:zlen_nas_clean
                            if(model_height_interp(x,z) > rad_height)
                                rad_dbz = linear_interpolate(model_height_interp(x,max(1,z-1)),model_height_interp(x,z),rad_height,data_interp(x,max(1,z-1)),data_interp(x,z));
                                refl_value_seen(rad_idx,x) = rad_dbz;
                                break;
                            end
                        end
                        if(rad_dbz > max_dbz)
                            max_dbz = rad_dbz;
                            max_dbz_idx = rad_idx;
                        end
                    end
                    height_values(x) = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,radar_height_grid(max_dbz_idx,idx_a,x),radar_height_grid(max_dbz_idx,idx_b,x));
                    pressure_found = 0;
                    for z = 1:zlen_nas_clean
                        slice_values(x,z) = linear_interpolate(lat_nas_clean(idx_a,x),lat_nas_clean(idx_b,x),slice_lat,data(idx_a,x,z),data(idx_b,x,z));
                        if(model_height_interp(x,z) >= height_values(x) && ~pressure_found)
                            intersection_pressure(x) = linear_interpolate(model_height_interp(x,max(z-1,1)),model_height_interp(x,z),height_values(x),pressure_slice_wrf(x,max(z-1,1)),pressure_slice_wrf(x,z));
                            pressure_found = 1;
                        end
                    end
                end
            end
        %end

        %% Plot vertical slices

        %beam_color = "#757575";
        %beam_color = "#0076ff";
        %beam_color = "#ff00f3";
        beam_color = "#000000";
        multibeam_color = "#ff00f3";

        plot_type = 'vsl';

        if(plot_all_beams)
            plot_type = strcat(plot_type,'-all');
        end
        if(plot_beam_line)
            plot_type = strcat(plot_type,'-base');
        end

        data_to_plot = slice_values;
        cmap = 'reflmap';
        output_format_slice = '%s/%s_%s_d0%d_sbd%d_%s_%s_%s_%s_s%d_%s.png'; % path, exp name a, exp name b, domain, subdomain, member string, bao, plot type, data type, slice idx, timestamp
        plot_filename = sprintf(output_format_slice,output_path_large,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,data_type,slice_idx,timestamp_file);
        plot_filename_small = sprintf(output_format_slice,output_path_small,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,data_type,slice_idx,timestamp_file);

        if(isfile(plot_filename) && ~remake_plots) % If shouldn't override existing plots, don't
            %continue;
        end

        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width_slice fig_height_slice]); % Create initial blank figure

        if(show_progress && announce_plots)
            toc
            fprintf('Plotting %s...\n',plot_type); % debug statement
        end


        if(y_axis)
            h = pcolor(slice_lat_grid,pressure_slice_wrf,data_to_plot); % Plot the data
        else
            h = pcolor(slice_lon_grid,pressure_slice_wrf,data_to_plot); % Plot the data
        end

        set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
        shading interp; % Smooth out plot from grid-boxes
        c = colorbar('FontSize',axes_font_size); % Make colorbar
        colormap(cmap); % Set colors
        ylabel(c,'Radar Reflectivity (dBZ)');
        caxis([clim_lower clim_upper]);
        set (gca,'YDir','reverse')

        hold on;

        if(plot_all_beams)
            if(y_axis)
                for rad_idx = 1:num_radars
                    h = plot(slice_lat_list,intersection_pressure_mult(rad_idx,:,:),'Color',multibeam_color,'LineWidth',2,'LineStyle','--');
                    %h = scatter(slice_lat_list,intersection_pressure_mult(rad_idx,:,:),5,'MarkerEdgeColor',multibeam_color,'MarkerFaceColor',multibeam_color);
                end
            else
                for rad_idx = 1:num_radars
                    h = plot(slice_lon_list,intersection_pressure_mult(rad_idx,:,:),'Color',multibeam_color,'LineWidth',2,'LineStyle','--');
                end
            end
        end
        if(plot_beam_line)
            if(y_axis)
                h = plot(slice_lat_list,intersection_pressure,'Color',beam_color,'LineWidth',2,'LineStyle','-');
            else
                h = plot(slice_lon_list,intersection_pressure,'Color',beam_color,'LineWidth',2,'LineStyle','-');
            end
        end
        if(limit_slice)
            if(y_axis)
                xlim([plot_s_slice_edge plot_n_slice_edge]);
                ylim([pressure_min pressure_max]);
            else
                xlim([plot_w_slice_edge plot_e_slice_edge]);
                ylim([pressure_min pressure_max]);
            end
        end

        % Apply labels

        if(y_axis)
            title(sprintf('[%s|d0%d|%s] Refl VSlice #%d (Lon %02.01f) | %s',exp_name_plot,domain,member_string,slice_idx,slice_lon,timestamp_title),'FontSize',title_font_size); 
            xlabel('Latitude','FontSize',label_font_size);
            ylabel('Pressure (mb)','FontSize',label_font_size);
        else
            title(sprintf('[%s|d0%d|%s] Refl VSlice #%d (Lat %02.01f) | %s',exp_name_plot,domain,member_string,slice_idx,slice_lat,timestamp_title),'FontSize',title_font_size); 
            xlabel('Longitude','FontSize',label_font_size);
            ylabel('Pressure (mb)','FontSize',label_font_size);
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
        title(sprintf('[%s|%s] Refl VSL #%d',exp_name_plot,member_string,slice_idx),'FontSize',title_font_size); 

        saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
        close('all');

        clearvars f h data_to_plot;

    end

end

% Get total runtime and print to stdout
runtime = toc;
hours = floor(runtime/3600);
mins = floor((runtime/60) - (hours*60));
secs = toc - (hours*3600) - (mins*60);
fprintf('Done. Total script runtime = %02.f:%02.f:%02.f\n',hours,mins,secs)

% END      