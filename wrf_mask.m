%% WRF_mask.m
% 
% Purpose: to create a binary mask on the WRF and GIS grids with which to
% cut off extraneous noncomparable data
% 
% Author(s): Jon Seibert
% Last updated: 15 Sept 2024
% 
% Inputs: [WRF Outputs].nc
% Outputs: mask_wrf_on_wrf.mat, mask_wrf_on_gis.mat
% Dependencies: borders.m, latlon_[]_trim.mat (predefined domain boundaries)
%
% NOTES:
%  - Not fully up to date with current formatting
%
% TODO:
%  - 
script_name = 'WRF_mask.m';

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

%% 0. General Settings

% Date & time (for sampling purposes only)
year = 2020;
month = 2;
day = 7;
hour = 14;
minute = 0;

% Spatial limits of analysis and plot (degrees lat/lon)
w_lim = -79;
e_lim = -69.75;
s_lim = 36;
n_lim = 46;

limit_borders = false; % Whether to apply spatial limits to plot
trim_e = false; % Whether to trim the eastern edge

% Figure font sizes
title_font_size = 18;
label_font_size = 14;
axes_font_size = 12;

data_type = "refl"; % Variable being analyzed
units = "dBZ";

run_name = 'mask';

%% 1A. Settings- WRF Data (2022 Case)

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020';
    input_path_gis = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    input_path_nas = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/exrad_data/';
    intermediate_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/ptp';
    path_to_code = "C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code";
    path_to_extra_code = './downloaded_code';
    input_path_ecmwf = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Data/ECMWF';
    input_path_wrf = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/code/input/wrf_data/2020/AIR/202002071400/';
    input_path_gis = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/code/input/gis_data/2020';
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/';
    input_path_gis = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/GIS';
    intermediate_path = '/storage/home/jjs5895/projects/IMPACTS/intermediate';
    output_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/ptp';
    path_to_code = "/storage/home/jjs5895/projects/IMPACTS/code";
    path_to_extra_code = '/storage/home/jjs5895/projects/IMPACTS/code/downloaded_code'; % Specify path to borders.m
    input_path_ecmwf = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/ECMWF';
    input_path_wrf = '-';
    input_path_gis = '-';
end

addpath(path_to_extra_code);

% Basic details
domain = 2;
member = 1; % Only one is needed
data_src_wrf = "wrf"; % Dataset label

% NetCDF retrieval
lon_name_wrf = "XLONG"; 
lat_name_wrf = "XLAT"; 
data_name_wrf = "REFL_10CM"; % Variable being analyzed

% Note whether transposition or compositing is necessary
transpose_data = true;
flip_data = true;
composite_data = true;

filename_format_wrf = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #

%% 1B. Settings- GIS Data

data_src_gis = "gis";

file_format_gis = '%s/n0q_%04.f%02.f%02.f%02.f%02.f.png';

lat_gis_raw = 49.9975:-0.005:23.0025;
lon_gis_raw = -125.9975:0.005:-65.0025;
[lon_gis,lat_gis] = meshgrid(lon_gis_raw,lat_gis_raw);

% NOTES: Comes in image format

%%

% N0R product is angle 0.5, 16 levels/ 230 km (short range base refl.) N0Q
% is long range: 460 km, still 0.5. we're using N0Q. 8 bit Base Reflectivity 0.5dbz
beam_angle = 0.5;
max_dist = 230000; % 230 km

g = 9.8; % m/s^2;
Re = 6371000; % radius of Earth in m
        
%% 2. Preliminary Processing & Setup

output_format = '%s/wrf_%s_d0%d_%s_%s_%s_%s_%s.png'; % output path, data source 2, domain, mem/mean, bao_short, fit/dif, data type, timestamp
datetime_format_file = '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

% Define sample file to read from [deprecated]
timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
filename_wrf = sprintf(filename_format_wrf,domain,member);
filename_gis = sprintf(file_format_gis,year,month,day,hour,minute);

% Retrieve lat/lon data
load(sprintf('%s/%s',intermediate_path,'nexrad_stations.mat')); % wban, station_ids,_station_names,lon_deg,lat_deg,elevations,tower_heights
load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'lon_wrf_trim','lat_wrf_trim'); % Load in preestablished latlon values

lon_wrf = lon_wrf_trim;
lat_wrf = lat_wrf_trim;

lon_station = lon_deg;
lat_station = lat_deg;

dimensions_wrf = size(lon_wrf);

%% Read in GIS data
% (Not entirely sure this section is still needed. Ain't broke don't fix.)

% GIS lon/lat are predetermined- image files do not contain this data

% Find shared area between WRF and GIS latlon coverage
w_bound = max(min(min(lon_wrf)),min(min(lon_gis)));
e_bound = min(max(max(lon_wrf)),max(max(lon_gis)));
s_bound = max(min(min(lat_wrf)),min(min(lat_gis)));
n_bound = min(max(max(lat_wrf)),max(max(lat_gis)));

% IN THIS CASE, GIS ENTIRELY CONTAINS WRF: FIND EDGE POINTS IN GIS
% LATLON THAT FALL JUST INSIDE THE WRF BOUNDARIES
% Trim GIS to those points, regrid WRF to the trimmed GIS

for x=1:length(lon_gis_raw)
    if(lon_gis_raw(x) > w_bound && lon_gis_raw(x) < e_bound)
        w_idx = x;
        break;
    end
end

for x=length(lon_gis_raw):-1:1
    if(lon_gis_raw(x) > w_bound && lon_gis_raw(x) < e_bound)
        e_idx = x;
        break;
    end
end

for y = 1:length(lat_gis_raw)
    if(lat_gis_raw(y) > s_bound && lat_gis_raw(y) < n_bound)
        s_idx = y;
        break;
    end
end

for y = length(lat_gis_raw):-1:1
    if(lat_gis_raw(y) > s_bound && lat_gis_raw(y) < n_bound)
        n_idx = y;
        break;
    end
end

lon_gis_trim_raw = lon_gis_raw(w_idx:1:e_idx);
lat_gis_trim_raw = lat_gis_raw(s_idx:1:n_idx);

[lon_gis_trim,lat_gis_trim] = meshgrid(lon_gis_trim_raw,lat_gis_trim_raw);

% Read in
data_gis = double(imread(sprintf(file_format_gis,input_path_gis,year,month,day,hour,minute)))*0.5 - 32.5; 
data_gis_trim = data_gis(s_idx:1:n_idx,w_idx:1:e_idx); % Trim down to WRF domain
dimensions_gis_trim = size(data_gis_trim);

%%

% Define ocean exclusion line
eel_x1 = -76;
eel_x2 = -69.8725;
eel_y1 = 36;
eel_y2 = 42.3825;

% y = m*x + b;
% x = (y-b)/m;
% b = y - m*x;
m = (eel_y2 - eel_y1)/(eel_x2 - eel_x1);
b = eel_y1 - m*eel_x1;


%% Take edges of WRF latlon grid at each edge and define mask on GIS data

wrf_gis_mask = zeros(dimensions_gis_trim);

for y = 1:dimensions_gis_trim(1)
    for x = 1:dimensions_gis_trim(2)
        if(lon_gis_trim(y,x) < e_lim && lon_gis_trim(y,x) > w_lim && lat_gis_trim(y,x) < n_lim && lat_gis_trim(y,x) > s_lim)
            if(lon_gis_trim(y,x) < ((lat_gis_trim(y,x)-b)/m))
                for station_idx = 1:length(station_ids)
                    d = haversine_distance(lon_station(station_idx),lat_station(station_idx),lon_gis_trim(y,x),lat_gis_trim(y,x));
                    if(d <= max_dist)
                        wrf_gis_mask(y,x) = 1;
                        break;
                    end
                end
            end
        end
    end
end

save(sprintf('%s/mask_wrf_on_gis.mat',intermediate_path),'wrf_gis_mask');

%% Take edges of WRF latlon grid at each edge and define mask on WRF data

wrf_wrf_mask = zeros(dimensions_wrf);

for y = 1:dimensions_wrf(1)
    for x = 1:dimensions_wrf(2)
        if(lon_wrf(y,x) < e_lim && lon_wrf(y,x) > w_lim && lat_wrf(y,x) < n_lim && lat_wrf(y,x) > s_lim)
            if(lon_wrf(y,x) < ((lat_wrf(y,x)-b)/m))
                for station_idx = 1:length(station_ids)
                    d = haversine_distance(lon_station(station_idx),lat_station(station_idx),lon_wrf(y,x),lat_wrf(y,x));
                    if(d <= max_dist)
                        wrf_wrf_mask(y,x) = 1;
                        break;
                    end
                end
            end
        end
    end
end

save(sprintf('%s/mask_wrf_on_wrf.mat',intermediate_path),'wrf_wrf_mask');

%% Show mask [WIP]

% Plot LARGE
f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
h = pcolor(lon_gis_trim,lat_gis_trim,wrf_gis_mask); % Plot the data
set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
shading interp; % Smooth out plot from grid-boxes
%colorbar('FontSize',axes_font_size); % Make colorbar
%colormap(cmap); % Set colors


% Plot state borders
hold on;
borders('continental us','black','linewidth',1); 
hold off;

% Focus on desired area, remove whitespace
if(limit_borders)
    xlim([w_lim e_lim]);
    ylim([s_lim n_lim]);
end

% Apply labels

xlabel('Longitude','FontSize',label_font_size);
ylabel('Latitude','FontSize',label_font_size);
set(gca,'Fontsize',axes_font_size);




