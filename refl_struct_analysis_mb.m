% refl_struct_analysis (v2.1.2)
%
% Purpose: Generate 2D map of the vertical heights / pressure levels of the
% maximum radar reflectivity per column & vertical slices
% MB VERSION: WORKS IN PRESSURE COORDINATES FOR VERTICAL
% 
% Author(s): Jon Seibert
% Last updated: 15 Sept 2024
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
script_name = 'refl_struct_analysis.m';

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

run_name = 'vslice_8_v3';

% Local or server execution mode
%local = 1; % true = running on local device, false = running on server

% Which experiments to use (0 is obs, only for B)
% To compare background and analysis, use the same experiment # for both
% To compare to observations, use 0 for experiment B
exp_choice_a = 1; 
exp_choice_b = 2;

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

if(local)
    num_members = 1;
    %num_members = 40; % TESTING
else
    num_members = 40;
    %num_members = 1; % TESTING
end

domain = 2;
subdomain = 0; %1 for E Coast Only subdomain, 2 for N only, 0 for all d02 included

use_contours = 0;
smooth_data = 0; % If true AND the variable isn't already being smoothed, smooth it before plotting (gaussian)

% Vertical levels to use if is 3D
%pressure_levels = [460,500,700,850,900,925];
pressure_target = 0; %mb [Artifact of copied code from ptp_analysis, still needed for filename architecture

multi_height = 0;
pressure_target_list = [925,850,700]; % Overrides above value

% REFL_10CM/REFL (dBZ), T/T (deg C), P/MSLP (mb), GPH/GPH(m),
% VORT/VORT(10^-5 s^-1), U-V-W/WIND(kt)
is_3d = 1;
%use_base_ref = 0; % Not used in this script

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
pressure_max = 1025; % NO LONGER OVERRIDDEN

% Plot selections
remake_plots = 1; % If true, overrides existing outputs. If false, only plots 'new' files.

show_plots = 0;
show_progress = 1; % If true, prints out progress bar/debug messages
announce_plots = 0;
use_auto_run_name = 0;

gauss_sd = 2.75;
gauss_width =  2*ceil(2*gauss_sd)+1;

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
%plot_beam_width = 1;
plot_skewt = 1;

radar_heights_filename = 'NEXRAD_radar_heights_wrf_grid_trim.mat';
%wrf_heights_filename = 'wrf_base_refl_heights.mat';

max_dist = 230000; %(m): radar radius

st_name_list = ["KENX","KBUF","KGYX","W-Mobile"]; % KENX = KALY == ALB
st_lat_list = [0,0,0,43.103];
st_lon_list = [0,0,0,-76.193];
num_st_official = 3;
num_st = 4;

%skewt_lat = 43.103;
%skewt_lon = -76.193;
make_wrf_soundings = 0;
make_height_maps = 0;
test_lv_eqn = 0;

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
    input_path_st = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Data/IMPACTS_field_Catalog';
    % Experiment paths
    exp_path_1 = 'AIR'; % TESTING
    exp_path_2 = 'CONV/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/';
    input_path_gis = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/GIS';
    intermediate_path = '/storage/home/jjs5895/projects/IMPACTS/intermediate';
    output_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/ptp';
    path_to_code = "/storage/home/jjs5895/projects/IMPACTS/code";
    path_to_extra_code = '/storage/home/jjs5895/projects/IMPACTS/code/downloaded_code'; % Specify path to borders.m
    input_path_ecmwf = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/ECMWF';
    input_path_st = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/IMPACTS_field_catalog';
    % Experiment paths
    exp_path_1 = 'AIRCFT/fc';
    exp_path_2 = 'CONV/fc/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
end

addpath(path_to_extra_code);

% Experiment names
exp_name_1 = 'AIR';
exp_name_2 = 'CONV';
exp_name_3 = 'AIR-F';
exp_name_4 = 'CONV-F';
exp_name_5 = 'NODA';

alt_input_exps = [3 4 5]; % Input is set up differently for these experiments

hour_step = 1; % Hour increment size

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 925;
fig_height = 900;
fig_width_small = 350;
fig_height_small = 350;

fig_width_slice = 1225;
fig_height_slice = 500;
fig_width_slice_small = 600;
fig_height_slice_small = 300;

fig_width_txrh = 500;
fig_height_txrh = 500;

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

% Subdomain settings
bound_lon_d1 = [-78,-78,-74.5,e_lim,e_lim,-75,-78];
bound_lat_d1 = [40,44,n_lim,n_lim,41.5,38,40];
%inpolygon([-74,50],[42,50],bound_lon,bound_lat);
bound_lon_d2 = [-77.5,-77.5,-74,e_lim,e_lim,-73.5,-77.5];
bound_lat_d2 = [42,44,n_lim,n_lim,42.5,42,42];

%% 1A. Settings- WRF Data (2022 Case)

% Basic details

%num_members = 40;
%num_members = 2; % TESTING
data_src_wrf = "wrf"; % Dataset label

% NetCDF retrieval
lon_name_wrf = "XLONG"; 
lat_name_wrf = "XLAT"; 

time_name_wrf = "Times";

% wrf_enkf_output_d02_001
% wrfinput_d02_2020-02-07_14:00:00_017
    % Model forward step from 13->14 is in 1300 folder, labeled with 1400.
    % DA Analysis step from 14b->14a is in 1400 folder, labeled with only domain and member #
    
% For time t:
%   Background = wrfinput_d02_t (found in t-1 directory)
%   Analysis = wrf_enkf_output_d02 (found in t directory)
% For free runs:
%   wrfout_d01_2020-02-07_20:00:00
filename_format_analysis = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #
filename_format_background = 'wrfinput_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00_%03.f'; % domain, year, month, day, hour, member #
filename_format_background_alt = 'wrfinput_d%02.f_%04.f-%02.f-%02.f_%02.f-00-00_%03.f'; % domain, year, month, day, hour, member #
filename_format_free = 'wrfout_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00'; % domain, year, month, day, hour

PT_varnames = ["T","PB","P"];
data_name_temp = "T"; % PERTURBATION POTENTIAL TEMPERATURE
data_name_pres = "PB"; % Base Pressure 
data_name_pres_prime = "P"; % Perturbation pressure 
data_name_gp = "PHB"; % Base Geopotential (m^2/s^2) [Z Stagger]
data_name_gp_prime = "PH"; % Perturbation Geopotential (m^2/s^2) [Z Stagger]
data_name_terrain_height = "HGT"; % Terrain height above sea level (m)
data_name_t2m = "T2"; % 2m Temperature
data_name_mr = "QVAPOR"; % Water vapor mixing ratio, kg/kg [NOT specific humidity]: mr, not Q
data_name_u = "U"; % Wind [X Stagger]
data_name_v = "V"; % Wind [Y Stagger]
data_name_w = "W"; % Wind [Z Stagger]

Re = 6.3781e6; % Radius pf Earth in m
g = 9.80665; % m/s^2 standard gravity at sea level
R = 287.0600676; % J/kg/K, specific gas constant

%% 1C. Extra Model information

% ECMWF reanalysis filepath

ecmwf_filename = 'ecmwf_mslp_d01_20200207.nc';
       
beam_width = 0.9; % Degrees, DIAMETER
beam_angle = 0.5;
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
title_format_fit =   '[%s-%s|d0%d|%s] %s: %s(%s)|%s'; % exp name 1, exp name 2, domain, member #/mean, background/analysis/obs (bao), bao, plot_type, data type, units, timestamp
title_format_dif =   '[%s-%s|d0%d|%s] %s: %s(%s) [%smb|%s]\nRMSE = %0.3f, Bias = %0.3f'; % exp name 1, exp name 2, domain, member #/mean, plot_type, data type, units, height timestamp, err, bias
title_format_trim =  '[%s] %s: %s(%s)|%s'; % data source, plot_type, data type, units, timestamp
title_format_raw =   '[%s|d0%d|%s] %s(%s) [%smb|%s]'; % exp name, domain, member #/mean, data type, units, height timestamp
title_format_inc =   '[%s|d0%d|%s] Analysis Inc: %s(%s) [%smb|%s]'; % exp name, domain, member #/mean, data type, units, height timestamp
title_format_sd =    '[%s|d0%d] %s: %s(%s)|%s'; % exp name, domain, sd_type, data type, units, timestamp
title_format_small = '[%s-%s] %s-%s-%s'; % exp name 1, exp name 2,plot_type,data type,timestamp
output_format =         '%s/%s_%s_d0%d_sbd%d_%s_%s_%s_%s_%s.png'; % output path, exp name 1, exp name 2, domain, subdomain, mem/mean, bao_short, plot_type, data type, timestamp
output_format_var =     '%s/%s_%s_d0%d_sbd%d_%s_%s_%s_%s_%dmb_%s.png'; % output path, exp name 1, exp name 2, domain, subdomain, mem/mean, bao_short, plot_type, data type, height, timestamp
output_format_small =   '%s/%s_%s_d0%d_sbd%d_%s_%s_%s_%s_%s_small.png'; % output path, exp name 1, exp name 2, domain, subdomain, mem/mean, bao_short, plot_type, data type, timestamp
output_format_var_small='%s/%s_%s_d0%d_sbd%d_%s_%s_%s_%s_p%d_%s_small.png'; % output path, exp name 1, exp name 2, domain, subdomain, mem/mean, bao_short, plot_type, data type, timestamp
input_format_base_refl ='%s/%s_%s_d%02.f_%s_%s_%s.mat'; % output path, data source, plot type, domain, member #, data type, datetime
datetime_format_file =  '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]
error_table_format = '%s/ev_table_%s_%dmb_%s.txt';

% Define run names and flags
bao_short = 'a';
bao_1 = 'Analysis';
compare_to_obs = false;

switch exp_choice_a
    case 1
        input_path_a = sprintf('%s/%s',input_path_base,exp_path_1);
        exp_name_a = exp_name_1;
    case 2
        input_path_a = sprintf('%s/%s',input_path_base,exp_path_2);
        exp_name_a = exp_name_2;
    case 3
        input_path_a = sprintf('%s/%s',input_path_base,exp_path_3);
        exp_name_a = exp_name_3;
    case 4
        input_path_a = sprintf('%s/%s',input_path_base,exp_path_4);
        exp_name_a = exp_name_4;
    case 5
        input_path_a = sprintf('%s/%s',input_path_base,exp_path_5);
        exp_name_a = exp_name_5;
    otherwise
        error('Unrecognized experiment choice');
end


if(exp_choice_a == exp_choice_b) % Note: Pointless to compare obs to obs
    bao_2 = 'Background';
    bao_add = 'b';
    compare_to_background = true;
else
    bao_2 = 'Analysis';
    bao_add = 'a';
    compare_to_background = false;
end
switch exp_choice_b
    case 0
        compare_to_obs = true;
        bao_2 = 'Obs';
        bao_short = strcat(bao_short,'o');
        input_path_b = input_path_gis;
        exp_name_b = 'OBS';
    case 1
        bao_short = strcat(bao_short,bao_add);
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_1);
        exp_name_b = exp_name_1;
    case 2
        bao_short = strcat(bao_short,bao_add);
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_2);
        exp_name_b = exp_name_2;
    case 3
        bao_short = strcat(bao_short,bao_add);
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_3);
        exp_name_b = exp_name_3;
    case 4
        bao_short = strcat(bao_short,bao_add);
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_4);
        exp_name_b = exp_name_4;
    case 5
        bao_short = strcat(bao_short,bao_add);
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_5);
        exp_name_b = exp_name_5;
    otherwise
        error('Unrecognized experiment choice');
end

if(compare_to_background)
    exp_name_a_plot = sprintf('%s-A',exp_name_a);
    exp_name_b_plot = sprintf('%s-B',exp_name_b);
else
    exp_name_a_plot = exp_name_a;
    exp_name_b_plot = exp_name_b;
end

if(use_auto_run_name)
    run_name = sprintf('%s_%s_%s',lower(exp_name_a),lower(exp_name_b),bao_short);
end
fprintf('Run designation: %s.\n',run_name);
if(show_progress)
    toc
    if(remake_plots)
        fprintf('WARNING: Overriding existing designated plots.\n');
    end
    fprintf('Performing initial setup.\n');
end

max_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 30, 31, 31]; % Number of days in each month

% Time counters [NOTE: NOT ROBUST TO CROSSING MONTH BOUNDARIES! Hard-code if needed]
time_count = (end_day - start_day + 1)*24 - (start_hour) - (24 - end_hour) + 1;

year = start_year;
month = start_month;
day = start_day;
hour = start_hour;
minute = 0;
member_idx = 1;

% Define sample file to read latlon data from
timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
filename_wrf = sprintf(filename_format_analysis,domain,member_idx);

% Create output folders
% Directory structure:
% output/ptp/(run_name)/(variable_name)/(plot_size)
output_path_large = sprintf('%s/%s/%s/%s/%s',output_path_base,run_name,data_type,'large');
output_path_small = sprintf('%s/%s/%s/%s/%s',output_path_base,run_name,data_type,'small');
output_path_st = sprintf('%s/%s/%s/%s',output_path_base,run_name,"ST");

% If output folders do not exist, create them
if(~isfolder(output_path_large) || ~isfolder(output_path_small))
    mkdir(output_path_large);
    mkdir(output_path_small);
    mkdir(sprintf('%s/%s/%s/%s',output_path_base,run_name,"ST/custom"));
end

% REFL_10CM/REFL (dBZ), T/T (deg C), P/MSLP (mb), GPH/GPH(m),
% VORT/VORT(10^-5 s^-1), U-V-W/WIND(kt)
% Define units
contour_steps = 6; % P: 4, T: 2, GPH: 6 (default to 6 for GPH)
plot_fl_np = false;
switch data_type
    case "GPH"
        units = "dam";
        data_name = ""; % WRF Referent
        contour_steps = 6; % P: 4, T: 2, GPH: 6
        clim_lower_dif = -4;
        clim_upper_dif = 4;
        use_contours = 1;
        limit_raw_colorbar = true;
    case "T"
        units = sprintf('%sC',char(176));
        data_name = "T";
        contour_steps = 2; % P: 4, T: 2, GPH: 6
        use_contours = 1;
        plot_fl_np = 1;
    case "MSLP"
        units = "mb";
        data_name = "P"; 
        contour_steps = 4; % P: 4, T: 2, GPH: 6
        use_contours = 1;
    case "REFL"
        units = "dBZ";
        data_name = "REFL_10CM"; 
        contour_steps = 2; %mb
        contour_units = "mb";
    case "WIND"
        units = "kt";
        data_name = "";
        contour_steps = 6; % P: 4, T: 2, GPH: 6 (plots with GPH)
        use_contours = 0;
    case "VORT"
        %units = "1/s";
        units = "s^{-1}";
        data_name = "";
        use_contours = 0;
        limit_raw_colorbar = true;
    case "Q"
        %units = "kg/kg";
        units = "kg kg^{-1}";
        data_name = data_name_mr;
        clim_lower = 0;
        clim_upper = 5e-3;
        clim_lower_dif = -2e-3;
        clim_upper_dif = 2e-3;
        use_contours = 0;
        limit_raw_colorbar = true;
    case "OMEGA"
        %units = "mb/s";
        units = "hPa s^{-1}";
        data_name = "";
        clim_lower = -0.3;
        clim_upper = 0.3;
        clim_lower_dif = -0.3;
        clim_upper_dif = 0.3;
        use_contours = 0;
        limit_raw_colorbar = true;
end

% Centered finite difference distances (cell 1 -> cell 3 across cell 2)
if(domain == 2) % 3km res
    dx = 6000;
    dy = 6000;
    load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'lon_wrf_trim','lat_wrf_trim','n_idx','s_idx','e_idx','w_idx'); % Load in preestablished latlon values
    lon_wrf = lon_wrf_trim;
    lat_wrf = lat_wrf_trim;
    clim_lower_vort = 0;
    clim_upper_vort = 100; % 10e-5
    
elseif(domain == 1) % 9km res
    dx = 18000;
    dy = 18000;
    w_lim = w_lim_d01;
    e_lim = e_lim_d01;
    s_lim = s_lim_d01;
    n_lim = n_lim_d01;
    fig_x = fig_x_d01;
    fig_y = fig_y_d01;
    fig_width = fig_width_d01;
    fig_height = fig_height_d01;
    fig_width_small = fig_width_d01_small;
    fig_height_small = fig_height_d01_small;
    clim_lower_vort = 0;
    clim_upper_vort = 50;
    
    load(sprintf('%s/latlon_wrf_trim_d01.mat',intermediate_path));
    lon_wrf = lon_wrf_trim;
    lat_wrf = lat_wrf_trim;
    
end

dimensions_wrf = size(lon_wrf);
dims = size(lon_wrf);
ylen_wrf = dimensions_wrf(1);
xlen_wrf = dimensions_wrf(2);
zlen_wrf = 50;
clearvars dimensions_wrf lon_wrf_trim lat_wrf_trim;

% Create cross-member storage objects
%data_storage_a = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
%data_storage_b = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);

% Existing model height field name is "model_height" (above ground)
% Load NEXRAD beam heights and WRF Base Refl height map
load(sprintf('%s/%s',intermediate_path,radar_heights_filename),'radar_height_grid');
radar_height_grid_full = radar_height_grid;
load(sprintf('%s/%s',intermediate_path,'nexrad_stations.mat')); % wban, station_ids,station_names,lon_deg,lat_deg,elevations,tower_heights

for st_idx = 1:num_st_official
    found_idx = find(station_ids == st_name_list(st_idx));
    st_lon_list(st_idx) = lon_deg(found_idx);
    st_lat_list(st_idx) = lat_deg(found_idx);
end


%% MAIN

% Identify filepaths

input_path_base_refl_a = sprintf('%s/BASE_REF/%s',input_path_base,exp_name_a);
input_path_base_refl_b = sprintf('%s/BASE_REF/%s',input_path_base,exp_name_b);

load(sprintf('%s/%s',intermediate_path,'latlon_gis_trim.mat'),'lon_gis_trim','lat_gis_trim'); % Load in preestablished latlon values
gauss_filt = fspecial('gaussian',gauss_width,gauss_sd);

if(show_progress)
    toc 
    fprintf('Comparison: %s-%s, %s\n',exp_name_a,exp_name_b,data_type);
end

%% 3B. For each timestamp being analyzed:
for time_idx = 1:time_count
   
    % Set timestamps and experiment names
    timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
    timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
    timestamp_file_background = sprintf(datetime_format_file,year,month,day,hour-1);
    data_src_a = lower(exp_name_a);
    data_src_b = lower(exp_name_b);
    
    % Define holding objects
    ens_mem_a = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_b = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_t_a = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_t_b = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_p_a = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_p_b = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_mr_a = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_mr_b = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_mh_a = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_mh_b = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);

    % Read in each member and compute the mean
    for member_idx = 1:(num_members+1)

        % If all members have been computed, take mean
        if(member_idx == (num_members+1)) 
            member_string = 'mean';
            data_a = mean(ens_mem_a,4,"omitnan");
            data_b = mean(ens_mem_b,4,"omitnan");
            data_mr_a = mean(ens_mem_mr_a,4,"omitnan");
            data_mr_b = mean(ens_mem_mr_b,4,"omitnan");
            data_t_a = mean(ens_mem_t_a,4,"omitnan");
            data_t_b = mean(ens_mem_t_b,4,"omitnan");
            pressure_a = mean(ens_mem_p_a,4,"omitnan");
            pressure_b = mean(ens_mem_p_b,4,"omitnan");
            model_height_a = mean(ens_mem_mh_a,4,"omitnan");
            model_height_b = mean(ens_mem_mh_b,4,"omitnan");
            
            if(show_progress)
                toc
                fprintf('Processing: %s | Ensemble Mean\n',timestamp_title);
            end
        else % Otherwise, do the next member
        
            member_string = sprintf('%03.f',member_idx);
            
            if(show_progress)
                toc
                fprintf('Processing: %s | Member %s\n',timestamp_title,member_string);
            end
            
            % Define full filepaths- exps 3,4,5 have dif format
            if(ismember(exp_choice_a,alt_input_exps))
                filename_a = sprintf(filename_format_free,domain,year,month,day,hour);
                full_filepath_a = sprintf('%s/%s/%s',input_path_a,member_string,filename_a);
            else
                filename_a = sprintf(filename_format_analysis,domain,member_idx);       
                full_filepath_a = sprintf('%s/%s/%s',input_path_a,timestamp_file,filename_a);
            end
            if(ismember(exp_choice_b,alt_input_exps))
                filename_b = sprintf(filename_format_free,domain,year,month,day,hour);
                full_filepath_b = sprintf('%s/%s/%s',input_path_b,member_string,filename_b);
            elseif(compare_to_background)
                if(local)
                    filename_b = sprintf(filename_format_background_alt,domain,year,month,day,hour,member_idx); % Local files can't have : in filename
                else
                    filename_b = sprintf(filename_format_background,domain,year,month,day,hour,member_idx);     
                end
                full_filepath_b = sprintf('%s/%s/%s',input_path_b,timestamp_file_background,filename_b);
            else
                filename_b = sprintf(filename_format_analysis,domain,member_idx);       
                full_filepath_b = sprintf('%s/%s/%s',input_path_b,timestamp_file,filename_b);
            end
            
            lon_wrf_raw = ncread(full_filepath_a,lon_name_wrf);
            lat_wrf_raw = ncread(full_filepath_a,lat_name_wrf);
            dims = size(lon_wrf_raw);
            xlen_wrf_raw = dims(1);
            ylen_wrf_raw = dims(2);
            %zlen_wrf_raw = dims(3);

            % Read in all necessary WRF fields
            data_pb_a = ncread(full_filepath_a,data_name_pres); % Base pressure (Pa)
            data_pp_a = ncread(full_filepath_a,data_name_pres_prime); % Perturbation pressure (Pa)
            data_ptp_a = ncread(full_filepath_a,data_name_temp); % Perturbation potential temperature (K)
            data_p_a = data_pb_a + data_pp_a; % 3D Pressure field (Pa)
            data_pt_a = data_ptp_a + 300; % Potential temperature (K)
            data_t_a = data_pt_a.*((data_p_a./100000).^(0.287)); % 3D Temperature field (K)
            data_gpb_raw_a = ncread(full_filepath_a,data_name_gp); % Base geopotential (m^2/s^2)
            data_gpp_raw_a = ncread(full_filepath_a,data_name_gp_prime); % Perturbation geopotential (m^2/s^2)
            data_p_a = data_p_a./100; % Convert P from Pascals (Pa) to Millibars (mb)
            data_mr_a = ncread(full_filepath_a,data_name_mr); % Water vapor (kg/kg)

            data_pb_b = ncread(full_filepath_b,data_name_pres);
            data_pp_b = ncread(full_filepath_b,data_name_pres_prime); 
            data_ptp_b = ncread(full_filepath_b,data_name_temp); 
            data_p_b = data_pb_b + data_pp_b; 
            data_pt_b = data_ptp_b + 300; 
            data_t_b = data_pt_b.*((data_p_b./100000).^(0.287));
            data_gpb_raw_b = ncread(full_filepath_b,data_name_gp); 
            data_gpp_raw_b = ncread(full_filepath_b,data_name_gp_prime); 
            data_p_b = data_p_b./100; 
            data_mr_b = ncread(full_filepath_b,data_name_mr); % Water vapor (kg/kg)

            if(member_idx == 1)
                terrain_height = ncread(full_filepath_a,data_name_terrain_height); % Height of physical terrain (m)
            end

            % Define holding containers
            data_gpb_a = zeros(size(data_t_a));
            data_gpb_b = zeros(size(data_t_b));
            data_gpp_a = zeros(size(data_t_a));
            data_gpp_b = zeros(size(data_t_b));

            % Fix staggered WRF grids for wind, geopotential      
            for z = 1:zlen_wrf
                data_gpb_a(:,:,z) = (data_gpb_raw_a(:,:,z) + data_gpb_raw_a(:,:,z+1))./2;
                data_gpp_a(:,:,z) = (data_gpp_raw_a(:,:,z) + data_gpp_raw_a(:,:,z+1))./2;
                data_gpb_b(:,:,z) = (data_gpb_raw_b(:,:,z) + data_gpb_raw_b(:,:,z+1))./2;
                data_gpp_b(:,:,z) = (data_gpp_raw_b(:,:,z) + data_gpp_raw_b(:,:,z+1))./2;
            end

            % Compute geopotential height
            data_gp_a = data_gpb_a + data_gpp_a; %m^2/s^2
            data_gph_a = data_gp_a./g; %m
            data_gp_b = data_gpb_b + data_gpp_b;
            data_gph_b = data_gp_b./g;
            data_gph_a = data_gph_a/10; % Convert from m to dam (decameters)
            data_gph_b = data_gph_b/10;  

            model_height_a = ((data_gpb_a+data_gpp_a)./g) - terrain_height; % model height above the *ground* (m)
            model_height_b = ((data_gpb_b+data_gpp_b)./g) - terrain_height; % model height above the *ground* (m)

            % Align all data fields with customized grid
            pressure_a = align_data(data_p_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            pressure_b = align_data(data_p_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_mr_a = align_data(data_mr_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_mr_b = align_data(data_mr_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_t_a = align_data(data_t_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_t_b = align_data(data_t_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            model_height_a = align_data(model_height_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            model_height_b = align_data(model_height_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);

            data_a = ncread(full_filepath_a,data_name); % Retrieve data
            data_b = ncread(full_filepath_b,data_name); % Retrieve data
            
            dims = [ylen_wrf xlen_wrf zlen_wrf];
            
            % All data must be 2D to be processed and plotted. If not
            % composite or interpolated already:
            data_a = align_data(data_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_b = align_data(data_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);

            % Store values for next loop
            ens_mem_a(:,:,:,member_idx) = data_a;
            ens_mem_b(:,:,:,member_idx) = data_b;
            ens_mem_t_a(:,:,:,member_idx) = data_t_a;
            ens_mem_t_b(:,:,:,member_idx) = data_t_b;
            ens_mem_p_a(:,:,:,member_idx) = pressure_a;
            ens_mem_p_b(:,:,:,member_idx) = pressure_b;
            ens_mem_mr_a(:,:,:,member_idx) = data_mr_a;
            ens_mem_mr_b(:,:,:,member_idx) = data_mr_b;
            ens_mem_mh_a(:,:,:,member_idx) = model_height_a;
            ens_mem_mh_b(:,:,:,member_idx) = model_height_b;

        end
        
        
        %% Compute height map
        if(make_height_maps)
            height_map_a = zeros(ylen_wrf,xlen_wrf);
            height_map_b = zeros(ylen_wrf,xlen_wrf);
            
            [~,max_idx] = max(data_a,[],3);
            for y = 1:ylen_wrf
                for x = 1:xlen_wrf
                    height_map_a(y,x) = pressure_a(y,x,max_idx(y,x));
                end
            end
        
            [~,max_idx] = max(data_b,[],3);
            for y = 1:ylen_wrf
                for x = 1:xlen_wrf
                    height_map_b(y,x) = pressure_b(y,x,max_idx(y,x));
                end
            end   
        end

        %% Plot vertical slices: x and y-aligned only
        
        for exp_idx = 1:2
            switch exp_idx
                case 1
                    exp_name_plot = exp_name_a;
                    data = data_a;
                    model_height = model_height_a;
                case 2
                    exp_name_plot = exp_name_b;
                    data = data_b;
                    model_height = model_height_b;
            end
            for slice_idx = 1:3
                switch slice_idx
                    case 1
                        % Slice 1: Y-Z, along latitude:
                        slice_lon = -74.4;
                        y_axis = 1;
                    case 2
                        % Slice 2: X-Z, along longitude:
                        slice_lat = 42.5;
                        y_axis = 0;
                    case 3
                        % Slice 3: X-Z, along longitude:
                        slice_lat = 44;
                        y_axis = 0;
                    otherwise
                end

                if(y_axis) % Slice is along y-axis
                    slice_lon_list = ones(xlen_wrf,1).*slice_lon;
                    slice_lat_list = zeros(ylen_wrf,1);
                    slice_values = zeros(ylen_wrf,zlen_wrf);
                    idx_a = 0; idx_b = 0;
                    y_mid = floor(ylen_wrf/2);
                    for x = 1:xlen_wrf
                        if(lon_wrf(y_mid,x) == slice_lon)
                            idx_a = x; idx_b = x;
                            break;
                        elseif(lon_wrf(y_mid,x) > slice_lon)
                            idx_a = max(1,x-1); idx_b = x;
                            break;
                        end
                    end
                    if(exp_idx == 1)
                        pressure_slice = squeeze(pressure_a(:,idx_a,:));
                        %pressure_max = max(pressure_a,[],"all");
                    else
                        pressure_slice = squeeze(pressure_b(:,idx_a,:));
                        %pressure_max = max(pressure_b,[],"all");
                    end
                    for y = 1:ylen_wrf
                        slice_lat_list(y) = (lat_wrf(y,idx_a)+lat_wrf(y,idx_b))/2;
                    end
                else % Slice is along x-axis
                    slice_lon_list = zeros(xlen_wrf,1);
                    slice_lat_list = ones(ylen_wrf,1).*slice_lat;
                    slice_values = zeros(xlen_wrf,zlen_wrf);
                    
                    idx_a = 0; idx_b = 0;
                    x_mid = floor(xlen_wrf/2);
                    for y = ylen_wrf:-1:1
                        if(lat_wrf(y,x_mid) == slice_lat)
                            idx_a = y; idx_b = y; 
                            break;
                        elseif(lat_wrf(y,x_mid) > slice_lat)
                            idx_a = max(1,y-1); idx_b = y;
                            break;
                        end
                    end
                    if(exp_idx == 1)
                        pressure_slice = squeeze(pressure_a(idx_a,:,:));
                        %pressure_max = max(pressure_a,[],"all");
                    else
                        pressure_slice = squeeze(pressure_b(idx_a,:,:));
                        %pressure_max = max(pressure_b,[],"all");
                    end
                    for x = 1:xlen_wrf
                        slice_lon_list(x) = (lon_wrf(idx_a,x)+lon_wrf(idx_b,x))/2;
                    end
                end



                slice_lat_grid = meshgrid(slice_lat_list);
                slice_lat_grid = slice_lat_grid(1:zlen_wrf,:)';
                slice_lon_grid = meshgrid(slice_lon_list);
                slice_lon_grid = slice_lon_grid(1:zlen_wrf,:)';

                %% Compute radar beam intersections

                % Existing model height field name is "model_height" (above ground)
                %load(sprintf('%s/%s',intermediate_path,radar_heights_filename),'radar_height_grid');
                %load(sprintf('%s/%s',intermediate_path,wrf_heights_filename),'wrf_height_map');
                %beam_angle_l = beam_angle - (beam_width/2);
                %beam_angle_u = beam_angle + (beam_width/2);

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
                radar_height_grid = radar_height_grid_full(radars_used,:,:);

                % Find each intersection
                if(y_axis)
                    intersection_pressure_mult = zeros(num_radars,ylen_wrf);
                    refl_value_seen = zeros(size(intersection_pressure_mult));
                    model_height_interp = zeros(size(slice_values));
                    data_interp = zeros(size(slice_values));
                    for rad_idx = 1:num_radars
                        for y = 1:ylen_wrf
                            rad_height = linear_interpolate(lon_wrf(y,idx_a),lon_wrf(y,idx_b),slice_lon,radar_height_grid(rad_idx,y,idx_a),radar_height_grid(rad_idx,y,idx_b));
                            for z = 1:zlen_wrf
                                model_height_interp(y,z) = linear_interpolate(lon_wrf(y,idx_a),lon_wrf(y,idx_b),slice_lon,model_height(y,idx_a,z),model_height(y,idx_b,z));
                                data_interp(y,z) = linear_interpolate(lon_wrf(y,idx_a),lon_wrf(y,idx_b),slice_lon,data(y,idx_a,z),data(y,idx_b,z));
                                if(model_height_interp(y,z) > rad_height)
                                    intersection_pressure_mult(rad_idx,y) = linear_interpolate(model_height_interp(y,max(z-1,1)),model_height_interp(y,z),rad_height,pressure_slice(y,max(z-1,1)),pressure_slice(y,z));
                                    break;
                                end
                            end
                        end  
                    end              
                else
                    intersection_pressure_mult = zeros(num_radars,xlen_wrf);
                    refl_value_seen = zeros(size(intersection_pressure_mult));
                    model_height_interp = zeros(size(slice_values));
                    data_interp = zeros(size(slice_values));
                    for rad_idx = 1:num_radars
                        for x = 1:xlen_wrf
                            rad_height = linear_interpolate(lat_wrf(idx_a,x),lat_wrf(idx_b,x),slice_lat,radar_height_grid(rad_idx,idx_a,x),radar_height_grid(rad_idx,idx_b,x));
                            for z = 1:zlen_wrf
                                model_height_interp(x,z) = linear_interpolate(lat_wrf(idx_a,x),lat_wrf(idx_b,x),slice_lat,model_height(idx_a,x,z),model_height(idx_b,x,z));
                                data_interp(x,z) = linear_interpolate(lat_wrf(idx_a,x),lat_wrf(idx_b,x),slice_lat,data(idx_a,x,z),data(idx_b,x,z));
                                if(model_height_interp(x,z) > rad_height)
                                    intersection_pressure_mult(rad_idx,x) = linear_interpolate(model_height_interp(x,max(z-1,1)),model_height_interp(x,z),rad_height,pressure_slice(x,max(z-1,1)),pressure_slice(x,z));
                                    break;
                                end
                            end
                        end  
                    end    
                end

                %% Compute real base reflectivity heights
                
                % Determine which radar refl value is greatest, take that
                % beam's intersection height and compute pressure
                if(y_axis) % Slice is along y-axis
                    height_values = zeros(ylen_wrf,1);
                    intersection_pressure = zeros(ylen_wrf,1);
                    intersection_pressure(:) = NaN;
                    for y = 1:ylen_wrf
                        max_dbz_idx = 0;
                        max_dbz = -999;
                        for rad_idx = 1:num_radars
                            rad_dbz = -999;
                            rad_height = linear_interpolate(lon_wrf(y,idx_a),lon_wrf(y,idx_b),slice_lon,radar_height_grid(rad_idx,y,idx_a),radar_height_grid(rad_idx,y,idx_b));
                            for z = 1:zlen_wrf
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
                        if(max_dbz_idx > 0)
                            height_values(y) = linear_interpolate(lon_wrf(y,idx_a),lon_wrf(y,idx_b),slice_lon,radar_height_grid(max_dbz_idx,y,idx_a),radar_height_grid(max_dbz_idx,y,idx_b));
                            pressure_found = 0;
                            for z = 1:zlen_wrf
                                slice_values(y,z) = linear_interpolate(lon_wrf(y,idx_a),lon_wrf(y,idx_b),slice_lon,data(y,idx_a,z),data(y,idx_b,z));
                                if(model_height_interp(y,z) >= height_values(y) && ~pressure_found)
                                    intersection_pressure(y) = linear_interpolate(model_height_interp(y,max(z-1,1)),model_height_interp(y,z),height_values(y),pressure_slice(y,max(z-1,1)),pressure_slice(y,z));
                                    pressure_found = 1;
                                end
                            end
                        end
                    end
                else % Slice is along x-axis
                    height_values = zeros(xlen_wrf,1);
                    intersection_pressure = zeros(xlen_wrf,1);
                    intersection_pressure(:) = NaN;
                    for x = 1:xlen_wrf
                        max_dbz_idx = 0;
                        max_dbz = -999;
                        for rad_idx = 1:num_radars
                            rad_dbz = -999;
                            rad_height = linear_interpolate(lat_wrf(idx_a,x),lat_wrf(idx_b,x),slice_lat,radar_height_grid(rad_idx,idx_a,x),radar_height_grid(rad_idx,idx_b,x));
                            for z = 1:zlen_wrf
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
                        if(max_dbz_idx > 0)
                            height_values(x) = linear_interpolate(lat_wrf(idx_a,x),lat_wrf(idx_b,x),slice_lat,radar_height_grid(max_dbz_idx,idx_a,x),radar_height_grid(max_dbz_idx,idx_b,x));
                            pressure_found = 0;
                            for z = 1:zlen_wrf
                                slice_values(x,z) = linear_interpolate(lat_wrf(idx_a,x),lat_wrf(idx_b,x),slice_lat,data(idx_a,x,z),data(idx_b,x,z));
                                if(model_height_interp(x,z) >= height_values(x) && ~pressure_found)
                                    intersection_pressure(x) = linear_interpolate(model_height_interp(x,max(z-1,1)),model_height_interp(x,z),height_values(x),pressure_slice(x,max(z-1,1)),pressure_slice(x,z));
                                    pressure_found = 1;
                                end
                            end
                        end
                    end
                end

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
                    h = pcolor(slice_lat_grid,pressure_slice,data_to_plot); % Plot the data
                else
                    h = pcolor(slice_lon_grid,pressure_slice,data_to_plot); % Plot the data
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
                            h = plot(slice_lat_list,intersection_pressure_mult(rad_idx,:),'Color',multibeam_color,'LineWidth',2,'LineStyle','--');
                        end
                    else
                        for rad_idx = 1:num_radars
                            h = plot(slice_lon_list,intersection_pressure_mult(rad_idx,:),'Color',multibeam_color,'LineWidth',2,'LineStyle','--');
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
                        %fprintf('Pmin: %f | Pmax: %f\n',pressure_min,pressure_max); % WIP DEBUG
                        ylim([pressure_min pressure_max]); % WIP LOCATION: PRESSURE_MAX IS 8.1 (TOO SMALL) ON MEM 27 TESTING
                    else
                        xlim([plot_w_slice_edge plot_e_slice_edge]);
                        ylim([pressure_min pressure_max]);
                    end
                end

                % Apply labels

                if(y_axis)
                    title(sprintf('[%s|d0%d|%s] Refl VSlice #%d (Lon %02.02f) | %s',exp_name_plot,domain,member_string,slice_idx,slice_lon,timestamp_title),'FontSize',title_font_size); 
                    xlabel('Latitude','FontSize',label_font_size);
                    ylabel('Pressure (mb)','FontSize',label_font_size);
                else
                    title(sprintf('[%s|d0%d|%s] Refl VSlice #%d (Lat %02.02f) | %s',exp_name_plot,domain,member_string,slice_idx,slice_lat,timestamp_title),'FontSize',title_font_size); 
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
                hold off;
                close('all');

                clearvars f h data_to_plot;

            end

            %% Plot soundings/vertical profiles

            if(make_wrf_soundings)
                for st_idx = 1:num_st

                    lon_st = st_lon_list(st_idx);
                    lat_st = st_lat_list(st_idx);
                    edge_buffer = 42;

                    rh_min = 0;
                    rh_max = 114.2;
                    t_max = 10;
                    t_min = -70;

                    load(sprintf('%s/lv_table.mat',intermediate_path),'T','LV5'); % C, J/kg
                    Lv_chart = LV5;
                    T_chart = T;
    
                    if(exp_idx == 1)
                        data_mr = data_mr_a;
                        data_t = data_t_a;
                        data_p = pressure_a;
                        plot_filename = sprintf("%s/skewT_%s_%s_%s_%s.png",output_path_st,exp_name_a,member_string,st_name_list(st_idx),timestamp_file); % 
                    else
                        data_mr = data_mr_b;
                        data_t = data_t_b;
                        data_p = pressure_b;
                        plot_filename = sprintf("%s/skewT_%s_%s_%s_%s.png",output_path_st,exp_name_b,member_string,st_name_list(st_idx),timestamp_file); % 
                    end

                    % Interpolate from grid to column at current location
                    mR = zeros([zlen_wrf 1]);
                    T = zeros([zlen_wrf 1]);
                    P = zeros([zlen_wrf 1]);
                    [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE] = find_nn_idx_irregular(lon_st,lat_st,lon_wrf,lat_wrf,edge_buffer);
                    nn_points = [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE];    
                    if(any(isnan(nn_points)) || any(nn_points == 0) || any(nn_points(1:4) > xlen_wrf) || any(nn_points(5:8) > ylen_wrf)) % If indices could not be found / reference points are NaN/ goes off edge
                        fprintf("Could not interpolate to Skew-T location #%d\n",st_idx);
                        continue;
                    end
                    for z = 1:zlen_wrf
                        P(z) = bilinear_interpolate_pt(lat_st,lon_st,lat_wrf,lon_wrf,data_p(:,:,z),nn_points);
                        mR(z) = bilinear_interpolate_pt(lat_st,lon_st,lat_wrf,lon_wrf,data_mr(:,:,z),nn_points);
                        T(z) = bilinear_interpolate_pt(lat_st,lon_st,lat_wrf,lon_wrf,data_t(:,:,z),nn_points);
                    end

                    fig_width_st = 700;
                    fig_height_st = 700;

                    % Constants
                    T0 = 273.15; % K
                    c = 243.04;
                    b = 17.625;
                    Rv = 461.52; % J/kgK; Lv = J/kg
                    eS0 = 611; % Pa

                    Lv = compute_Lv(T,Lv_chart,T_chart);

                    %-------------
                    if(test_lv_eqn)
                        test_P = [1004 1000 992 977 925 913 911 891 877 868]; % hPa
                        test_T = [4.8 5.6 6.4 6.0 2.6 1.8 1.8 1.4 3.7 5.4]; % C
                        test_TD = [0.2 -1.4 -0.6 -2.0 -3.4 -4.2 -4.3 -5.6 -15.9 -23.6]; % C
                        test_MR = [3.88 3.47 3.71 3.39 3.23 3.08 3.07 2.84 1.27 0.66];
    
                        test_Tk = test_T + T0;
                        test_Tc = test_T;test_P = [1004 1000 992 977 925 913 911 891 877 868]; % hPa
                        test_T = [4.8 5.6 6.4 6.0 2.6 1.8 1.8 1.4 3.7 5.4]; % C
                        test_TD = [0.2 -1.4 -0.6 -2.0 -3.4 -4.2 -4.3 -5.6 -15.9 -23.6]; % C
                        test_MR = [3.88 3.47 3.71 3.39 3.23 3.08 3.07 2.84 1.27 0.66];
    
                        test_Tk = test_T + T0;
                        test_Tc = test_T;
                        test_P = test_P.*100; % Pa
                        test_Lv = compute_Lv(test_Tk,Lv_chart,T_chart);
    
                        eS_a = 611*exp((17.67*(test_Tk-T0))./(test_Tk-29.65)); % Pa
                        eS_b = eS0*exp((test_Lv./Rv).*((1/T0)-(1./test_Tk)));
                        eS = eS_b;
                        mS = 621.97*(eS./(test_P-eS)); % g/kg???
                        RH = test_MR./mS; % decimal, not x/100
                        lam = log(RH) + b*test_Tc./(c + test_Tc);
                        test_TD_eqn = c*lam./(b-lam);
                        test_P = test_P.*100; % Pa
                        test_Lv = compute_Lv(test_Tk,Lv_chart,T_chart);
    
                        eS_a = 611*exp((17.67*(test_Tk-T0))./(test_Tk-29.65)); % Pa
                        %eS_b = eS0*exp((test_Lv./Rv).*((1/T0)-(1./test_Tk)));
                        eS = eS_a;
                        mS = 621.97*(eS./(test_P-eS)); % g/kg???
                        RH = test_MR./mS; % decimal, not x/100
                        lam = log(RH) + b*test_Tc./(c + test_Tc);
                        test_TD_eqn = c*lam./(b-lam);
                    end
                    %------------
                  
                    % Unit conversions
                    Tk = T; % Temperature in Kelvin
                    Tc = T-T0; % Temperature in Celsius
                    P_Pa = P.*100; % Air Pressure hPa -> Pa
                    mR_gkg = mR.*1000; % kg/kg -> g/kg

                    use_long_eqn = 1;

                    eS = 611*exp((17.67*(Tk-T0))./(Tk-29.65)); % Pa
                    %eS = eS0*exp((Lv./Rv).*((1/T0)-(1./T)));
                    mS = 621.97*(eS./(P_Pa-eS)); % g/kg???
                    %mS = 0.622*eS./P_Pa;
                    RH = mR_gkg./mS;
                    %lam = log(RH/100) + b*Tc./(c + Tc);

                    if(use_long_eqn)
                        % a = 6.1121 mbar, b = 18.678, c = 257.14 °C, d = 234.5 °C.
                        a = 6.1121;
                        b = 18.678;
                        c = 257.14;
                        d = 234.5;
                        lam = log(RH.*exp( (b-(Tc/d)).*(Tc./(c+Tc)) ));
                    else
                        lam = log(RH) + b*Tc./(c + Tc);
                    end

                    TD = c*lam./(b-lam);
                    

                    % DEBUG
                    %for idx = 1:length(TD)
                    %    fprintf('P = %04.03f, mS = %01.03f\n',P(idx),mS(idx))
                    %end
        
                    f = skewt(P,Tc,TD,'TempOpts',{'color','red','linewidth',2},'DwptOpts',{'color','green','linewidth',2}); % DOES NOT WORK
                    f.Position = [fig_x fig_y fig_width_st fig_height_st]; % Shrink figure
                    title(sprintf('[%s|d%d|%s] %s | %dz',exp_name_plot,domain,member_string,st_name_list(st_idx),hour),'FontSize',title_font_size); 
                    xlabel('Temperature','FontSize',label_font_size);
                    ylabel('Pressure (mb)','FontSize',label_font_size);
                    set(gca,'Fontsize',axes_font_size);

                    saveas(f,plot_filename); % Save as .png
                    close('all');

                    % Custom sounding plots
                    plot_filename = sprintf("%s/custom/TxRH_%s_%s_%s_%s.png",output_path_st,exp_name_plot,member_string,st_name_list(st_idx),timestamp_file); % 
                    f = figure('Position',[fig_x fig_y fig_width_txrh fig_height_txrh]); % Create initial blank figure
                    yyaxis left;
                    z = plot([0 1100],[0 0],'Color','#3b3b3b','LineStyle','--','LineWidth',1);
                    hold on;
                    plot(P,Tc,'r','LineStyle','-','LineWidth',3);
                    ylabel(sprintf('Temperature (%sC)',char(176)));
                    xlabel('Pressure(mb)')
                    ylim([t_min t_max]);
                    set(gca,'YColor','r');

                    yyaxis right;
                    plot(P,RH*100,'Color','#0270d1','LineWidth',3);
                    ylabel('Relative Humidity (%)');
                    ylim([rh_min rh_max]);
                    set(gca,'YColor','#0270d1');
                    xlim([50 1000]);
                    set(gca, 'XDir','reverse')
                    title(sprintf('[%s|d%d|%s] %s | %dz',exp_name_plot,domain,member_string,st_name_list(st_idx),hour),'FontSize',title_font_size); 
                    set(gca,'Fontsize',axes_font_size);
                    
                    saveas(gcf,plot_filename); % Save as .png
                    close('all');

                    clearvars Tc RH P TD;

                    %% Test IMPACTs sounding data

                    if(member_idx == 1) % Only once per timestamp
                        exp_name_o = 'OBS';
                        
                        t_name_obs = 'temp';
                        rh_name_obs = 'relh';
                        p_name_obs = 'pres';
    
                        switch st_idx
                            case 1
                                filename = sprintf('%s/IMPACTS_sounding_20200207_150000_ALB.nc',input_path_st);
                            case 2
                                filename = sprintf('%s/IMPACTS_sounding_20200207_150000_BUF.nc',input_path_st);
                            case 3
                                filename = sprintf('%s/IMPACTS_sounding_20200207_150000_GYX.nc',input_path_st);
                            case 4
                                filename = sprintf('%s/IMPACTS_UIUC_Mobile_research_sounding_20200207_1500.nc',input_path_st);
                                t_name_obs = 'TC';
                                rh_name_obs = 'RH';
                                p_name_obs = 'PRESS';
                        end
        
    
            
                        TC = ncread(filename,t_name_obs);
                        RH = ncread(filename,rh_name_obs);
                        P = ncread(filename,p_name_obs);
            
                        plot_filename = sprintf("%s/custom/TxRH_%s_%s_%s.png",output_path_st,exp_name_o,st_name_list(st_idx),timestamp_file); % 
                        f = figure('Position',[fig_x fig_y fig_width_txrh fig_height_txrh]); % Create initial blank figure
                        yyaxis left;
                        z = plot([0 1100],[0 0],'Color','#3b3b3b','LineStyle','--','LineWidth',1);
                        hold on;
                        plot(P,TC,'r','LineStyle','-','LineWidth',3);
                        ylabel(sprintf('Temperature (%sC)',char(176)));
                        xlabel('Pressure(mb)')
                        ylim([t_min t_max]);
                        set(gca,'YColor','r');
    
                        yyaxis right;
                        plot(P,RH,'Color','#0270d1','LineWidth',3);
                        ylabel('Relative Humidity (%)');
                        ylim([rh_min rh_max]);
                        set(gca,'YColor','#0270d1');
                        xlim([50 1000]);
                        set(gca, 'XDir','reverse')
                        title(sprintf('[%s] %s | %dz',exp_name_o,st_name_list(st_idx),hour),'FontSize',title_font_size); saveas(f,plot_filename); % Save as .png
                        set(gca,'Fontsize',axes_font_size);
    
                        
                        saveas(gcf,plot_filename);
                        close('all');
                        clearvars TC RH P;
                    end

                end
            end
            

        end
        %% Clear big variables between members
        clearvars rmse_err bias pressure* data* -except pressure_min pressure_max data_storage* data_type data_name* data_gis_trim data_src* data_wrf_mean_ref;
        
        
    end

    clearvars ens_mem* obs_* member_sum* data_wrf_mean_ref model_height* terrain_height;
    
    %% Time increment
    hour = hour + hour_step;

    if(hour > 23)
        hour = 0;
        day = day + 1;
    end

    if(day > max_day(month))
        day = 1;
        month = month + 1;
    end

    if (month > 12)
        month = 1;
        year = year + 1;
    end

end

% Get total runtime and print to stdout
runtime = toc;
hours = floor(runtime/3600);
mins = floor((runtime/60) - (hours*60));
secs = toc - (hours*3600) - (mins*60);
fprintf('Done. Total script runtime = %02.f:%02.f:%02.f\n',hours,mins,secs)

% END


%% Local functions

% Process input data to align with customized 2D grid
function out_data = align_data(in_data,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx)

    temp_data = in_data;

    % If data needs to be composited down to 2D, do so
    if(z_composite_data)
        temp_data = max(temp_data,[],3);
    end

    % If data needs to be transposed and/or flipped to align with lat/lon grid, do so
    if(transpose_data) 
        temp_data = permute(temp_data,[2 1 3]);
    end

    if(flip_data)
        temp_data = flip(temp_data);
    end

    % Trim data to desired spatial limits
    out_data = temp_data(n_idx:s_idx,w_idx:e_idx,:);
end
                