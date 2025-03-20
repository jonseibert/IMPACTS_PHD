% refl_struct_analysis (v2.0.8)
%
% Purpose: Generate 2D map of the vertical heights / pressure levels of the
% maximum radar reflectivity per column & vertical slices
% 
% Author(s): Jon Seibert
% Last updated: 15 August 2024
% 
% Inputs: [WRF Outputs].nc, [Observations].png
% Outputs: 
% Dependencies: reflmap.m, borders.m
%
% NOTES:
%  - TIME_COUNT is not robust when crossing month boundaries.
%  - Not robust for use in southern hemisphere (lat < 0)
%
% TODO:
%  - 
script_name = 'refl_struct_analysis.m';
   
%% 0A. Script Controls

run_name = 'vslice';

% Local or server execution mode
local = 1; % true = running on local device, false = running on server

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
start_hour = 15;
end_hour = 15;

domain = 2;
subdomain = 0; %1 for E Coast Only subdomain, 2 for N only, 0 for all d02 included

smooth_data = 0; % If true AND the variable isn't already being smoothed, smooth it before plotting (gaussian)

% Vertical levels to use if is 3D
%pressure_levels = [460,500,700,850,900,925];
pressure_target = 500; %mb

multi_height = 0;
pressure_target_list = [925,850,700]; % Overrides above value

% REFL_10CM/REFL (dBZ), T/T (deg C), P/MSLP (mb), GPH/GPH(m),
% VORT/VORT(10^-5 s^-1), U-V-W/WIND(kt)
is_3d = 1;
use_base_ref = 0;
use_contours = 1;

plot_raw = 1;
plot_dif = 1;

label_contours = 1;
gauss_sd = 2.75;
gauss_width =  2*ceil(2*gauss_sd)+1;

% Note whether transposition or compositing is necessary
transpose_data = true;
flip_data = true;
z_composite_data = false; % Unused if use_base_ref is true

limit_raw_colorbar = false;
clim_lower = -30; % Colorbar limits
clim_upper = 75;
clim_lower_dif = -10; %  Overridden by some choices of data variables
clim_upper_dif = 10;

% Plot selections
remake_plots = 1; % If true, overrides existing outputs. If false, only plots 'new' files.

show_plots = 1;
show_progress = 1; % If true, prints out progress bar/debug messages
announce_plots = 0;
use_auto_run_name = 0;

%num_members = 40;
num_members = 2; % TESTING

% Vertical slice boundaries
SW_lon = -75.7; 
SW_lat = 42; 
NE_lon = -73.4;
NW_lat = 43.8;
min_z = 1;
max_z = 50;

%% 0B. General Settings

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020';
    input_path_gis = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    intermediate_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/testing';
    path_to_code = "C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code";
    path_to_extra_code = './downloaded_code';
    input_path_ecmwf = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Data/ECMWF';
    % Experiment paths
    exp_path_1 = ''; % TESTING
    exp_path_2 = 'CONV/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
    exp_path_6 = 'GEFS/';
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/';
    input_path_gis = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/GIS';
    intermediate_path = '/storage/home/jjs5895/projects/IMPACTS/intermediate';
    output_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/ptp';
    path_to_code = "/storage/home/jjs5895/projects/IMPACTS/code";
    path_to_extra_code = '/storage/home/jjs5895/projects/IMPACTS/code/downloaded_code'; % Specify path to borders.m
    input_path_ecmwf = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/ECMWF';
    % Experiment paths
    exp_path_1 = 'AIRCFT/fc';
    exp_path_2 = 'CONV/fc/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
    exp_path_6 = 'GEFS/';
end

addpath(path_to_extra_code);

% Experiment names
exp_name_1 = 'AIR';
exp_name_2 = 'CONV';
exp_name_3 = 'AIR-F';
exp_name_4 = 'CONV-F';
exp_name_5 = 'NODA';
%exp_name_6 = 'GEFS';

alt_input_exps = [3 4 5]; % Input is set up differently for these experiments

hour_step = 1; % Hour increment size

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 925;
fig_height = 900;
fig_width_small = 350;
fig_height_small = 350;

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
trim_e = false; % Whether to trim the eastern edge

% Figure font sizes
title_font_size = 18;
label_font_size = 18;
axes_font_size = 16;
contour_font_size = 14;

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
data_name_q = "QVAPOR"; % Water vapor mixing ratio, kg/kg
data_name_u = "U"; % Wind [X Stagger]
data_name_v = "V"; % Wind [Y Stagger]
data_name_w = "W"; % Wind [Z Stagger]

Re = 6.3781e6; % Radius pf Earth in m
g = 9.80665; % m/s^2 standard gravity at sea level
R = 287.0600676; % J/kg/K, specific gas constant

%% 1C. GEFS Model information

%%% WIP %%%

% ECMWF reanalysis filepath

ecmwf_filename = 'ecmwf_mslp_d01_20200207.nc';
        
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

% If output folders do not exist, create them
if(~isfolder(output_path_large) || ~isfolder(output_path_small))
    mkdir(output_path_large);
    mkdir(output_path_small);
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
        data_name = data_name_q;
        if(~multi_height && pressure_target == 500)
            clim_lower = 0;
            clim_upper = 5e-3;
        elseif(~multi_height && pressure_target > 500)
            clim_lower = 0;
            clim_upper = 1.5e-2;
        elseif(multi_height && any(pressure_target_list > 500))
            clim_lower = 0;
            clim_upper = 1.5e-2;
        else
            clim_lower = 0;
            clim_upper = 5e-3;
        end
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
    
    if(isfile(sprintf('%s/latlon_wrf_trim_d01.mat',intermediate_path)))
        load(sprintf('%s/latlon_wrf_trim_d01.mat',intermediate_path));
    else
    
        lon_wrf_raw = ncread(sprintf('%s/%s/%s',input_path_a,timestamp_file,filename_wrf),lon_name_wrf);
        lat_wrf_raw = ncread(sprintf('%s/%s/%s',input_path_a,timestamp_file,filename_wrf),lat_name_wrf);

        if(transpose_data)
            lon_wrf_raw = lon_wrf_raw';
            lat_wrf_raw = lat_wrf_raw';
        end
        if(flip_data)
            lon_wrf_raw = flip(lon_wrf_raw);
            lat_wrf_raw = flip(lat_wrf_raw);
        end

        dims = size(lon_wrf_raw); % Determine grid dimensions

        % Define lat/lon boundaries to use when trimming (refer to Settings)
        w_bound = w_lim_d01;
        e_bound = e_lim_d01;
        s_bound = s_lim_d01;
        n_bound = n_lim_d01;

        % IN THIS DOMAIN, WITH MY ORIENTATION:
        % WRF Lon primarily increases with column, secondarily decreases with row
        % WRF Lat primarily decreases with row, secondarily decreases with column
        % Matlab is [row,col] ~ [Y,X]/[lat,lon]

        % Find W border index
        % Theoretically, Max row = furthest S- should contain the first instance of within-domain lon
        % Check both ends to be safe
        row = dims(1); 
        for col = 1:dims(2)
            if(lon_wrf_raw(row,col) > w_bound)
                w_idx_a = col;
                break;
            end
        end 

        row = 1;
        for col = 1:dims(2)
            if(lon_wrf_raw(row,col) > w_bound)
                w_idx_b = col;
                break;
            end
        end 

        w_idx = min(w_idx_a,w_idx_b); % Whichever index is further W

        % Find E border index
        row = dims(1); % 
        for col = dims(2):-1:1
            if(lon_wrf_raw(row,col) < e_bound)
                e_idx_a = col;
                break;
            end
        end

        row = 1;
        for col = dims(2):-1:1
            if(lon_wrf_raw(row,col) < e_bound)
                e_idx_b = col;
                break;
            end
        end

        e_idx = min(e_idx_a,e_idx_b);

        % Find N border index
        col = dims(2); 
        for row = 1:dims(1)
            if(lat_wrf_raw(row,col) < n_bound)
                n_idx_a = row;
                break;
            end
        end

        col = 1;
        for row = 1:dims(1)
            if(lat_wrf_raw(row,col) < n_bound)
                n_idx_b = row;
                break;
            end
        end 

        n_idx = min(n_idx_a,n_idx_b);

        % Find S border index
        col = dims(2); 
        for row = dims(1):-1:1
            if(lat_wrf_raw(row,col) > s_bound)
                s_idx_a = row;
                break;
            end
        end 

        col = 1;
        for row = dims(1):-1:1
            if(lat_wrf_raw(row,col) > s_bound)
                s_idx_b = row;
                break;
            end
        end 

        s_idx = max(s_idx_a,s_idx_b);

        if(~trim_e) % If Eastern cutoff is not desired
            e_idx = dims(2);
        end

        lon_wrf_trim = lon_wrf_raw(n_idx:s_idx,w_idx:e_idx);
        lat_wrf_trim = lat_wrf_raw(n_idx:s_idx,w_idx:e_idx);

        save(sprintf('%s/latlon_wrf_trim_d01.mat',intermediate_path),'n_idx','s_idx','e_idx','w_idx','lon_wrf_trim','lat_wrf_trim');
    end
    
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
data_storage_a = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
data_storage_b = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);


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

    % Read in each member and compute the mean
    for member_idx = 1:(num_members+1)
        
        % If all members have been computed, take mean
        if(member_idx == (num_members+1)) 
            member_string = 'mean';
            data_a = mean(ens_mem_a,4,"omitnan");
            data_b = mean(ens_mem_b,4,"omitnan");
            
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

            data_pb_b = ncread(full_filepath_b,data_name_pres);
            data_pp_b = ncread(full_filepath_b,data_name_pres_prime); 
            data_ptp_b = ncread(full_filepath_b,data_name_temp); 
            data_p_b = data_pb_b + data_pp_b; 
            data_pt_b = data_ptp_b + 300; 
            data_t_b = data_pt_b.*((data_p_b./100000).^(0.287));
            data_gpb_raw_b = ncread(full_filepath_b,data_name_gp); 
            data_gpp_raw_b = ncread(full_filepath_b,data_name_gp_prime); 
            data_p_b = data_p_b./100; 

            terrain_height = ncread(full_filepath_a,data_name_terrain_height); % Height of physical terrain (m)

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

            model_height = ((data_gpb_a+data_gpp_a)./g) - terrain_height; % model height above the *ground* (m)

            % Align all data fields with customized grid
            pressure_a = align_data(data_p_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            pressure_b = align_data(data_p_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            
            % If using base reflectivity
            if(use_base_ref)
                load(sprintf(input_format_base_refl,input_path_base_refl_a,data_src_a,'base',domain,member_string,data_type,timestamp_file),'wrf_base_refl');
                data_a = wrf_base_refl;
                clearvars wrf_base_refl;
                load(sprintf(input_format_base_refl,input_path_base_refl_b,data_src_b,'base',domain,member_string,data_type,timestamp_file),'wrf_base_refl');
                data_b = wrf_base_refl;
                clearvars wrf_base_refl;
            else     
                data_a = ncread(full_filepath_a,data_name); % Retrieve data
                data_b = ncread(full_filepath_b,data_name); % Retrieve data
            end
            
            dims = [ylen_wrf xlen_wrf zlen_wrf];
            
            % All data must be 2D to be processed and plotted. If not
            % composite or interpolated already:
            if(~use_base_ref)
                data_a = align_data(data_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_b = align_data(data_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            end

            % Store values for next loop
            ens_mem_a(:,:,:,member_idx) = data_a;
            ens_mem_b(:,:,:,member_idx) = data_b;

        end
        
        
        %% Compute height map
        height_map_a = zeros(ylen_wrf,xlen_wrf);
        height_map_b = zeros(ylen_wrf,xlen_wrf);
        
        [max_values,max_idx] = max(data_a,[],3);
        for y = 1:ylen_wrf
            for x = 1:xlen_wrf
                height_map_a(y,x) = pressure_a(y,x,max_idx(y,x));
            end
        end
    
        [max_values,max_idx] = max(data_b,[],3);
        for y = 1:ylen_wrf
            for x = 1:xlen_wrf
                height_map_b(y,x) = pressure_b(y,x,max_idx(y,x));
            end
        end
        
       %% Plot vertical slices: find line
       

       %% Apply area mask
        
        if(show_progress)
            toc
            fprintf('Masking...\n'); % debug statement
        end

        load('mask_wrf_on_wrf','wrf_wrf_mask');
        if(subdomain == 0)
            wrf_sbd_mask = ones(size(lon_wrf));
            gis_sbd_mask = ones(size(lon_gis_trim));
        elseif(subdomain == 1)
            wrf_sbd_mask = inpolygon(lon_wrf,lat_wrf,bound_lon_d1,bound_lat_d1);
            gis_sbd_mask = inpolygon(lon_gis_trim,lat_gis_trim,bound_lon_d1,bound_lat_d1);
        elseif(subdomain == 2)
            wrf_sbd_mask = inpolygon(lon_wrf,lat_wrf,bound_lon_d2,bound_lat_d2);
            gis_sbd_mask = inpolygon(lon_gis_trim,lat_gis_trim,bound_lon_d2,bound_lat_d2);
        else
            fprintf('ERROR: Invalid subdomain.\n');
            exit(1);
        end
            
        data_a_masked = data_a;
        data_a_masked(data_a_masked == 0) = 42069;
        data_a_masked = data_a_masked.*wrf_sbd_mask;
        data_a_masked(data_a_masked == 0) = nodata_val_gis;
        data_a_masked(data_a_masked == 42069) = 0;
        data_a_compare = data_a_masked;
        
        data_b_masked = data_b;
        data_b_masked(data_b_masked == 0) = 42069;
        data_b_masked = data_b_masked.*wrf_sbd_mask;
        data_b_masked(data_b_masked == 0) = nodata_val_gis;
        data_b_masked(data_b_masked == 42069) = 0;
        data_b_compare = data_b_masked;
       
        clearvars wrf_gis_mask wrf_wrf_mask;

        %% Make dif, compute RMSE, bias, ETS
        
        if(show_progress)
            toc
            fprintf('Computing error values...\n'); % debug statement
        end

        data_dif = data_a_compare - data_b_compare;
        data_dif_clean = data_dif(~isnan(data_dif));
        num_clean_points = numel(data_dif_clean);
        num_points = numel(data_dif);
        bias = mean(data_dif_clean,'all');
        rmse_err = sqrt(sum(data_dif_clean.^2)/num_clean_points);

        %% Store for specific plots
        data_storage_a(:,:,:,member_idx) = data_a_compare;
        data_storage_b(:,:,:,member_idx) = data_b_compare;

        %% WRF-GIS.2: Plot

        if(show_progress)
            toc
            fprintf('Plotting...\n'); % debug statement
        end
        
        gauss_filt = fspecial('gaussian',gauss_width,gauss_sd);
        
        for plot_idx = 1:3
            if(plot_idx == 1)
                if(~plot_raw) continue; end
                plot_type = 'raw';  % Specify fit, dif, trim, raw
                data_to_plot = height_map_a;
                exp_name_plot = exp_name_a_plot;
                cmap = 'jet';
                plot_filename = sprintf(output_format_var,output_path_large,lower(exp_name_a_plot),lower(exp_name_a_plot),domain,subdomain,member_string,bao_short,plot_type,data_type,pressure_target,timestamp_file);
                plot_filename_small = sprintf(output_format_var,output_path_small,lower(exp_name_a_plot),lower(exp_name_a_plot),domain,subdomain,member_string,bao_short,plot_type,data_type,pressure_target,timestamp_file);
            elseif(plot_idx == 2)
                if(~plot_raw) continue; end
                plot_type = 'raw';  % Specify fit, dif, trim, raw
                data_to_plot = height_map_b;
                exp_name_plot = exp_name_b_plot;
                cmap = 'jet';
                plot_filename = sprintf(output_format_var,output_path_large,lower(exp_name_b_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,data_type,pressure_target,timestamp_file);
                plot_filename_small = sprintf(output_format_var,output_path_small,lower(exp_name_b_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,data_type,pressure_target,timestamp_file);
            elseif(plot_idx == 3)
                if(~plot_dif) continue; end
                if(compare_to_background)
                    plot_type = 'inc';
                else
                    plot_type = 'dif';
                end
                data_to_plot = height_map_a - height_map_b;
                cmap = 'redblue';
                plot_filename = sprintf(output_format_var,output_path_large,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,data_type,pressure_target,timestamp_file);
                plot_filename_small = sprintf(output_format_var,output_path_small,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,data_type,pressure_target,timestamp_file);
            else
                fprintf('ERROR: You should not ever see this message. You have a bug in the plot code.\n');
            end

            if(isfile(plot_filename) && ~remake_plots) % If shouldn't override existing plots, don't
                continue;
            end

            % Plot LARGE
            f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
            
            if(show_progress && announce_plots)
                toc
                fprintf('%s...\n',plot_type); % debug statement
            end

            %if(smooth_data)
                %data_to_plot = nanconv(data_to_plot,gauss_filt,'nanout','edge');
                %clim_mod = 1/gauss_sd;
            %else
                clim_mod = 1;
            %end

            if(use_contours)
                data_smooth = nanconv(data_to_plot,gauss_filt,'nanout','edge'); % Keep contours nice and clean
                [M,h] = contour(lon_wrf,lat_wrf,data_smooth,'EdgeColor',"#0031d9","LineWidth",3,"LevelStep",40); % Plot the data
                if(label_contours)
                    clabel(M,h,'FontSize',contour_font_size,'LabelSpacing',250);
                end
            else
                h = pcolor(lon_wrf,lat_wrf,data_to_plot); % Plot the data
                set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                shading interp; % Smooth out plot from grid-boxes
                c = colorbar('FontSize',axes_font_size); % Make colorbar
                colormap(cmap); % Set colors
                ylabel(c,'Pressure at height of maximum dBZ (mb)');
            end
            

            if(plot_type == "dif" || plot_type == "inc")
                caxis([clim_lower_dif*clim_mod clim_upper_dif*clim_mod]);% Plotting dif
            else
                if(limit_raw_colorbar)
                    caxis([clim_lower*clim_mod clim_upper*clim_mod]);% Plotting dif
                end
            end

            if(show_progress && announce_plots)
                toc
                fprintf('Borders...\n'); % debug statement
            end

            % Plot state borders
            borders('continental us','black','linewidth',1); 
            hold off;

            if(show_progress && announce_plots)
                toc
                fprintf('Limits...\n'); % debug statement
            end

            if(limit_borders)
                xlim([w_lim e_lim]);
                ylim([s_lim n_lim]);
            end

            if(show_progress && announce_plots)
                toc
                fprintf('Labels...\n'); % debug statement
            end

            % Apply labels
            if(plot_type == "raw") %title_format_raw =   '[%s|d0%d|%s] Raw: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
                title(sprintf('[%s|d0%d|%s] Height of Maximum Refl (mb) | %s',exp_name_plot,domain,member_string,timestamp_title),'FontSize',title_font_size); 
            elseif(plot_type == "dif")
                title(sprintf('[%s-%s|d0%d|%s] Dif: Heights of Maximum Refl (mb) | %s',exp_name_a,exp_name_b,domain,member_string,timestamp_title),'FontSize',title_font_size); 
            end
            xlabel('Longitude','FontSize',label_font_size);
            ylabel('Latitude','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);

            if(show_progress && announce_plots)
                toc
                fprintf('Saving...\n'); % debug statement
            end
            if(smooth_data)
                plot_filename = sprintf('%s_smoothed.png',extractBefore(plot_filename,".png"));
            end
            saveas(gcf,plot_filename); % Save as .png

            if(show_progress && announce_plots)
                toc
                fprintf('Making small...\n'); % debug statement
            end

            % Plot SMALL
            f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
            title(sprintf(title_format_small,upper(sprintf('%s-%s',data_src_a,data_src_b)),plot_type,data_type,timestamp_file),'FontSize',title_font_size); % Replace title
            if(show_progress && announce_plots)
                toc
                fprintf('Saving small...\n'); % debug statement
            end
            saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
            close('all');
            
            clearvars f h data_to_plot;
        
        end
        
        %% Clear big variables between members
        clearvars rmse_err bias data* -except data_storage* data_type data_name* data_gis_trim data_src* data_wrf_mean_ref;
        
    end

    clearvars ens_mem* obs_* member_sum* data_wrf_mean_ref;
    
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

% Interpolate a 3D variable to a 2D field at a specified pressure level
% MUST ALREADY HAVE BEEN RUN THROUGH ALIGN_DATA!
function data_interp = interpolate_to_pl(data_raw,dims_interp,pressure_3d,pressure_target)
    ylen = dims_interp(1);
    xlen = dims_interp(2);
    zlen = dims_interp(3);
    data_interp = zeros(ylen,xlen);
    for y = 1:ylen
        for x = 1:xlen
            p_idx = [0,0];
            for z = 1:zlen % Height DOES go up with Z here (so P goes down)
                if(pressure_3d(y,x,z) <= pressure_target)
                    if(pressure_3d(y,x,z) == pressure_target || z == 1)
                        p_idx(:) = z;
                    else
                        p_idx = [z-1 z];
                    end
                    break;
                end
            end
            if(~ismember(0,p_idx))
                data_interp(y,x) = linear_interpolate(pressure_3d(y,x,p_idx(1)),pressure_3d(y,x,p_idx(2)),pressure_target,data_raw(y,x,p_idx(1)),data_raw(y,x,p_idx(2)));
            else
                data_interp(y,x) = NaN;
            end
        end
    end
end


                