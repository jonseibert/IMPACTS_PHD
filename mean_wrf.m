% mean_wrf.m (v2.1.0)
%
% Purpose: Generate point-to-point (PTP) diagnostic comparison plots of 
% outputs from IMPACTS WRF to NEXRAD observations and each other. This 
% version computes outputs for all ensemble members, as well as their mean,
% and allows for comparison of multiple forecast variables.
% 
% Author(s): Jon Seibert
% Last updated: 9 September 2024
% 
% Inputs: [WRF Outputs].nc, [Observations].png, [Base Reflectivity].mat
% Outputs: [WRF Raw].png, [Obs Trim].png, [WRF-Obs Dif].png,
%   [WRF-WRF Analysis Dif].png, [WRF-WRF Background Dif].png, 
%   [WRF Analysis Increment Raw/Dif].png, [Error Table].txt.
%   
% Dependencies: reflmap.m, borders.m, find_nn_idx_irregular_exp.m, RMSE.m,
%       bilinear_interpolate_irregular_to_grid.m, wrf_3d_ref_to_base.m; 
%       latlon_[]_trim.mat (predefined domain boundaries), 
%       [WRF/GIS mask].mat, [Regridded WRF Data].mat (if available)
%
% Local functions:
%       - align_data: Performs rigid transformations to align input data
%           with standard grid, returned aligned subsection that falls within
%           the input borders
%           Parameters: (in_data,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx)
%       - interpolate_to_pl: Vertically interpolates 3D input data to the
%           desired height in pressure coordinates
%           Parameters: (data_raw,dims_interp,pressure_3d,pressure_target)
%
% NOTES:
%  - TIME_COUNT is not robust when crossing month boundaries.
%  - Not robust for use in southern hemisphere (lat < 0)
%
% TODO:
%  - Finish updating documentation
script_name = 'mean_wrf.m';


   
%% 0A. Script Controls
% Meant to be changed often

run_name = 'mean_test'; % Name of output folder

% Local or server execution mode
%local = 1; % true = running on local device, false = running on server

% Which experiments to use (0 is obs, only for B)
% To compare background and analysis, use the same experiment # for both
% To compare to observations, use 0 for experiment B
exp_choice_a = 2; 
exp_choice_b = 2;

data_type = "T"; % Working name of variable being analyzed (REFL, T, MSLP, GPH, VORT, WIND, Q, OMEGA)
% FOR THIS SCRIPT, DATA TYPE IS IRRELEVANT AND ONLY ONE EXP IS USED


% Date & time
% Range: 2020-02-07-1400 : 2020-02-07-1800; 14-00, 18-00
% 1400 is identical across all experiments
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 7;
start_hour = 14;
end_hour = 18;

domain = 2; % 2: d02, E coast region; 1: d01, eastern CONUS
subdomain = 0; %1 for E Coast Only subdomain, 2 for N only, 0 for all d02 included

num_members = 40; % Number of ensemble members to consider
%num_members = 2; % TESTING ONLY

% Vertical levels to use if is 3D
%pressure_levels = [460,500,700,850,900,925];
is_3d = 1; % Is the input data 3-dimensional?
pressure_target = 500; %mb
multi_height = 0; % Whether to loop over multiple target pressure values
pressure_target_list = [925,850,700]; % Overrides above pressure_target value if multi_height true

% Note whether transposition or compositing are necessary to align data with standard grid
transpose_data = true;
flip_data = true;
z_composite_data = false; % Unused if use_base_ref is true
smooth_data = 0; % If true AND the variable isn't already being smoothed, smooth it before plotting (gaussian)

% Smoothing settings
gauss_sd = 2.75; % Standard deviation to use for gaussian smoothing kernels
gauss_sd_pva = 4; % SD for PVA smoothing
gauss_width =  2*ceil(2*gauss_sd)+1; % Standard formula from SD for gaussian smoothing kernels
gauss_width_pva = 2*ceil(2*gauss_sd_pva)+1;

% Plot suite selections
remake_plots = 0; % If true, overrides existing outputs. If false, only plots 'new' files.
plot_mean_only = 1; % Suppresses individual member plots
plot_fit = 0; % 1: make plot, 0: do not
plot_dif = 0;
plot_trim = 0;
plot_raw = 1;
plot_sd = 0;
compute_error_table = 0; % Whether to compute and output error suite values

% Plot settings
use_base_ref = 0; % Whether to override REFL inputs with precalculated Base Reflectivity values
combo_gph_full = 1; % Whether to plot wind and vorticity on GPH plots
include_gph = 1; % Whether to plot GPH on wind plots
include_pva = 0; % Whether to also plot PVA when plotting VORT
label_contours = 1; % Whether to label contour plots with their values

sd_max = 3; % Observed maximum spatial SD value
thinning_factor = 15; % How many of the wind arrows to plot (1/this)

% Colorbar settings: overridden by some choices of data variables
limit_raw_colorbar = 0;
contour_L_coloring = 0; % If true, use experimental contour line coloring

show_plots = 0; % 0: Suppress plot popups (REQUIRES MATLAB RESTART TO UNDO)
show_progress = 1; % If true, prints out progress bar/debug messages
announce_plots = 0; % If true, prints out individual plot messages

ecmwf_filename = 'ecmwf_mslp_d01_20200207.nc';

%% 0B. General Settings
% Should not be changed without good reason

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
    exp_path_1 = 'AIR'; % TESTING
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

addpath(path_to_extra_code); % Make sure Matlab can see the addon code

%ecmwf_filename = 'ecmwf_mslp_d01_20200207.nc';

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
% These are overridden by latlon_wrf_trim.mat
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

% Subdomain settings
bound_lon_d1 = [-78,-78,-74.5,e_lim,e_lim,-75,-78];
bound_lat_d1 = [40,44,n_lim,n_lim,41.5,38,40];
%inpolygon([-74,50],[42,50],bound_lon,bound_lat);
bound_lon_d2 = [-77.5,-77.5,-74,e_lim,e_lim,-73.5,-77.5];
bound_lat_d2 = [42,44,n_lim,n_lim,42.5,42,42];

process_all_members = true; % Whether to compute error table values for each ensemble member
use_auto_run_name = false; % Whether to override manual run name

tight_L = true; % Compress ensemble L spread plot for readability

% Constants
Re = 6.3781e6; % Radius pf Earth in m
g = 9.80665; % m/s^2 standard gravity at sea level
R = 287.0600676; % J/kg/K, specific gas constant

input_format_base_refl ='%s/%s_%s_d%02.f_%s_%s_%s.mat'; % output path, data source, plot type, domain, member #, data type, datetime
datetime_format_file =  '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

%% 0C. Settings- WRF Data (2020 Case)

% Basic details
%num_members = 40;
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
% Set WRF NetCDF variable names for data retrieval
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
data_name_refl = "REFL_10CM";
        
%% 1. Preliminary Processing & Setup

if(~show_plots)
    set(groot,'DefaultFigureVisible','off') % Turn off figure popups for local
end

% Announce operating mode
fprintf('Starting %s in %s mode.\n',script_name,upper(mode));
tic; % Start script timer

% Check for correct execution mode
working_dir = strrep(pwd,"\","/");
if(working_dir ~= path_to_code)
    error("Operating mode does not match execution directory.")
end

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

% Define sample file to read latlon data from
timestamp_file = sprintf(datetime_format_file,year,month,day,hour);

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

% Get dimensions
dimensions_wrf = size(lon_wrf);
dims = size(lon_wrf);
ylen_wrf = dimensions_wrf(1);
xlen_wrf = dimensions_wrf(2);
zlen_wrf = 50;
clearvars dimensions_wrf lon_wrf_trim lat_wrf_trim;


%%

exp_choice = exp_choice_a;
exp_name = exp_name_a;
input_path = input_path_a;

%% 3B. For each timestamp being analyzed:
for time_idx = 1:time_count
   
    % Set timestamps and experiment names
    timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
    timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
    timestamp_file_background = sprintf(datetime_format_file,year,month,day,hour-1);
    data_src = lower(exp_name);
    
    % Define holding objects
    ens_mem_u =   zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_v =   zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_w =   zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_gph = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_p =   zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_t =   zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_omega = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_refl = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);
    ens_mem_mslp = zeros(ylen_wrf,xlen_wrf,num_members);

    ens_mem_mh = zeros(ylen_wrf,xlen_wrf,zlen_wrf,num_members);

    L_val = zeros(num_members+1);
    L_dist = zeros(num_members+1);
    L_loc = zeros(num_members+1,2); % member, lon/lat
    
    % Read in ECMWF obs data to compare with
    % Oriented as [lon,lat,time]
    obs_full_filename = sprintf('%s/%s',input_path_ecmwf,ecmwf_filename);
    obs_lon_raw = ncread(obs_full_filename,'longitude');
    obs_lat_raw = ncread(obs_full_filename,'latitude');
    obs_time = ncread(obs_full_filename,'time');
    obs_mslp_raw = ncread(obs_full_filename,'msl');
    [obs_lat_full,obs_lon_full] = meshgrid(obs_lat_raw,obs_lon_raw);
    obs_lat_full = obs_lat_full';
    obs_lon_full = obs_lon_full';
    %[obs_lon_grid,obs_lat_grid] = meshgrid(obs_lon,obs_lat);
    % obs_mslp = obs_mslp(:,:,hour+1)/100; % Convert from Pa to hPa/mb
    obs_mslp = obs_mslp_raw(85:122,13:53,hour+1)/100; % Convert from Pa to hPa/mb
    obs_mslp_full = obs_mslp_raw(:,:,hour+1)/100;
    obs_lon = obs_lon_raw(85:122);
    obs_lat = obs_lat_raw(13:53);
    [obs_lat,obs_lon] = meshgrid(obs_lat,obs_lon);
    obs_mslp_full = obs_mslp_full';
    obs_lon = obs_lon';
    obs_lat = obs_lat';
    obs_mslp = obs_mslp';
    % Left idx = 85
    % Right idx = 122
    % Top_idx = 13
    % Bot_idx = 53
    [L_val_o,p_i] = min(obs_mslp,[],'all','linear');
    p_min_row_o = mod(p_i,size(obs_mslp,1));
    p_min_col_o = 1 + floor(p_i/size(obs_mslp,1));
    L_lon_o = obs_lon(p_min_row_o,p_min_col_o);
    L_lat_o = obs_lat(p_min_row_o,p_min_col_o);

%% 3C: Member loop
    for member_idx = 1:(num_members+1)
        
        % If all members have been computed, take mean
        if(member_idx == (num_members+1)) 
            member_string = 'mean';

            data_u = mean(ens_mem_u,4,"omitnan");
            data_v = mean(ens_mem_v,4,"omitnan");
            data_w = mean(ens_mem_w,4,"omitnan");
            data_gph = mean(ens_mem_gph,4,"omitnan");
            data_p = mean(ens_mem_p,4,"omitnan");
            data_t = mean(ens_mem_t,4,"omitnan");
            data_omega = mean(ens_mem_omega,4,"omitnan");
            data_refl = mean(ens_mem_refl,4,"omitnan");
            model_height = mean(ens_mem_mh,4,"omitnan");

            data_mslp = mean(ens_mem_mslp,3,"omitnan");

            L_loc(num_members+1,:) = mean(L_loc(1:num_members,:));
            L_val(num_members+1) = mean(L_val(1:num_members));
            L_dist(num_members+1) = mean(L_dist(1:num_members));

            save(sprintf('%s/mean_storage_%s_%s.mat',intermediate_path,exp_name,timestamp_file),'data_refl','data_u','data_v','data_w','data_gph','data_p','data_t','data_mslp','data_omega','model_height','terrain_height','L_loc','L_val');

        else % Otherwise, do the next member
        
            member_string = sprintf('%03.f',member_idx);
            
            if(show_progress)
                toc
                fprintf('Processing: %s | Member %s\n',timestamp_title,member_string);
            end
            
            % Define full filepaths- exps 3,4,5 have dif format
            if(ismember(exp_choice,alt_input_exps))
                filename = sprintf(filename_format_free,domain,year,month,day,hour);
                full_filepath = sprintf('%s/%s/%s',input_path,member_string,filename);
            else
                filename = sprintf(filename_format_analysis,domain,member_idx);       
                full_filepath = sprintf('%s/%s/%s',input_path,timestamp_file,filename);
            end
    
            lon_wrf_raw = ncread(full_filepath,lon_name_wrf);
            lat_wrf_raw = ncread(full_filepath,lat_name_wrf);
            dims = size(lon_wrf_raw);
            xlen_wrf_raw = dims(1);
            ylen_wrf_raw = dims(2);
            %zlen_wrf_raw = dims(3);
            
            % Read in all necessary WRF fields
            data_pb = ncread(full_filepath,data_name_pres); % Base pressure (Pa)
            data_pp = ncread(full_filepath,data_name_pres_prime); % Perturbation pressure (Pa)
            data_ptp = ncread(full_filepath,data_name_temp); % Perturbation potential temperature (K)
            data_p = data_pb + data_pp; % 3D Pressure field (Pa)
            data_pt = data_ptp + 300; % Potential temperature (K)
            data_t = data_pt.*((data_p./100000).^(0.287)); % 3D Temperature field (K)
            data_gpb_raw = ncread(full_filepath,data_name_gp); % Base geopotential (m^2/s^2)
            data_gpp_raw = ncread(full_filepath,data_name_gp_prime); % Perturbation geopotential (m^2/s^2)
            data_p = data_p./100; % Convert P from Pascals (Pa) to Millibars (mb)
            data_t2m = ncread(full_filepath,data_name_t2m); % 2m Temperature (K)
            data_q = ncread(full_filepath,data_name_q); % Water vapor (kg/kg)
            data_u_raw = ncread(full_filepath,data_name_u); % Wind (m/s)
            data_v_raw = ncread(full_filepath,data_name_v); % Wind (m/s)
            data_w_raw = ncread(full_filepath,data_name_w); % Wind (m/s)
            data_refl = ncread(full_filepath,data_name_refl);

            terrain_height = ncread(full_filepath,data_name_terrain_height); % Height of physical terrain (m)
            
            % Define holding containers
            data_u = zeros(size(data_t));
            data_v = zeros(size(data_t));
            data_w = zeros(size(data_t));
            data_gpb = zeros(size(data_t));
            data_gpp = zeros(size(data_t));
            
            % Fix staggered WRF grids for wind, geopotential
            for y = 1:ylen_wrf_raw
                for x = 1:xlen_wrf_raw % NOT ALIGNED TO (Y,X) YET
                    data_u(x,y,:) = (data_u_raw(x,y,:) + data_u_raw(x+1,y,:))./2;
                    data_v(x,y,:) = (data_v_raw(x,y,:) + data_v_raw(x,y+1,:))./2;
                end
            end
            
            for z = 1:zlen_wrf
                data_gpb(:,:,z) = (data_gpb_raw(:,:,z) + data_gpb_raw(:,:,z+1))./2;
                data_gpp(:,:,z) = (data_gpp_raw(:,:,z) + data_gpp_raw(:,:,z+1))./2;
                data_w(:,:,z) = (data_w_raw(:,:,z) + data_w_raw(:,:,z+1))./2;
            end
            
            % Compute geopotential height
            data_gp = data_gpb + data_gpp; %m^2/s^2
            data_gph = data_gp./g; %m
            data_gph = data_gph/10; % Convert from m to dam (decameters)
            
            % Compute MSLP: P1 = P2*e^(((g/(R*Tv))*(z2-z1))
            z1 = 0; % m
            h1 = 1;
            h2 = 10;
            
            P2 = data_p(:,:,h2); % Lowest level of P
            Tv_1 = data_t(:,:,h1).*(1 + 0.608*data_q(:,:,h1));
            Tv_2 = data_t(:,:,h2).*(1 + 0.608*data_q(:,:,h2));
            Tv = (Tv_1 + Tv_2)/2;
            z2 = (data_gp(:,:,h2).*Re)./((g*Re) - data_gp(:,:,h2)); % Height above mean sea level
            data_mslp = P2.*exp((g./(R.*Tv)).*(z2-z1)); % MSLP (mb)
            
            % Compute statistics relative to observed sfc Low
            [p_min,p_i] = min(data_mslp,[],'all','linear');
            p_min_row = mod(p_i,size(data_mslp,1));
            p_min_col = 1 + floor(p_i/size(data_mslp,1));
            L_val(member_idx) = p_min;
            L_loc(member_idx,1) = lon_wrf_raw(p_min_row,p_min_col);
            L_loc(member_idx,2) = lat_wrf_raw(p_min_row,p_min_col);
            L_dist(member_idx) = haversine_distance(lon_wrf_raw(p_min_row,p_min_col),lat_wrf_raw(p_min_row,p_min_col),L_lon_o,L_lat_o);
            
            model_height = ((data_gpb+data_gpp)./g) - terrain_height; % model height above the *ground* (m)
            
            data_t = data_t - 273.15; % Convert T from K to C
            
            % Align all data fields with customized grid
            data_p = align_data(data_p,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                
            data_u = align_data(data_u,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_v = align_data(data_v,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_w = align_data(data_w,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_gph = align_data(data_gph,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_t = align_data(data_t,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_mslp = align_data(data_mslp,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            data_refl = align_data(data_refl,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
            model_height = align_data(model_height,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);

            data_omega = zeros(size(data_t));

            % Also calculate omega (uplift, dp/dt, vertical wind in pressure coords)
            for z = 1:zlen_wrf
                if(z == 1)
                    dp = (data_p(:,:,z+1) - data_p(:,:,z));
                    dz = (model_height(:,:,z+1) - model_height(:,:,z));
                elseif(z == zlen_wrf)
                    dp = (data_p(:,:,z) - data_p(:,:,z-1));
                    dz = (model_height(:,:,z) - model_height(:,:,z-1));
                else
                    dp = (data_p(:,:,z+1) - data_p(:,:,z-1));
                    dz = (model_height(:,:,z+1) - model_height(:,:,z-1));
                end
                data_omega(:,:,z) = data_w(:,:,z).*(dp./dz); %dp/dt = (dz/dt)*(dp/dz)
            end
                
            dims = [ylen_wrf xlen_wrf zlen_wrf];

            % Store for mean
            ens_mem_u(:,:,:,member_idx) = data_u;
            ens_mem_v(:,:,:,member_idx) = data_v;
            ens_mem_w(:,:,:,member_idx) = data_w;
            ens_mem_gph(:,:,:,member_idx) = data_gph;
            ens_mem_p(:,:,:,member_idx) = data_p;
            ens_mem_t(:,:,:,member_idx) = data_t;
            ens_mem_omega(:,:,:,member_idx) = data_omega;
            ens_mem_refl(:,:,:,member_idx) = data_refl;
            ens_mem_mh(:,:,:,member_idx) = model_height;
            ens_mem_mslp(:,:,member_idx) = data_mslp;
        end
        
        % Compute thin indices for wind barbs
        thin_y_idxs = 1:thinning_factor:ylen_wrf;
        thin_x_idxs = 1:thinning_factor:xlen_wrf;
    end

    clearvars ens_mem* obs_* member_sum* data_wrf_mean_ref;
    clearvars rmse_err bias data* -except data_storage* data_type data_name* data_gis_trim data_src* data_wrf_mean_ref;


    %% 3J. Time increment
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

toc
fprintf('Done.\n');


%% Appendix: Local functions

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
