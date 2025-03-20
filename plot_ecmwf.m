% ptp_analysis_ens_full.m (v2.1.9)
%
% Purpose: Generate point-to-point (PTP) diagnostic comparison plots of 
% outputs from IMPACTS WRF to NEXRAD observations and each other. This 
% version computes outputs for all ensemble members, as well as their mean,
% and allows for comparison of multiple forecast variables.
% 
% Author(s): Jon Seibert
% Last updated: 27 Nov 2024
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
script_name = 'ptp_analysis_ens_full.m';
   
%% 0A. Script Controls
% Meant to be changed often

run_name = 'burst_testing'; % Name of output folder

% Which experiments to use (0 is obs, only for B)
% To compare background and analysis, use the same experiment # for both
% To compare to observations, use 0 for experiment B
exp_choice_a = 3; 
exp_choice_b = 4;

data_type = "SNOWH"; % Working name of variable being analyzed (REFL, T, MSLP, GPH, VORT, WIND, Q, OMEGA) [+T2M,SNOWH(depth),SNOWNC(grid total per timestep)]
use_base_refl = 0; % Whether to override REFL inputs with precalculated Base Reflectivity values

% Date & time
% Range: 2020-02-07-1400 : 2020-02-07-1800; 14-00, 18-00
% 1400 is identical across all experiments
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 8;
start_hour = 18;
end_hour = 0;

domain = 2; % 2: d02, E coast region; 1: d01, eastern CONUS
subdomain = 0; %1 for E Coast Only subdomain, 2 for N only, 0 for all d02 included

% Vertical levels to use if is 3D
is_3d = 0; % Is the input data 3-dimensional?
pressure_target = 1000; %mb [0 = vc, 1000 = sfc]

z_composite_data = 0; % Unused if use_base_ref is true
multi_height = 0; % Whether to loop over multiple target pressure values
pressure_target_list = [1000,925]; % [925,850,700]; Overrides above pressure_target value if multi_height true

% Plot suite selections
smooth_data = 0; % If true AND the variable isn't already being smoothed, smooth it before plotting (gaussian)
remake_wrf_regrid = 0; % Whether to override existing saved WRF regrids
compute_error_table = 0; % Whether to compute and output error suite values

remake_plots = 1; % If true, overrides existing outputs. If false, only plots 'new' files.
plot_mean_only = 1; % Suppresses individual member plots
plot_raw = 1; % Base values
plot_dif = 1; % Difference plots
plot_trim = 0; % Base values, but for obs

% Progress messages
show_plots = 0; % 0: Suppress plot popups [REQUIRES MATLAB RESTART TO UNDO]
show_progress = 1; % If true, prints out progress bar/debug messages
announce_plots = 0; % If true, prints out individual plot messages

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
    exp_path_1 = 'AIR';
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

if(local)
    num_members = 2;
else
    num_members = 40; % Number of ensemble members to consider
    %num_members = 2; % TESTING ONLY
end

ecmwf_filename = 'ecmwf_mslp_d01_20200207.nc'; % [TESTING]
ecmwf_filename_qt_a = 'ecmwf_qt_d02_20200207.nc';
ecmwf_filename_qt_b = 'ecmwf_qt_d02_20200208.nc';
ecmwf_filename_t2m = 'ecmwf_t2m_d02_20200207.nc';

% Experiment names
exp_name_1 = 'AIR';
exp_name_2 = 'CONV';
exp_name_3 = 'AIR-F';
exp_name_4 = 'CONV-F';
exp_name_5 = 'NODA';

alt_input_exps = [3 4 5]; % Input is set up differently for these experiments

% Note whether transposition or compositing are necessary to align data with standard grid
transpose_data = 1;
flip_data = 1;
%z_composite_data = 0; % Unused if use_base_ref is true

hour_step = 1; % Hour increment size

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 925;
fig_height = 900;
%fig_width_small = 350;
%fig_height_small = 350;
fig_width_small = 525;
fig_height_small = 500;

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

plot_sd = 0; % Standard deviations [NON-FUNCTIONAL]
plot_original_wrf = 0; % Base values, but without masking and interpolation
plot_eqt = 0; % ECMWF ERA-5 pressure

% Smoothing settings
gauss_sd = 2.75; % Standard deviation to use for gaussian smoothing kernels
gauss_sd_pva = 4; % SD for PVA smoothing
gauss_width =  2*ceil(2*gauss_sd)+1; % Standard formula from SD for gaussian smoothing kernels
gauss_width_pva = 2*ceil(2*gauss_sd_pva)+1;

% Plot settings
combo_gph_full = 1; % Whether to plot wind and vorticity on GPH plots
include_gph = 1; % Whether to plot GPH on wind plots
include_pva = 0; % Whether to also plot PVA when plotting VORT
label_contours = 1; % Whether to label contour plots with their values

sd_max = 3; % Observed maximum spatial SD value
thinning_factor = 15; % How many of the wind arrows to plot (1/this)

% Filler values
nodata_val = NaN;
nodata_val_gis = NaN;

edge_buffer = 42; % Minimum number of points to search from theoretical min in each direction (WARNING: Artifacting may occur at values lower than 36.)

ets_threshold_list = [20,35]; % dBZ threshold for Equitable Threat Scores (ETS)
ets_radius = 9000; % Neighborhood radius for "hits" in threat score calculations

error_dec = 3; % Decimal places in error table

% Subdomain settings
bound_lon_d1 = [-78,-78,-74.5,e_lim,e_lim,-75,-78];
bound_lat_d1 = [40,44,n_lim,n_lim,41.5,38,40];
bound_lon_d2 = [-77.5,-77.5,-74,e_lim,e_lim,-73.5,-77.5];
bound_lat_d2 = [42,44,n_lim,n_lim,42.5,42,42];

process_all_members = true; % Whether to compute error table values for each ensemble member
use_auto_run_name = false; % Whether to override manual run name
tight_L = true; % Compress ensemble L spread plot for readability

% More colorbar settings
clim_lower = -30; % Colorbar limits (REFL)
clim_upper = 75;
clim_lower_dif = -10; 
clim_upper_dif = 10;
clim_lower_vort = 0;
clim_upper_vort = 50;
clim_upper_pva = 0.20;
clim_lower_pva = 0;

% Constants
Re = 6.3781e6; % Radius pf Earth in m
g = 9.80665; % m/s^2 standard gravity at sea level
R = 287.0600676; % J/kg/K, specific gas constant

% Colorbar settings: overridden by some choices of data variables
limit_raw_colorbar = 1; % Whether to apply specified min and max colorbar values
jet_blue_percent = 0.2; % Cut off coldest (Darkest blues) _% of jet colorbar for ease of analysis

% Custom colormaps
%[100 0 0] -> [255 115 115] (reds)
reds = [100:((255-100)/(num_members-1)):255; 0:((115-0)/(num_members-1)):115; 0:((115-0)/(num_members-1)):115];
reds = (reds./255)';
%[0 100 0] -> [115, 255, 115] (greens)
greens = [0:((115-0)/(num_members-1)):115; 100:((255-100)/(num_members-1)):255; 0:((115-0)/(num_members-1)):115];
greens = (greens./255)';
%[0 43 100] -> [115, 175, 255] (blues)
blues = [0:((115-0)/(num_members-1)):115; 43:((175-43)/(num_members-1)):175; 100:((255-100)/(num_members-1)):255];
blues = (blues./255)';
%[18 0 100] -> [140 115 255] (purples)
purples = [18:((140-18)/(num_members-1)):140; 0:((100-0)/(num_members-1)):100; 100:((255-100)/(num_members-1)):255];
purples = (purples./255)';

% Chop down custom version of "jet" colormap
jet_dims = size(jet);
jet_lims = [round(jet_dims(1)*jet_blue_percent) jet_dims(1)];
jet_modded = jet;
jet_modded = jet_modded(jet_lims(1):jet_lims(2),:);

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

%% 0C. Settings- WRF Data (2020 Case)

% Basic details
%num_members = 40; [Moved to section 0A]
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
data_name_t2m = "T2"; % 2m Temperature


%%



ecmwf_filename = 'ecmwf_mslp_d01_20200207.nc'; % [TESTING]
ecmwf_filename_qt_a = 'ecmwf_qt_d02_20200207.nc';
ecmwf_filename_qt_b = 'ecmwf_qt_d02_20200208.nc';




L_val = zeros(num_members+1,2);
        L_dist = zeros(num_members+1,2);
        L_loc = zeros(num_members+1,2,2); % member, lon/lat, A/B
        
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
    
        if(plot_eqt)
    
            % Read in ECMWF obs data to compare with
            % Oriented as [lon,lat,time]
            eqt_full_filename = sprintf('%s/%s',input_path_ecmwf,ecmwf_filename_qt_a);
            %eqt_full_filename_b = sprintf('%s/%s',input_path_ecmwf,ecmwf_filename_qt_b);
            eqt_lon_raw = ncread(eqt_full_filename,'longitude');
            eqt_lat_raw = ncread(eqt_full_filename,'latitude');
            eqt_time_sec = ncread(eqt_full_filename,'valid_time');
            eqt_time = 14:23;
            eqt_time = eqt_time'; % 2
            eqt_q = ncread(eqt_full_filename,'q');
            eqt_t = ncread(eqt_full_filename,'t');
            [eqt_lat_full,eqt_lon_full] = meshgrid(eqt_lat_raw,eqt_lon_raw);
            eqt_lat = eqt_lat_full;
            eqt_lon = eqt_lon_full;
            eqt_p = ncread(eqt_full_filename,'pressure_level'); % 9
            eqt_time_choice = 2;
            eqt_p_choice = 9;
        
        end