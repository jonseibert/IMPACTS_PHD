% ptp_analysis_ens_full.m (v2.2.0)
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

run_name = 'burst_testing/snowcheck'; % Name of output folder

% Which experiments to use (0 is obs, only for B)
% To compare background and analysis, use the same experiment # for both
% To compare to observations, use 0 for experiment B
exp_choice_a = 2; 
exp_choice_b = 0;

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
end_day = 7;
start_hour = 4;
end_hour = 17;

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
plot_mean_only = 0; % Suppresses individual member plots
plot_first_and_mean_only = 1;
plot_raw = 1; % Base values
plot_dif = 1; % Difference plots
plot_trim = 0; % Base values, but for obs

% Progress messages
show_plots = 0; % 0: Suppress plot popups [REQUIRES MATLAB RESTART TO UNDO]
show_progress = 1; % If true, prints out progress bar/debug messages
announce_plots = 0; % If true, prints out individual plot messages

transpose_obs = 1;
flip_obs = 0;

s_box = 43;
n_box = 45;
w_box = -76;
e_box = -73;
box_only = 0;
apply_mask = 0;
make_mask = 0;

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
    input_path_snow = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/nws_data';
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
    input_path_snow = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/NOAA';
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

noaa_snow_24h_filename = 'sfav2_CONUS_24h_2020020800.nc';

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
nodata_val_obs = NaN;

edge_buffer = 42; % Minimum number of points to search from theoretical min in each direction (WARNING: Artifacting may occur at values lower than 36.)

ets_threshold_list = [20,35]; % dBZ threshold for Equitable Threat Scores (ETS)
ets_radius = 9000; % Neighborhood radius for "hits" in threat score calculations

error_dec = 3; % Decimal places in error table

% Subdomain settings
bound_lon_d1 = [-78,-78,-74.5,e_lim,e_lim,-75,-78];
bound_lat_d1 = [40,44,n_lim,n_lim,41.5,38,40];
bound_lon_d2 = [-77.5,-77.5,-74,e_lim,e_lim,-73.5,-77.5];
bound_lat_d2 = [42,44,n_lim,n_lim,42.5,42,42];
bound_lon_box = [-76.5,-76.5,-73,-73,-76.5];
bound_lat_box = [42,45,45,42,42];

s_box = 43;
n_box = 45;
w_box = -76;
e_box = -73;

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

cm_to_in = 0.3937008;
use_snow_inches = true;

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
title_format_fit =   '[%s|d0%d|%s] %s: %s(%s)|%s'; % exp name 1, domain, member #/mean, background/analysis/obs (bao), bao, plot_type, data type, units, timestamp
title_format_dif =   '[%s-%s|d0%d|%s] %s: %s(%s) |%s\nRMSE = %0.3f, Bias = %0.3f'; % exp name 1, exp name 2, domain, member #/mean, plot_type, data type, units, height, timestamp, err, bias
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

%% 0D. Settings- NEXRAD Image Composites via GIS, Iowa Env. Mesonet

data_src_gis = "nex";
file_format_gis = '%s/n0q_%04.f%02.f%02.f%02.f%02.f.png';

lat_gis_raw = 49.9975:-0.005:23.0025;
lon_gis_raw = -125.9975:0.005:-65.0025;
[lon_gis,lat_gis] = meshgrid(lon_gis_raw,lat_gis_raw);

% NOTES: Comes in image format, does not contain own latlon data
% Conversion formula: data_gis = double(imread())*0.5 - 32.5;
        
%% 1. Preliminary Processing & Setup

if(~show_plots)
    set(groot,'DefaultFigureVisible','off') % Turn off figure popups for local
end

% Announce operating mode
fprintf('#--------------------------------------------------#\n');
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
        input_path_b = input_path_snow;
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
    if(remake_wrf_regrid)
        fprintf('WARNING: Overriding existing WRF regrids.\n');
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

if(multi_height)
    num_heights = length(pressure_target_list);
else
    num_heights = 1;
end

% Create output folders
% Directory structure:
% output/ptp/(run_name)/(variable_name)/(plot_size)
output_path_large = sprintf('%s/%s/%s/%s/%s',output_path_base,run_name,upper(data_type),'large');
output_path_small = sprintf('%s/%s/%s/%s/%s',output_path_base,run_name,upper(data_type),'small');

% If output folders do not exist, create them
if(~isfolder(output_path_large) || ~isfolder(output_path_small))
    mkdir(output_path_large);
    mkdir(output_path_small);
end

% Define table to hold all the error numbers and later write to a file
err_rowNames = 1:(num_members+6);
err_rowNames = string(err_rowNames);
err_rowNames(num_members+1) = "EMean";
err_rowNames(num_members+2) = "Median";
err_rowNames(num_members+3) = "MeanEV";
err_rowNames(num_members+4) = "Min";
err_rowNames(num_members+5) = "Max";
err_rowNames(num_members+6) = "StdEV";

% Set appropriate values for the input data type; define units
% REFL_10CM/REFL (dBZ), T/T (deg C), P/MSLP (mb), GPH/GPH(m),
% VORT/VORT(10^-5 s^-1), U-V-W/WIND(kt)
contour_steps = 6; % P: 4, T: 2, GPH: 6 (default to 6 for GPH)
plot_fl_np = false; % Whether to plot freezing line contours
switch data_type
    case "GPH"
        units = "dam";
        data_name = ""; % WRF Referential
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
    case "T2M" % 2m Temperature
        units = sprintf('%sC',char(176));
        data_name = "T2";
        contour_steps = 2; % P: 4, T: 2, GPH: 6
        use_contours = 1;
        plot_fl_np = 0;
    case "MSLP"
        units = "mb";
        data_name = "P"; 
        contour_steps = 4; % P: 4, T: 2, GPH: 6
        use_contours = 1;
    case "REFL"
        units = "dBZ";
        data_name = "REFL_10CM"; 
        use_contours = 0;
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
    case "SNOWNC"
        units = "mm";
        data_name = "SNOWNC";
        clim_lower = 0;
        %clim_upper = 360; % ~14 inches in mm
        clim_upper = 8;
        clim_lower_dif = -8;
        clim_upper_dif = 8;
        use_contours = 0;
        limit_raw_colorbar = 1;
    case "SNOWH"
        data_name = "SNOWH";
        units = "cm"; % Converted from m
        clim_lower = 0;
        clim_upper = 7e1;
        clim_lower_dif = -40;
        clim_upper_dif = 40;
        use_contours = 0;
        limit_raw_colorbar = 1;
        if(use_snow_inches)
            units = "in";
            clim_lower = 0;
            clim_upper = clim_upper*cm_to_in;
            clim_lower_dif = clim_lower_dif*cm_to_in;
            clim_upper_dif = clim_upper_dif*cm_to_in;
        end
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

% Create cross-member storage objects
data_storage_a = zeros(ylen_wrf,xlen_wrf,num_members,num_heights);
data_storage_b = zeros(ylen_wrf,xlen_wrf,num_members,num_heights);

% Split into two sections: compare to obs, compare to another WRF output

%% 2A. Main Loop: If comparing WRF to OBS
%if(compare_to_obs)

    input_path_wrf = input_path_a;    
    exp_choice = exp_choice_a;
    exp_name = exp_name_a;
    
    
    if(show_progress)
        toc
        fprintf('Comparison: %s-obs, %s\n',exp_name_a,upper(data_type));
    end
    
    if(upper(data_type) ~= "SNOWH")
        error(sprintf("Data Type: %s. Only snow depth observations are supported.",upper(data_type)));
    end
    
    %noaa_snow_24h_filename = 'sfav2_CONUS_24h_2020020800.nc';
    noaa_24h_snow = ncread(sprintf('%s/%s',input_path_snow,noaa_snow_24h_filename),'Data'); %m
    lon_obs_raw = ncread(sprintf('%s/%s',input_path_snow,noaa_snow_24h_filename),'lon');
    lat_obs_raw = ncread(sprintf('%s/%s',input_path_snow,noaa_snow_24h_filename),'lat');
    %latb_obs = ncread(sprintf('%s/%s',input_path_snow,noaa_snow_24h_filename),'lat_bounds');
    [lon_obs, lat_obs] = meshgrid(lon_obs_raw,lat_obs_raw); % latlon natural alignment is good, data needs transpose (+ flip?)

    if(transpose_obs)
        noaa_24h_snow = noaa_24h_snow';
    end
    if(flip_obs)
        noaa_24h_snow = flip(noaa_24h_snow);
    end
    

    if(use_snow_inches)
        noaa_24h_snow = noaa_24h_snow.*100*cm_to_in; %in
    else
        noaa_24h_snow = noaa_24h_snow.*100; % m-> cm
    end

    
    % Find shared area between WRF and GIS latlon coverage
    w_bound = max(min(min(lon_wrf)),min(min(lon_obs)));
    e_bound = min(max(max(lon_wrf)),max(max(lon_obs)));
    s_bound = max(min(min(lat_wrf)),min(min(lat_obs)));
    n_bound = min(max(max(lat_wrf)),max(max(lat_obs)));
    
    % IN THIS CASE, GIS ENTIRELY CONTAINS WRF: FIND EDGE POINTS IN GIS
    % LATLON THAT FALL JUST INSIDE THE WRF BOUNDARIES
    % Trim GIS to those points, regrid WRF to the trimmed GIS
    
    for x=1:length(lon_obs_raw)
        if(lon_obs_raw(x) > w_bound && lon_obs_raw(x) < e_bound)
            w_idx = x;
            break;
        end
    end
    
    for x=length(lon_obs_raw):-1:1
        if(lon_obs_raw(x) > w_bound && lon_obs_raw(x) < e_bound)
            e_idx = x;
            break;
        end
    end
    
    for y = 1:length(lat_obs_raw)
        if(lat_obs_raw(y) > s_bound && lat_obs_raw(y) < n_bound)
            s_idx = y;
            break;
        end
    end
    
    for y = length(lat_obs_raw):-1:1
        if(lat_obs_raw(y) > s_bound && lat_obs_raw(y) < n_bound)
            n_idx = y;
            break;
        end
    end
    
    lon_obs_trim_raw = lon_obs_raw(w_idx:1:e_idx);
    lat_obs_trim_raw = lat_obs_raw(s_idx:1:n_idx);
    data_obs = noaa_24h_snow(s_idx:1:n_idx,w_idx:1:e_idx);
    
    [lon_obs_trim,lat_obs_trim] = meshgrid(lon_obs_trim_raw,lat_obs_trim_raw);
    clearvars *obs*raw x y;
    
    % Define empty grid to fit wrf data into
    dimensions_obs_trim = size(lon_obs_trim);
    ylen_obs = dimensions_obs_trim(1);
    xlen_obs = dimensions_obs_trim(2);
    
    % Determine short grid box to search in for neighborhood version of ETS/CSI
    x_center = round(xlen_obs/2);
    y_center = round(ylen_obs/2);
    grid_dist_radius = 0;
    for x = x_center:xlen_obs
        d = haversine_distance(lon_obs_trim(x,y_center),lat_obs_trim(x,y_center),lon_obs_trim(x_center,y_center),lat_obs_trim(x_center,y_center));
        if(d > ets_radius)
            grid_dist_radius = abs(x_center - x);
            break;
        end
    end
    clearvars x_center y_center;
    
    %% 2B. For each timestamp being analyzed:
    %for time_idx = 1:time_count
       
        timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
        timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
        data_src = lower(exp_name);
        
        RMSE_list = zeros(size(err_rowNames));
        bias_list = zeros(size(err_rowNames));
        SDE_list = zeros(size(err_rowNames));
        SD_list = zeros(size(err_rowNames));
        CSI_list_a = zeros(size(err_rowNames));
        CSI_list_b = zeros(size(err_rowNames));
        ETS_list_a = zeros(size(err_rowNames));
        ETS_list_b = zeros(size(err_rowNames));
        brier_list_a = zeros(size(err_rowNames));
        brier_list_b = zeros(size(err_rowNames));
        sd_grid = zeros(ylen_wrf,xlen_wrf);
        sde_grid = zeros(ylen_obs,xlen_obs);
        
        % Regrid WRF data onto trimmed obs grid
        % Brute force bilinear interp to make the WRF data match the shape of the obs data
        % NOTE: Assumes that the WRF data fully encompasses the trimmed obs region, so that
        % edge cases can work as normal!
    
        member_sum = zeros(ylen_wrf,xlen_wrf);
        member_sum_bin_a = zeros(ylen_obs,xlen_obs);
        member_sum_bin_b = zeros(ylen_obs,xlen_obs); 
        
        data_3z = zeros(ylen_wrf,xlen_wrf);
        data_0z = zeros(ylen_wrf,xlen_wrf);

        clearvars n_idx s_idx w_idx e_idx;
        load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'n_idx','s_idx','e_idx','w_idx'); % Load in preestablished latlon values
        

        for time_idx = 1:2

            if(time_idx == 1)
                hour = start_hour;
            else
                hour = end_hour;
                day = end_day;
            end

            timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
            timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
            timestamp_file_background = sprintf(datetime_format_file,year,month,day,hour-1);
    
            % Read in each member and compute the mean
            for member_idx = 1:(num_members+1)
            %% 2C: Member loop
                if(member_idx == (num_members+1))
                    data_wrf = member_sum/num_members;
                    member_string = 'mean';
                    if(show_progress)
                        toc
                        fprintf('Processing: %s | Ensemble Mean\n',timestamp_title);
                    end
                else
                
                    member_string = sprintf('%03.f',member_idx);
                    if(show_progress)
                        toc
                        fprintf('Processing: %s | Member %s\n',timestamp_title,member_string);
                    end
        
                    if(time_idx == 1)
                        if(local)
                            filename = sprintf(filename_format_background_alt,domain,year,month,day,hour,member_idx); % Local files can't have : in filename
                        else
                            filename = sprintf(filename_format_background,domain,year,month,day,hour,member_idx);     
                        end     
                        input_path_alt = sprintf('%s/%s',input_path_base,exp_path_2);
                        data_wrf = ncread(sprintf('%s/%s/%s',input_path_alt,timestamp_file_background,filename),data_name); % Retrieve data
                    elseif(ismember(exp_choice,alt_input_exps))
                        filename = sprintf(filename_format_free,domain,year,month,day,hour);
                        data_wrf = ncread(sprintf('%s/%s/%s',input_path_wrf,member_string,filename),data_name); % Retrieve data
                    else
                        filename = sprintf(filename_format_analysis,domain,member_idx);       
                        data_wrf = ncread(sprintf('%s/%s/%s',input_path_wrf,timestamp_file,filename),data_name); % Retrieve data
                    end
        
                    data_wrf = align_data(data_wrf,transpose_data,flip_data,0,n_idx,s_idx,e_idx,w_idx);
        
                    member_sum = member_sum + data_wrf;
                
                end
            end
            if(time_idx == 1)
                data_3z = data_wrf;
            else
                data_0z = data_wrf;
            end
        end

        data_wrf = 100*(data_0z - data_3z); % Subtract snow from before this storm (m -> cm)

        if(use_snow_inches)
            data_wrf = data_wrf*cm_to_in;
        end
            
        % Determine whether individual member will be fully processed
        process_this_member = (process_all_members || member_string == "mean");
        
        data_wrf_regridded_filename = sprintf('%s/snow_%s_%s_d0%d_%s_%s_%s_%s.mat',intermediate_path,exp_name_a,exp_name_b,domain,member_string,bao_short,upper(data_type),timestamp_file);
         
        % Don't recompute if we don't need to
        if(~remake_wrf_regrid && isfile(data_wrf_regridded_filename))
            load(data_wrf_regridded_filename,'data_wrf_regridded');
        else
            %% 2D: Loop through each obs point and find corresponding WRF value 
            data_wrf_regridded = zeros(ylen_obs,xlen_obs); % Composite anyway, no need for vertical
            
            if(show_progress)
                toc
                fprintf('Regridding...\n'); % debug statement
            end
            for y = 1:ylen_obs
                for x = 1:xlen_obs
                    [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE] = find_nn_idx_irregular(lon_obs_trim(y,x),lat_obs_trim(y,x),lon_wrf,lat_wrf,edge_buffer); % Get nn index values
                    nn_points = [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE];        
                    if(any(nn_points == 0) || any(nn_points(1:4) > xlen_wrf) || any(nn_points(5:8) > ylen_wrf)) % If indices could not be found / reference points are NaN/ goes off edge
                        data_wrf_regridded(y,x) = nodata_val; % Lowest dBZ value used on color scale
                    else
                        data_wrf_regridded(y,x) = bilinear_interpolate_irregular_to_grid(y,x,lon_obs_trim,lat_obs_trim,data_wrf,lon_wrf,lat_wrf,nn_points);
                    end
                end
            end

            % Save for exterior use
            save(data_wrf_regridded_filename,'data_wrf_regridded');
            %'%s/%s_%s_d0%d_%s_%s_%s_%s_%s.png'; % output path, exp name 1, exp name 2, domain, mem/mean, bao_short, plot_type, data type, timestamp
            
            clearvars *_data_lon *_data_lat;
        
        end

        data_wrf_compare = data_wrf_regridded;
        data_obs_compare = data_obs;

        if(make_mask)
            box_mask = zeros(size(data_wrf_regridded)); %lon_obs_trim,lat_obs_trim
            for y = 1:ylen_obs
                for x = 1:xlen_obs
                    if((lon_obs_trim(y,x) < w_box) || (lon_obs_trim(y,x) > e_box) || (lat_obs_trim(y,x) < s_box) || (lat_obs_trim(y,x) > n_box))
                        box_mask(y,x) = 0;
                    else
                        box_mas(y,x) = 1;
                    end
                end
            end
            save(sprintf('%s/%s',intermediate_path,'box_mask.mat'),'box_mask');
        end

        if(apply_mask) % Apply box mask
            %load(sprintf('%s/%s',intermediate_path,'box_mask.mat'),'box_mask');
            if(box_only)
                box_mask = inpolygon(lon_obs_trim,lat_obs_trim,bound_lon_box,bound_lat_box);
            elseif(subdomain == 0)
                box_mask = ones(size(lon_obs_trim));
            elseif(subdomain == 1)
                box_mask = inpolygon(lon_obs_trim,lat_obs_trim,bound_lon_d1,bound_lat_d1);
            elseif(subdomain == 2)
                box_mask = inpolygon(lon_obs_trim,lat_obs_trim,bound_lon_d2,bound_lat_d2);
            else
                fprintf('ERROR: Invalid subdomain.\n');
                exit(1);
            end
                
            data_wrf_masked = data_wrf_regridded;
            data_wrf_masked(data_wrf_masked == 0) = 42069;
            data_wrf_masked = data_wrf_masked.*box_mask;
            data_wrf_masked(data_wrf_masked == 0) = nodata_val;
            data_wrf_masked(data_wrf_masked == 42069) = 0;
            data_wrf_compare = data_wrf_masked;
    
            data_obs_masked = data_obs;
            data_obs_masked(data_obs_masked == 0) = 42069;
            data_obs_masked = data_obs_masked.*box_mask;
            data_obs_masked(data_obs_masked == 0) = nodata_val;
            data_obs_masked(data_obs_masked == 42069) = 0;
            data_obs_compare = data_obs_masked;
        end

        %% 2F: Make dif, compute RMSE, bias, ETS
        
        if(show_progress)
            toc
            fprintf('Computing error values...\n'); % debug statement
        end

        data_dif = data_wrf_compare - data_obs_compare;
        data_dif_clean = data_dif(~isnan(data_dif));
        num_clean_points = numel(data_dif_clean);
        num_points = numel(data_dif);
        bias = mean(data_dif_clean,'all');
        rmse_err = sqrt(sum(data_dif_clean.^2)/num_clean_points);

            
        % Brier Score
        % For each gridpoint, squared dif of (ensemble probability of >=
        % threshold) to (0/1 obs >= threshold)
        % BS = mean of the above
        % BS = (1/n)*sum((EPi - Oi)^2)
        % Need to save binary of each ensemble member... member_sum_bin?
        % Divide mem_sum_bin by num_members (k) [n = num gridpoints]
        % Then take sqdif of mem_sum_bin to obs at each grid point, take mean
        
        if(compute_error_table)

            for ets_idx = 1:2

                ets_threshold = ets_threshold_list(ets_idx);

                % Equitable Threat Score (ETS)
                data_wrf_thresh_list = find(data_wrf_compare>=ets_threshold);
                data_obs_thresh_list = find(data_obs_compare>=ets_threshold);
                data_wrf_nan_list = find(isnan(data_wrf_compare));
                data_obs_nan_list = find(isnan(data_obs_compare));
                data_wrf_bin = zeros(size(data_wrf_compare));
                data_obs_bin = zeros(size(data_obs_compare));
                data_wrf_bin(data_wrf_thresh_list) = 1;
                data_obs_bin(data_obs_thresh_list) = 1;
                data_wrf_bin(data_wrf_nan_list) = NaN;
                data_obs_bin(data_obs_nan_list) = NaN;

                if(process_this_member)

                    a = 0; % Hit (event in both obs and model)
                    b = 0; % Miss (event in obs, not in model)
                    c = 0; % False alarm (event in model, not in obs)
                    d = 0; % Correct miss (neither)

                    for y = 1:ylen_obs
                        for x = 1:xlen_obs
                            if(isnan(data_obs_compare(y,x)) || isnan(data_wrf_compare(y,x))) % Do not count NaNs
                                continue;
                            elseif(data_obs_bin(y,x) == 1) % If event occurs in obs:
                                % 1) Move out in 1D to find 9-km distance in grid
                                % points, 2) within that box, check if each point is a)
                                % within the radius, b) == 1
                                event_found = any(data_wrf_bin(max(y-grid_dist_radius,1):min(y+grid_dist_radius,ylen_obs),max(x-grid_dist_radius,1):min(x+grid_dist_radius,xlen_obs)),'all');
                                if(event_found) % Hit (a)
                                    a = a+1;
                                else % Miss (b)
                                    b = b+1;
                                end
                            elseif(data_wrf_bin(y,x) == 1) % If event occurs in model, but not right here in obs:
                                event_found = any(data_obs_bin(max(y-grid_dist_radius,1):min(y+grid_dist_radius,ylen_obs),max(x-grid_dist_radius,1):min(x+grid_dist_radius,xlen_obs)),'all');
                                if(event_found)
                                    c = c+1; % False alarm (c)
                                else % Correct miss (d)
                                    d = d+1;
                                end
                            end
                        end
                    end

                    Ar = ((a + b)*(a + c))/(a + b + c + d);
                    ETS = (a - Ar)/(a + b + c - Ar); % Equitable Threat Score (compared against hits expected of a random guess) [-1/3 : +1]
                    CSI = a/(a + b + c); % Critical Success Index (aka Threat Score) [0:1]

                else % If not computing for this member
                    ETS = 0;
                    CSI = 0;
                end

                if(ets_idx == 1)
                    ETS_a = ETS;
                    CSI_a = CSI;
                    member_sum_bin_a = member_sum_bin_a + data_wrf_bin;
                else
                    ETS_b = ETS;
                    CSI_b = CSI;
                    member_sum_bin_b = member_sum_bin_b + data_wrf_bin;
                end
            end

            % Save error scores in matrix for output at end of script
            if(member_string == "mean")
                idx = num_members+1;
                % Compute Brier Score [1 (worst) : 0 (perfect)], penalizes large
                % errors more than small ones
                ens_prob_a = member_sum_bin_a/num_members;
                ens_prob_b = member_sum_bin_b/num_members;
                brier_list_a(idx) = mean(power((ens_prob_a - data_obs_bin),2),'all',"omitnan");
                brier_list_b(idx) = mean(power((ens_prob_b - data_obs_bin),2),'all',"omitnan");
                sde_grid = sqrt(sde_grid./(num_clean_points-1));
                sde_grid_clean = sde_grid(~isnan(sde_grid));
                SDE_list(idx) = mean(sde_grid_clean,'all',"omitnan");
                sd_grid = sqrt(sd_grid./(num_clean_points-1));
                sd_grid_clean = sd_grid(~isnan(sd_grid));
                SD_list(idx) = mean(sd_grid_clean,'all',"omitnan");
            else
                idx = str2num(member_string);
            end

            RMSE_list(idx) = rmse_err;
            bias_list(idx) = bias;
            CSI_list_a(idx) = CSI_a;
            CSI_list_b(idx) = CSI_b;
            ETS_list_a(idx) = ETS_a;
            ETS_list_b(idx) = ETS_b;
        
        end
        
        clearvars sde CSI_a CSI_b ETS_a ETS_b idx ens_prob* sd_grid_clean Ar ETS CSI a b c d ets_threshold ets_idx;

        %% 2G: Plot

        if(show_progress)
            toc
            fprintf('Plotting...\n'); % debug statement
        end
        
        for plot_idx = 1:3
            if(plot_idx == 1)
                plot_type = 'fit';  % Specify fit, dif, trim, raw
                data_to_plot = data_wrf_compare;
                cmap = jet_modded;
            elseif((plot_idx == 2))
                plot_type = 'obs';
                data_to_plot = data_obs_compare;
                cmap = jet_modded;
            elseif(plot_idx == 3)
                plot_type = 'dif';
                data_to_plot = data_dif;
                cmap = 'redblue';
            else
                fprintf('ERROR: You should not ever see this message. You have a bug in the plot code.\n');
            end
            
            plot_filename = sprintf(output_format,output_path_large,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),timestamp_file);
            if(isfile(plot_filename) && ~remake_plots) % If shouldn't override existing plots, don't
                continue;
            end
            
            % Plot LARGE
            f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
            % Focus on desired area, remove whitespace
            
            if(show_progress && announce_plots)
                toc
                fprintf('%s...\n',plot_type); % debug statement
            end
            
            h = pcolor(lon_obs_trim,lat_obs_trim,data_to_plot); % Plot the data
            set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
            shading interp; % Smooth out plot from grid-boxes
            colorbar('FontSize',axes_font_size); % Make colorbar
            colormap(cmap); % Set colors
            
            if(plot_idx == 3) % Plotting dif
                caxis([clim_lower_dif clim_upper_dif]);
            else
                caxis([clim_lower clim_upper]);
            end

            if(show_progress && announce_plots)
                toc
                fprintf('Borders...\n'); % debug statement
            end
            
            % Plot state borders
            hold on;
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
            if(plot_type == "fit")
                title(sprintf(title_format_fit,exp_name_a_plot,domain,member_string,plot_type,upper(data_type),units,timestamp_title),'FontSize',title_font_size); 
            elseif(plot_type == "dif")
                title(sprintf(title_format_dif,exp_name_a_plot,exp_name_b_plot,domain,member_string,plot_type,upper(data_type),units,timestamp_title,rmse_err,bias),'FontSize',title_font_size); 
            else
                title(sprintf(title_format_trim,exp_name_b_plot,plot_type,upper(data_type),units,timestamp_title),'FontSize',title_font_size); 
            end
            xlabel('Longitude','FontSize',label_font_size);
            ylabel('Latitude','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);
            
            if(show_progress && announce_plots)
                toc
                fprintf('Saving...\n'); % debug statement
            end
            
            saveas(h,plot_filename); % Save as .png

            if(show_progress && announce_plots)
                toc
                fprintf('Making small...\n'); % debug statement
            end
            % Plot SMALL
            f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
            title(sprintf(title_format_small,upper(data_src),plot_type,upper(data_type),timestamp_file),'FontSize',title_font_size); % Replace title
            if(show_progress && announce_plots)
                toc
                fprintf('Saving small...\n'); % debug statement
            end
            saveas(h,sprintf(output_format_small,output_path_small,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),timestamp_file)); % Save as .png
            
            close('all');
            clearvars data_to_plot f h;
        end
        
        clearvars f h data_to_plot rmse_err bias;
        
        %% Clear big variables between members
        clearvars data* -except data_type data_name* data_obs_trim data_src* data_*_mean_ref;
        
    %end

    %% 2J. Write out error value table

    if(compute_error_table)
        if(show_progress)
            fprintf('Saving error table...\n'); % debug statement
            toc;
        end

        % Compute mean and median values across ensemble
        RMSE_list(num_members+2) = median(RMSE_list(1:num_members));
        RMSE_list(num_members+3) = mean(RMSE_list(1:num_members));
        RMSE_list(num_members+4) = min(RMSE_list(1:num_members));
        RMSE_list(num_members+5) = max(RMSE_list(1:num_members));
        RMSE_list(num_members+6) = std(RMSE_list(1:num_members));
        bias_list(num_members+2) = median(bias_list(1:num_members));
        bias_list(num_members+3) = mean(bias_list(1:num_members));
        bias_list(num_members+4) = min(bias_list(1:num_members));
        bias_list(num_members+5) = max(bias_list(1:num_members));
        bias_list(num_members+6) = std(bias_list(1:num_members));
        CSI_list_a(num_members+2) = median(CSI_list_a(1:num_members));
        CSI_list_a(num_members+3) = mean(CSI_list_a(1:num_members));
        CSI_list_a(num_members+4) = min(CSI_list_a(1:num_members));
        CSI_list_a(num_members+5) = max(CSI_list_a(1:num_members));
        CSI_list_a(num_members+6) = std(CSI_list_a(1:num_members));
        CSI_list_b(num_members+2) = median(CSI_list_b(1:num_members));
        CSI_list_b(num_members+3) = mean(CSI_list_b(1:num_members));
        CSI_list_b(num_members+4) = min(CSI_list_b(1:num_members));
        CSI_list_b(num_members+5) = max(CSI_list_b(1:num_members));
        CSI_list_b(num_members+6) = std(CSI_list_b(1:num_members));
        ETS_list_a(num_members+2) = median(ETS_list_a(1:num_members));
        ETS_list_a(num_members+3) = mean(ETS_list_a(1:num_members));
        ETS_list_a(num_members+4) = min(ETS_list_a(1:num_members));
        ETS_list_a(num_members+5) = max(ETS_list_a(1:num_members));
        ETS_list_a(num_members+6) = std(ETS_list_a(1:num_members));
        ETS_list_b(num_members+2) = median(ETS_list_b(1:num_members));
        ETS_list_b(num_members+3) = mean(ETS_list_b(1:num_members));
        ETS_list_b(num_members+4) = min(ETS_list_b(1:num_members));
        ETS_list_b(num_members+5) = max(ETS_list_b(1:num_members));
        ETS_list_b(num_members+6) = std(ETS_list_b(1:num_members));

        % Make the numbers pretty
        RMSE_list = round(RMSE_list,error_dec);
        bias_list = round(bias_list,error_dec);
        SDE_list = round(SDE_list,error_dec);
        SD_list = round(SD_list,error_dec);
        CSI_list_a = round(CSI_list_a,error_dec);
        CSI_list_b = round(CSI_list_b,error_dec);
        ETS_list_a = round(ETS_list_a,error_dec);
        ETS_list_b = round(ETS_list_b,error_dec);
        brier_list_a = round(brier_list_a,error_dec);
        brier_list_b = round(brier_list_b,error_dec);

        % Avoid hardcoding the threshold numbers
        CSI_colnames = [sprintf("CSI%d",ets_threshold_list(1)),sprintf("CSI%d",ets_threshold_list(2))];
        ETS_colnames = [sprintf("ETS%d",ets_threshold_list(1)),sprintf("ETS%d",ets_threshold_list(2))];
        Brier_colnames = [sprintf("BS%d",ets_threshold_list(1)),sprintf("BS%d",ets_threshold_list(2))];

        % Write
        error_q_table = table(RMSE_list',bias_list',ETS_list_a',ETS_list_b',CSI_list_a',CSI_list_b',SD_list',SDE_list',brier_list_a',brier_list_b','RowNames',err_rowNames,'VariableNames',["RMSE","Bias",ETS_colnames(1),ETS_colnames(2),CSI_colnames(1),CSI_colnames(2),"SD","SDE",Brier_colnames(1),Brier_colnames(2)]);
        writetable(error_q_table,sprintf('%s/ev_table_snow_%s_%s.txt',output_path_large,run_name,timestamp_file),'Delimiter','\t','WriteRowNames',true);
    end
    
    clearvars error_q_table RMSE_list bias_list SDE_list CSI_list* ETS_list* brier_list* RPS_list* member_sum* *_colnames data_wrf_mean_ref sd_grid sde_grid;
        



% Get total runtime and print to stdout
runtime = toc;
hours = floor(runtime/3600);
mins = floor((runtime/60) - (hours*60));
secs = toc - (hours*3600) - (mins*60);
fprintf('Done. Total script runtime = %02.f:%02.f:%02.f\n',hours,mins,secs);
fprintf('#---------------------------------------------------------------------#\n');

% END

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


                