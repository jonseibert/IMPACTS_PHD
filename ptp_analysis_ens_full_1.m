% ptp_analysis_ens_full.m (v2.2.2)
%
% Purpose: Generate point-to-point (PTP) diagnostic comparison plots of 
% outputs from IMPACTS WRF to NEXRAD observations and each other. This 
% version computes outputs for all ensemble members, as well as their mean,
% and allows for comparison of multiple forecast variables.
% 
% Author(s): Jon Seibert
% Last updated: 5 Dec 2024
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

run_name = 'extra_plots'; % Name of output folder

% Which experiments to use (0 is obs, only for B)
% To compare background and analysis, use the same experiment # for both
% To compare to observations, use 0 for experiment B
exp_choice_a = 1; 
exp_choice_b = 1;

data_type = "T"; % Working name of variable being analyzed (REFL, T, MSLP, GPH, VORT, WIND, Q, OMEGA) [+T2M,SNOWH(depth),SNOWNC(grid total per timestep)]
use_base_refl = 1; % Whether to override REFL inputs with precalculated Base Reflectivity values

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

% Vertical levels to use if is 3D
is_3d = 1; % Is the input data 3-dimensional?
pressure_target = 800; %mb [0 = vc, 1 = no height, 1000 = sfc, other = x]

z_composite_data = 0; % Unused if use_base_ref is true
multi_height = 0; % Whether to loop over multiple target pressure values
pressure_target_list = [1000,925]; % [925,850,700]; Overrides above pressure_target value if multi_height true

% Plot suite selections
smooth_data = 0; % If true AND the variable isn't already being smoothed, smooth it before plotting (gaussian)
remake_wrf_regrid = 0; % Whether to override existing saved WRF regrids
compute_error_table = 0; % Whether to compute and output error suite values

remake_plots = 1; % If true, overrides existing outputs. If false, only plots 'new' files.
plot_mean_only = 0; % Suppresses individual member plots
plot_first_and_mean_only = 0;
plot_raw = 0; % Base values
plot_dif = 1; % Difference plots
plot_trim = 0; % Base values, but for obs

% Progress messages
show_plots = 0; % 0: Suppress plot popups [REQUIRES MATLAB RESTART TO UNDO]
show_progress = 1; % If true, prints out progress bar/debug messages
announce_plots = 0; % If true, prints out individual plot messages

%% 0B. General Settings
% Should not be changed without good reason

local_wd = "C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code";
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
    input_path_base = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020';
    input_path_gis = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    intermediate_path = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/output/testing';
    path_to_code = "C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code";
    path_to_extra_code = './downloaded_code';
    input_path_ecmwf = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Data/ECMWF';
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
        clim_lower_dif = -1.2e1;
        clim_upper_dif = 1.2e1;
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
if(compare_to_obs)

    input_path_wrf = input_path_a;    
    exp_choice = exp_choice_a;
    exp_name = exp_name_a;
    input_path_base_refl = sprintf('%s/BASE_REF/%s',input_path_base,exp_name);
    
    if(show_progress)
        toc
        fprintf('Comparison: %s-obs, %s\n',exp_name_a,upper(data_type));
    end
    
    if(upper(data_type) ~= "REFL")
        error(sprintf("Data Type: %s. Only reflectivity observations are supported.",upper(data_type)));
    end
    
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
    clearvars *gis*raw x y;
    
    % Define empty grid to fit wrf data into
    dimensions_gis_trim = size(lon_gis_trim);
    ylen_gis = dimensions_gis_trim(1);
    xlen_gis = dimensions_gis_trim(2);
    
    % Determine short grid box to search in for neighborhood version of ETS/CSI
    x_center = round(xlen_gis/2);
    y_center = round(ylen_gis/2);
    grid_dist_radius = 0;
    for x = x_center:xlen_gis
        d = haversine_distance(lon_gis_trim(x,y_center),lat_gis_trim(x,y_center),lon_gis_trim(x_center,y_center),lat_gis_trim(x_center,y_center));
        if(d > ets_radius)
            grid_dist_radius = abs(x_center - x);
            break;
        end
    end
    clearvars x_center y_center;
    
    %% 2B. For each timestamp being analyzed:
    for time_idx = 1:time_count
       
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
        sde_grid = zeros(ylen_gis,xlen_gis);
        
        % Read in GIS data
        % Convert from greyscale to dBZ values: dBZ = grey*0.5 - 32.5; (Still a bit unclear on WHY this is the conversion.)
        data_gis = double(imread(sprintf(file_format_gis,input_path_gis,year,month,day,hour,minute)))*0.5 - 32.5; 
        data_gis_trim = data_gis(s_idx:1:n_idx,w_idx:1:e_idx); % Trim down to WRF domain
        clearvars data_gis;
        
        % Regrid WRF data onto trimmed GIS grid
        % Brute force bilinear interp to make the WRF data match the shape of the GIS data
        % NOTE: Assumes that the WRF data fully encompasses the trimmed GIS region, so that
        % edge cases can work as normal!
    
        member_sum = zeros(ylen_wrf,xlen_wrf);
        member_sum_bin_a = zeros(ylen_gis,xlen_gis);
        member_sum_bin_b = zeros(ylen_gis,xlen_gis); 
        
        % Load in mean for SD reference
        load(sprintf(input_format_base_refl,input_path_base_refl,data_src,'base',domain,'mean',upper(data_type),timestamp_file),'wrf_base_refl');
        data_a_mean_ref = wrf_base_refl;
    
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
    
                if(use_base_refl)
                    load(sprintf(input_format_base_refl,input_path_base_refl,data_src,'base',domain,member_string,upper(data_type),timestamp_file),'wrf_base_refl');
                    data_wrf = wrf_base_refl;
                    clearvars wrf_base_refl;
                elseif(ismember(exp_choice,alt_input_exps))
                    filename = sprintf(filename_format_free,domain,year,month,day,hour);
                    data_wrf = ncread(sprintf('%s/%s/%s',input_path_wrf,member_string,filename),data_name); % Retrieve data
                else
                    filename = sprintf(filename_format_analysis,domain,member_idx);       
                    data_wrf = ncread(sprintf('%s/%s/%s',input_path_wrf,timestamp_file,filename),data_name); % Retrieve data
                end
    
                if(~use_base_refl)
                    % If data needs to be composited down to 2D, do so
                    if(z_composite_data)
                        data_wrf = max(data_wrf,[],3);
                    end
    
                    % If data needs to be transposed to align with lat/lon grid, do so
                    if(transpose_data) 
                        data_wrf = data_wrf';
                    end
    
                    if(flip_data)
                        data_wrf = flip(data_wrf);
                    end
                end
    
                member_sum = member_sum + data_wrf;
            
            end
            
            % Determine whether individual member will be fully processed
            process_this_member = (process_all_members || member_string == "mean");
            
            data_wrf_regridded_filename = sprintf('%s/data_%s_%s_d0%d_%s_%s_%s_%s.mat',intermediate_path,exp_name_a,exp_name_b,domain,member_string,bao_short,upper(data_type),timestamp_file);
             
            % Don't recompute if we don't need to
            if(~remake_wrf_regrid && isfile(data_wrf_regridded_filename))
                load(data_wrf_regridded_filename,'data_wrf_regridded');
            else
                %% 2D: Loop through each GIS point and find corresponding WRF value 
                data_wrf_regridded = zeros(ylen_gis,xlen_gis); % Composite anyway, no need for vertical
                
                if(show_progress)
                    toc
                    fprintf('Regridding...\n'); % debug statement
                end
                for y = 1:ylen_gis
                    for x = 1:xlen_gis
                        [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE] = find_nn_idx_irregular(lon_gis_trim(y,x),lat_gis_trim(y,x),lon_wrf,lat_wrf,edge_buffer); % Get nn index values
                        nn_points = [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE];        
                        if(any(nn_points == 0) || any(nn_points(1:4) > xlen_wrf) || any(nn_points(5:8) > ylen_wrf)) % If indices could not be found / reference points are NaN/ goes off edge
                            data_wrf_regridded(y,x) = nodata_val; % Lowest dBZ value used on color scale
                        else
                            data_wrf_regridded(y,x) = bilinear_interpolate_irregular_to_grid(y,x,lon_gis_trim,lat_gis_trim,data_wrf,lon_wrf,lat_wrf,nn_points);
                        end
                    end
                end
    
                % Calculate the bounds of the real data for proper comparison with
                % GIS / count how many NaN rows and columns there are
                [first_data_lat,last_data_lat,first_data_lon,last_data_lon] = find_data_edges_nas(ylen_gis,xlen_gis,data_gis_trim);
    
                % Fill in NaNs on the border where they are in EXRAD
                if(first_data_lat > 1)
                    data_wrf_regridded(1:first_data_lat,:) = NaN;
                end
                if(last_data_lat < ylen_gis)
                    data_wrf_regridded(last_data_lat:ylen_gis,:) = NaN;
                end
                if(first_data_lon > 1)
                    data_wrf_regridded(:,1:first_data_lon) = NaN;
                end
                if(last_data_lon < xlen_gis)
                    data_wrf_regridded(:,last_data_lon:xlen_gis) = NaN;
                end
    
                % Save for exterior use
                save(data_wrf_regridded_filename,'data_wrf_regridded');
                %'%s/%s_%s_d0%d_%s_%s_%s_%s_%s.png'; % output path, exp name 1, exp name 2, domain, mem/mean, bao_short, plot_type, data type, timestamp
                
                clearvars *_data_lon *_data_lat;
            
            end
            
            %% 2E: Apply mask
            
            if(show_progress)
                toc
                fprintf('Masking...\n'); % debug statement
            end
    
            load('mask_wrf_on_gis','wrf_gis_mask');
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
                
            data_a_masked = data_wrf_regridded;
            data_a_masked(data_a_masked == 0) = 42069;
            data_a_masked = data_a_masked.*wrf_gis_mask;
            data_a_masked = data_a_masked.*gis_sbd_mask;
            data_a_masked(data_a_masked == 0) = nodata_val_gis;
            data_a_masked(data_a_masked == 42069) = 0;
            data_wrf_compare = data_a_masked;
    
            data_gis_masked = data_gis_trim;
            data_gis_masked(data_gis_masked == 0) = 42069;
            data_gis_masked = data_gis_masked.*wrf_gis_mask;
            data_gis_masked = data_gis_masked.*gis_sbd_mask;
            data_gis_masked(data_gis_masked == 0) = nodata_val_gis;
            data_gis_masked(data_gis_masked == 42069) = 0;
            data_gis_compare = data_gis_masked;
            
            if(member_idx == 1)
                data_a_mean_ref(data_a_mean_ref == 0) = 42069;
                data_a_mean_ref = data_a_mean_ref.*wrf_wrf_mask;
                data_a_mean_ref(data_a_mean_ref == 0) = nodata_val_gis;
                data_a_mean_ref(data_a_mean_ref == 42069) = 0;
            end  
            
            data_wrf_raw_masked = data_wrf;
            data_wrf_raw_masked(data_wrf_raw_masked == 0) = 42069;
            data_wrf_raw_masked = data_wrf_raw_masked.*wrf_wrf_mask;
            data_wrf_raw_masked = data_wrf_raw_masked.*wrf_sbd_mask;
            data_wrf_raw_masked(data_wrf_raw_masked == 0) = nodata_val_gis;
            data_wrf_raw_masked(data_wrf_raw_masked == 42069) = 0;
           
            clearvars wrf_gis_mask wrf_wrf_mask;
    
            %% 2F: Make dif, compute RMSE, bias, ETS
            
            if(show_progress)
                toc
                fprintf('Computing error values...\n'); % debug statement
            end
    
            data_dif = data_wrf_compare - data_gis_compare;
            data_variance = data_wrf_raw_masked - data_a_mean_ref;
            data_dif_clean = data_dif(~isnan(data_dif));
            num_clean_points = numel(data_dif_clean);
            num_points = numel(data_dif);
            bias = mean(data_dif_clean,'all');
            rmse_err = sqrt(sum(data_dif_clean.^2)/num_clean_points);
            sde_grid = sde_grid + (data_dif.^2);
            sd_grid = sd_grid + (data_variance.^2);
                
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
                    data_gis_thresh_list = find(data_gis_compare>=ets_threshold);
                    data_wrf_nan_list = find(isnan(data_wrf_compare));
                    data_gis_nan_list = find(isnan(data_gis_compare));
                    data_wrf_bin = zeros(size(data_wrf_compare));
                    data_gis_bin = zeros(size(data_gis_compare));
                    data_wrf_bin(data_wrf_thresh_list) = 1;
                    data_gis_bin(data_gis_thresh_list) = 1;
                    data_wrf_bin(data_wrf_nan_list) = NaN;
                    data_gis_bin(data_gis_nan_list) = NaN;
    
                    if(process_this_member)
    
                        a = 0; % Hit (event in both obs and model)
                        b = 0; % Miss (event in obs, not in model)
                        c = 0; % False alarm (event in model, not in obs)
                        d = 0; % Correct miss (neither)
    
                        for y = 1:ylen_gis
                            for x = 1:xlen_gis
                                if(isnan(data_gis_compare(y,x)) || isnan(data_wrf_compare(y,x))) % Do not count NaNs
                                    continue;
                                elseif(data_gis_bin(y,x) == 1) % If event occurs in obs:
                                    % 1) Move out in 1D to find 9-km distance in grid
                                    % points, 2) within that box, check if each point is a)
                                    % within the radius, b) == 1
                                    event_found = any(data_wrf_bin(max(y-grid_dist_radius,1):min(y+grid_dist_radius,ylen_gis),max(x-grid_dist_radius,1):min(x+grid_dist_radius,xlen_gis)),'all');
                                    if(event_found) % Hit (a)
                                        a = a+1;
                                    else % Miss (b)
                                        b = b+1;
                                    end
                                elseif(data_wrf_bin(y,x) == 1) % If event occurs in model, but not right here in obs:
                                    event_found = any(data_gis_bin(max(y-grid_dist_radius,1):min(y+grid_dist_radius,ylen_gis),max(x-grid_dist_radius,1):min(x+grid_dist_radius,xlen_gis)),'all');
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
                    brier_list_a(idx) = mean(power((ens_prob_a - data_gis_bin),2),'all',"omitnan");
                    brier_list_b(idx) = mean(power((ens_prob_b - data_gis_bin),2),'all',"omitnan");
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
    
            if(plot_mean_only && member_idx ~= num_members+1)
                if(show_progress)
                    toc
                    fprintf('Skipping member plots...\n'); % debug statement
                end
                continue;
            elseif(plot_first_and_mean_only && member_idx ~= 1 && member_idx ~= num_members+1)
                if(show_progress)
                    toc
                    fprintf('Skipping all but 1st member plots...\n'); % debug statement
                end
                continue;
            end
    
            if(show_progress)
                toc
                fprintf('Plotting...\n'); % debug statement
            end
            
            for plot_idx = 1:3
                if(plot_idx == 1)
                    if(~plot_raw) continue; end
                    plot_type = 'fit';  % Specify fit, dif, trim, raw
                    data_to_plot = data_wrf_compare;
                    cmap = 'reflmap';
                elseif((plot_idx == 2))
                    if(~plot_trim || member_idx ~= 1) continue; end
                    plot_type = 'trim';
                    data_to_plot = data_gis_compare;
                    cmap = 'reflmap';
                elseif(plot_idx == 3)
                    if(~plot_dif) continue; end
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
                
                h = pcolor(lon_gis_trim,lat_gis_trim,data_to_plot); % Plot the data
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
                    title(sprintf(title_format_fit,exp_name_a_plot,exp_name_b_plot,domain,member_string,plot_type,upper(data_type),units,timestamp_title),'FontSize',title_font_size); 
                elseif(plot_type == "dif")
                    title(sprintf(title_format_dif,exp_name_a_plot,exp_name_b_plot,domain,member_string,plot_type,upper(data_type),units,timestamp_title,rmse_err,bias),'FontSize',title_font_size); 
                else
                    title(sprintf(title_format_trim,exp_name_b_plot,plot_type,upper(data_type),units,timestamp_title),'FontSize',title_font_size); 
                end
                xlabel('Longitude (deg)','FontSize',label_font_size);
                ylabel('Latitude (deg)','FontSize',label_font_size);
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
    
            %% 2H: Plot WRF on WRF grid for reference
            
            plot_type = 'raw';
            plot_filename = sprintf(output_format,output_path_large,'wrf',lower(exp_name_a_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),timestamp_file);
            
            if(plot_original_wrf && (remake_plots || ~isfile(plot_filename)))     
                if(show_progress && announce_plots)
                    toc
                    fprintf('Plotting Original WRF...\n'); % debug statement
                end
       
                cmap = 'reflmap'; 
    
                % Plot LARGE
                f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                % Focus on desired area, remove whitespace
                
                h = pcolor(lon_wrf,lat_wrf,data_wrf); % Plot the data
                set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                shading interp; % Smooth out plot from grid-boxes
                colorbar('FontSize',axes_font_size); % Make colorbar
                colormap(cmap); % Set colors
                caxis([clim_lower clim_upper]);
    
                % Plot state borders
                hold on;
                borders('continental us','black','linewidth',1); 
                hold off;
                
                if(limit_borders)
                    xlim([w_lim e_lim]);
                    ylim([s_lim n_lim]);
                end
    
                % Apply labels
                title(sprintf(title_format_raw,exp_name_a_plot,domain,member_string,upper(data_type),units,timestamp_title),'FontSize',title_font_size); 
                xlabel('Longitude (deg)','FontSize',label_font_size);
                ylabel('Latitude (deg)','FontSize',label_font_size);
                set(gca,'Fontsize',axes_font_size);
                saveas(h,plot_filename); % Save as .png
                close('all');
            end 
            
            %% 2I. Plot Standard Deviation across Ensemble
            
            plot_type = 'sd';
            plot_filename = sprintf(output_format,output_path_large,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),timestamp_file);
            
            if(plot_sd && member_string == "mean" && (remake_plots || ~isfile(plot_filename)))
                if(~compute_error_table)
                    fprintf('Error: Cannot compute STDEV independently of err table. Skipping.\n');
                else
                    if(show_progress && announce_plots)
                        toc
                        fprintf('Plotting StDev...\n'); % debug statement
                    end
                    data_to_plot = sd_grid;
                    cmap = 'jet';
    
                    % Plot LARGE
                    f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                    h = pcolor(lon_wrf,lat_wrf,data_to_plot); % Plot the data
                    set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                    shading interp; % Smooth out plot from grid-boxes
                    colorbar('FontSize',axes_font_size); % Make colorbar
                    colormap(cmap); % Set colors
                    caxis([0 sd_max]);
    
                    % Plot state borders
                    hold on;
                    borders('continental us','color','#bcbcbc','linewidth',1); 
                    hold off;
    
                    if(limit_borders)
                        xlim([w_lim e_lim]);
                        ylim([s_lim n_lim]);
                    end
    
                    % Apply labels
                    title(sprintf(title_format_sd,exp_name_a_plot,domain,'SD',upper(data_type),units,timestamp_title),'FontSize',title_font_size); 
                    xlabel('Longitude (deg)','FontSize',label_font_size);
                    ylabel('Latitude (deg)','FontSize',label_font_size);
                    set(gca,'Fontsize',axes_font_size);
                    saveas(h,plot_filename); % Save as .png
    
                    % Plot SMALL
                    f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
                    title(sprintf(title_format_small,upper(data_src),upper(plot_type),upper(data_type),timestamp_file),'FontSize',title_font_size); % Replace title
                    saveas(h,sprintf(output_format_small,output_path_small,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),timestamp_file)); % Save as .png
                    close('all');
                    
                    %%% Plot SDE as well %%%
                    plot_type = 'sde';
                    plot_filename = sprintf(output_format,output_path_large,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),timestamp_file);
                    data_to_plot = sde_grid;
                    cmap = 'jet';
    
                    % Plot LARGE
                    f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                    h = pcolor(lon_gis_trim,lat_gis_trim,data_to_plot); % Plot the data
                    set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                    shading interp; % Smooth out plot from grid-boxes
                    colorbar('FontSize',axes_font_size); % Make colorbar
                    colormap(cmap); % Set colors
                    %caxis([0 sd_max]);
    
                    % Plot state borders
                    hold on;
                    borders('continental us','color','#bcbcbc','linewidth',1); 
                    hold off;
    
                    if(limit_borders)
                        xlim([w_lim e_lim]);
                        ylim([s_lim n_lim]);
                    end
    
                    % Apply labels
                    title(sprintf(title_format_sd,exp_name_a_plot,domain,'SDE',upper(data_type),units,timestamp_title),'FontSize',title_font_size); 
                    xlabel('Longitude (deg)','FontSize',label_font_size);
                    ylabel('Latitude (deg)','FontSize',label_font_size);
                    set(gca,'Fontsize',axes_font_size);
                    saveas(h,plot_filename); % Save as .png
    
                    % Plot SMALL
                    f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
                    title(sprintf(title_format_small,upper(data_src),upper(plot_type),upper(data_type),timestamp_file),'FontSize',title_font_size); % Replace title
                    saveas(h,sprintf(output_format_small,output_path_small,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),timestamp_file)); % Save as .png
                    close('all');
                end
            end
            
            clearvars f h data_to_plot rmse_err bias;
            
            %% Clear big variables between members
            clearvars data* -except data_type data_name* data_gis_trim data_src* data_*_mean_ref;
            
        end
    
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
            writetable(error_q_table,sprintf('%s/ev_table_%s_%s.txt',output_path_large,run_name,timestamp_file),'Delimiter','\t','WriteRowNames',true);
        end
        
        clearvars error_q_table RMSE_list bias_list SDE_list CSI_list* ETS_list* brier_list* RPS_list* member_sum* *_colnames data_wrf_mean_ref sd_grid sde_grid;
        
        %% 2K. Time increment
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
else
    %% 3A: MAIN LOOP: If comparing two WRF runs together
        
    % Identify filepaths
    %input_path_a = sprintf('%s/%s',input_path_base,exp_name_a);
    %input_path_b = sprintf('%s/%s',input_path_base,exp_name_b);
    input_path_base_refl_a = sprintf('%s/BASE_REF/%s',input_path_base,exp_name_a);
    input_path_base_refl_b = sprintf('%s/BASE_REF/%s',input_path_base,exp_name_b);
    
    load(sprintf('%s/%s',intermediate_path,'latlon_gis_trim.mat'),'lon_gis_trim','lat_gis_trim'); % Load in preestablished latlon values
    gauss_filt = fspecial('gaussian',gauss_width,gauss_sd);
    
    if(show_progress)
        toc
        fprintf('Comparison: %s-%s, %s\n',exp_name_a,exp_name_b,upper(data_type));
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
        RMSE_list = zeros(size(err_rowNames));
        bias_list = zeros(size(err_rowNames));
        ens_mem_a = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_b = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_u_a = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_u_b = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_v_a = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_v_b = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_vort_a = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_vort_b = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_pva_a = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_pva_b = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_gph_a = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_gph_b = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_omega_a = zeros(ylen_wrf,xlen_wrf,num_members);
        ens_mem_omega_b = zeros(ylen_wrf,xlen_wrf,num_members);
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
    
        %obs_hour = mod(obs_time,24);
         %obs_days = mod(floor(obs_time/24),365.24);
        %obs_day = ones(size(obs_time))*start_day;
        %obs_month = ones(size(obs_time))*start_month;
        %obs_year = 1900 + floor(obs_time/(24*365.24));
        
        for height_idx = 1:num_heights
            if(multi_height)
                pressure_target = pressure_target_list(height_idx);
            end
            
            % 500: 3.5e-3
            % 700: 0.8e-2
            % 800: 1.0e-2
            % 850: 1.0e-2
            % 925: 1.5e-2 
            if(upper(data_type) == "Q")
                switch pressure_target
                    case 500
                        clim_upper = 3.5e-3;
                    case 700
                        clim_upper = 1.0e-2;
                    case 800
                        clim_upper = 1.0e-2;
                    case 850
                        clim_upper = 1.0e-2;
                    case 925
                        clim_upper = 1.0e-2;
                    otherwise
                        clim_upper = 3.5e-3;
                end
            end
    
        % Read in each member and compute the mean
        %% 3C: Member loop
        for member_idx = 1:(num_members+1)
            
            % If all members have been computed, take mean
            if(member_idx == (num_members+1)) 
                member_string = 'mean';
                data_a = mean(ens_mem_a,3,"omitnan");
                data_b = mean(ens_mem_b,3,"omitnan");
                data_u_a = mean(ens_mem_u_a,3,"omitnan");
                data_v_a = mean(ens_mem_v_a,3,"omitnan");
                data_u_b = mean(ens_mem_u_b,3,"omitnan");
                data_v_b = mean(ens_mem_v_b,3,"omitnan");
                data_vort_a = mean(ens_mem_vort_a,3,"omitnan");
                data_vort_b = mean(ens_mem_vort_b,3,"omitnan");
                data_pva_a = mean(ens_mem_pva_a,3,"omitnan");
                data_pva_b = mean(ens_mem_pva_b,3,"omitnan");
                data_gph_a = mean(ens_mem_gph_a,3,"omitnan");
                data_gph_b = mean(ens_mem_gph_b,3,"omitnan");
                data_omega_a = mean(ens_mem_omega_a,3,"omitnan");
                data_omega_b = mean(ens_mem_omega_b,3,"omitnan");
                L_loc(num_members+1,:,1) = mean(L_loc(1:num_members,:,1));
                L_loc(num_members+1,:,2) = mean(L_loc(1:num_members,:,2));
                L_val(num_members+1,1) = mean(L_val(1:num_members,1));
                L_val(num_members+1,2) = mean(L_val(1:num_members,2));
                L_dist(num_members+1,1) = mean(L_dist(1:num_members,1));
                L_dist(num_members+1,2) = mean(L_dist(1:num_members,2));
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

                flattened = false;
                aligned = false;
                
                %if(upper(data_type) ~= "REFL") % If working with non-reflectivity fields
                
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
                data_t2m_a = ncread(full_filepath_a,data_name_t2m); % 2m Temperature (K)
                data_q_a = ncread(full_filepath_a,data_name_q); % Water vapor (kg/kg)
                data_u_raw_a = ncread(full_filepath_a,data_name_u); % Wind (m/s)
                data_v_raw_a = ncread(full_filepath_a,data_name_v); % Wind (m/s)
                data_w_raw_a = ncread(full_filepath_a,data_name_w); % Wind (m/s)

                data_pb_b = ncread(full_filepath_b,data_name_pres);
                data_pp_b = ncread(full_filepath_b,data_name_pres_prime); 
                data_ptp_b = ncread(full_filepath_b,data_name_temp); 
                data_p_b = data_pb_b + data_pp_b; 
                data_pt_b = data_ptp_b + 300; 
                data_t_b = data_pt_b.*((data_p_b./100000).^(0.287));
                data_gpb_raw_b = ncread(full_filepath_b,data_name_gp); 
                data_gpp_raw_b = ncread(full_filepath_b,data_name_gp_prime); 
                data_p_b = data_p_b./100; 
                data_t2m_b = ncread(full_filepath_b,data_name_t2m);
                data_q_b = ncread(full_filepath_b,data_name_q); 
                data_u_raw_b = ncread(full_filepath_b,data_name_u);
                data_v_raw_b = ncread(full_filepath_b,data_name_v);
                data_w_raw_b = ncread(full_filepath_b,data_name_w);
                
                if(member_idx == 1)
                    terrain_height = ncread(full_filepath_a,data_name_terrain_height); % Height of physical terrain (m)
                end

                % Define holding containers
                data_u_a = zeros(size(data_t_a));
                data_v_a = zeros(size(data_t_a));
                data_u_b = zeros(size(data_t_b));
                data_v_b = zeros(size(data_t_b));
                data_w_a = zeros(size(data_t_b));
                data_w_b = zeros(size(data_t_b));
                data_gpb_a = zeros(size(data_t_a));
                data_gpb_b = zeros(size(data_t_b));
                data_gpp_a = zeros(size(data_t_a));
                data_gpp_b = zeros(size(data_t_b));
                
                
                % Fix staggered WRF grids for wind, geopotential
                for y = 1:ylen_wrf_raw
                    for x = 1:xlen_wrf_raw % NOT ALIGNED TO (Y,X) YET
                        data_u_a(x,y,:) = (data_u_raw_a(x,y,:) + data_u_raw_a(x+1,y,:))./2;
                        data_u_b(x,y,:) = (data_u_raw_b(x,y,:) + data_u_raw_b(x+1,y,:))./2;
                        data_v_a(x,y,:) = (data_v_raw_a(x,y,:) + data_v_raw_a(x,y+1,:))./2;
                        data_v_b(x,y,:) = (data_v_raw_b(x,y,:) + data_v_raw_b(x,y+1,:))./2;
                    end
                end
                
                for z = 1:zlen_wrf
                    data_gpb_a(:,:,z) = (data_gpb_raw_a(:,:,z) + data_gpb_raw_a(:,:,z+1))./2;
                    data_gpp_a(:,:,z) = (data_gpp_raw_a(:,:,z) + data_gpp_raw_a(:,:,z+1))./2;
                    data_gpb_b(:,:,z) = (data_gpb_raw_b(:,:,z) + data_gpb_raw_b(:,:,z+1))./2;
                    data_gpp_b(:,:,z) = (data_gpp_raw_b(:,:,z) + data_gpp_raw_b(:,:,z+1))./2;
                    data_w_a(:,:,z) = (data_w_raw_a(:,:,z) + data_w_raw_a(:,:,z+1))./2;
                    data_w_b(:,:,z) = (data_w_raw_b(:,:,z) + data_w_raw_b(:,:,z+1))./2;
                end
                
                % Compute geopotential height
                data_gp_a = data_gpb_a + data_gpp_a; %m^2/s^2
                data_gph_a = data_gp_a./g; %m
                data_gp_b = data_gpb_b + data_gpp_b;
                data_gph_b = data_gp_b./g;
                data_gph_a = data_gph_a/10; % Convert from m to dam (decameters)
                data_gph_b = data_gph_b/10;
                
                % Compute MSLP: P1 = P2*e^(((g/(R*Tv))*(z2-z1))
                z1 = 0; % m
                h1 = 1;
                h2 = 10;
                
                P2 = data_p_a(:,:,h2); % Lowest level of P
                Tv_1 = data_t_a(:,:,h1).*(1 + 0.608*data_q_a(:,:,h1));
                Tv_2 = data_t_a(:,:,h2).*(1 + 0.608*data_q_a(:,:,h2));
                Tv = (Tv_1 + Tv_2)/2;
                z2 = (data_gp_a(:,:,h2).*Re)./((g*Re) - data_gp_a(:,:,h2)); % Height above mean sea level
                mslp_a = P2.*exp((g./(R.*Tv)).*(z2-z1)); % MSLP (mb)
                
                % Compute statistics relative to observed sfc Low
                [p_min,p_i] = min(mslp_a,[],'all','linear');
                [p_min_row,p_min_col] = ind2sub(size(mslp_a),p_i);
                %p_min_row = mod(p_i,size(mslp_a,1));
                %p_min_col = 1 + floor(p_i/size(mslp_a,1));
                L_val(member_idx,1) = p_min;
                L_loc(member_idx,1,1) = lon_wrf_raw(p_min_row,p_min_col);
                L_loc(member_idx,2,1) = lat_wrf_raw(p_min_row,p_min_col);
                L_dist(member_idx,1) = haversine_distance(lon_wrf_raw(p_min_row,p_min_col),lat_wrf_raw(p_min_row,p_min_col),L_lon_o,L_lat_o);
                
                P2 = data_p_b(:,:,h2); % Lowest level of P
                Tv_1 = data_t_b(:,:,h1).*(1 + 0.608*data_q_b(:,:,h1));
                Tv_2 = data_t_b(:,:,h2).*(1 + 0.608*data_q_b(:,:,h2));
                Tv = (Tv_1 + Tv_2)/2;
                z2 = (data_gp_b(:,:,h2).*Re)./((g*Re) - data_gp_b(:,:,h2)); % Height above mean sea level
                mslp_b = P2.*exp((g./(R.*Tv)).*(z2-z1)); % MSLP
                
                [p_min,p_i] = min(mslp_b,[],'all','linear');
                [p_min_row,p_min_col] = ind2sub(size(mslp_b),p_i);
                %p_min_row = mod(p_i,size(mslp_b,1));
                %p_min_col = 1 + floor(p_i/size(mslp_b,1));
                L_val(member_idx,2) = p_min;
                L_loc(member_idx,1,2) = lon_wrf_raw(p_min_row,p_min_col);
                L_loc(member_idx,2,2) = lat_wrf_raw(p_min_row,p_min_col);
                L_dist(member_idx,2) = haversine_distance(lon_wrf_raw(p_min_row,p_min_col),lat_wrf_raw(p_min_row,p_min_col),L_lon_o,L_lat_o);
                
                model_height_a = ((data_gpb_a+data_gpp_a)./g) - terrain_height; % model height above the *ground* (m)
                model_height_b = ((data_gpb_b+data_gpp_b)./g) - terrain_height; % model height above the *ground* (m)
                
                data_t_a = data_t_a - 273.15; % Convert T from K to C
                data_t_b = data_t_b - 273.15;

                data_omega_a = zeros(size(data_t_a));
                data_omega_b = zeros(size(data_t_b));

                % Also calculate omega (uplift, dp/dt, vertical wind in pressure coords)
                for z = 1:zlen_wrf
                    if(z == 1)
                        dp = (data_p_a(:,:,z+1) - data_p_a(:,:,z));
                        dz = (model_height_a(:,:,z+1) - model_height_a(:,:,z));
                    elseif(z == zlen_wrf)
                        dp = (data_p_a(:,:,z) - data_p_a(:,:,z-1));
                        dz = (model_height_a(:,:,z) - model_height_a(:,:,z-1));
                    else
                        dp = (data_p_a(:,:,z+1) - data_p_a(:,:,z-1));
                        dz = (model_height_a(:,:,z+1) - model_height_a(:,:,z-1));
                    end
                    data_omega_a(:,:,z) = data_w_a(:,:,z).*(dp./dz); %dp/dt = (dz/dt)*(dp/dz)
                    if(z == 1)
                        dp = (data_p_b(:,:,z+1) - data_p_b(:,:,z));
                        dz = (model_height_b(:,:,z+1) - model_height_b(:,:,z));
                    elseif(z == zlen_wrf)
                        dp = (data_p_b(:,:,z) - data_p_b(:,:,z-1));
                        dz = (model_height_b(:,:,z) - model_height_b(:,:,z-1));
                    else
                        dp = (data_p_b(:,:,z+1) - data_p_b(:,:,z-1));
                        dz = (model_height_b(:,:,z+1) - model_height_b(:,:,z-1));
                    end
                    data_omega_b(:,:,z) = data_w_b(:,:,z).*(dp./dz); % dp is negative, so NEGATIVE omega is upward motion
                end

                 % Align all data fields with customized grid
                pressure_a = align_data(data_p_a,transpose_data,flip_data,false,n_idx,s_idx,e_idx,w_idx);
                pressure_b = align_data(data_p_b,transpose_data,flip_data,false,n_idx,s_idx,e_idx,w_idx);
                
                dims = size(pressure_a);
                
                data_t_a = align_data(data_t_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_t_b = align_data(data_t_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_u_a = align_data(data_u_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_v_a = align_data(data_v_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_w_a = align_data(data_w_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_u_b = align_data(data_u_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_v_b = align_data(data_v_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_w_b = align_data(data_w_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_gph_a = align_data(data_gph_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_gph_b = align_data(data_gph_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_omega_a = align_data(data_omega_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_omega_b = align_data(data_omega_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                data_t_a = interpolate_to_pl(data_t_a,dims,pressure_a,pressure_target);
                data_t_b = interpolate_to_pl(data_t_b,dims,pressure_b,pressure_target);
                data_u_a = interpolate_to_pl(data_u_a,dims,pressure_a,pressure_target);
                data_v_a = interpolate_to_pl(data_v_a,dims,pressure_a,pressure_target);
                data_w_a = interpolate_to_pl(data_w_a,dims,pressure_a,pressure_target);
                data_u_b = interpolate_to_pl(data_u_b,dims,pressure_b,pressure_target);
                data_v_b = interpolate_to_pl(data_v_b,dims,pressure_b,pressure_target);
                data_w_b = interpolate_to_pl(data_w_b,dims,pressure_b,pressure_target);
                data_gph_a = interpolate_to_pl(data_gph_a,dims,pressure_a,pressure_target);
                data_gph_b = interpolate_to_pl(data_gph_b,dims,pressure_b,pressure_target);
                data_omega_a = interpolate_to_pl(data_omega_a,dims,pressure_a,pressure_target);
                data_omega_b = interpolate_to_pl(data_omega_b,dims,pressure_b,pressure_target);

                model_height_a = align_data(model_height_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                model_height_b = align_data(model_height_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);



                % Compute 2D/vertical vorticity: (dv/dx - du/dy)
                data_vort_a = zeros(size(data_u_a));
                data_vort_b = zeros(size(data_u_b));
                data_pva_a = zeros(size(data_vort_a));
                data_pva_b = zeros(size(data_vort_b));
                    
                %if(ismember(upper(data_type),["GPH","WIND","VORT"]))
                    for y = 1:ylen_wrf
                        for x = 1:xlen_wrf
                            if(y == 1 || y == ylen_wrf || x == 1 || x == xlen_wrf)
                                data_vort_a(y,x) = NaN;
                                data_vort_b(y,x) = NaN;
                            else
                                dv = data_v_a(y,x+1) - data_v_a(y,x-1); % In my reference frame, data(y,x), N->S along y+, W->E along x+
                                du = data_u_a(y-1,x) - data_u_a(y+1,x); % We want du and dv to move N and E, so subtract higher y (further S) from lower y (further N)
                                data_vort_a(y,x) = (dv/dx)-(du/dy); % Length scales defined by domain (3km, 9km, x2 for centered finite dif)
                                dv = data_v_b(y,x+1) - data_v_b(y,x-1); 
                                du = data_u_b(y-1,x) - data_u_b(y+1,x); 
                                data_vort_b(y,x) = (dv/dx)-(du/dy);
                            end
                        end
                    end
                    
                    % Compute PVA (Positive Vorticity Advection)
                    gauss_filt_pva = fspecial('gaussian',gauss_width_pva,gauss_sd_pva);
                    data_vort_a_smooth = nanconv(data_vort_a,gauss_filt_pva,'nanout','edge'); 
                    data_vort_b_smooth = nanconv(data_vort_b,gauss_filt_pva,'nanout','edge');
                    data_u_a_smooth = nanconv(data_u_a,gauss_filt_pva,'nanout','edge'); 
                    data_u_b_smooth = nanconv(data_u_b,gauss_filt_pva,'nanout','edge'); 
                    data_v_a_smooth = nanconv(data_v_a,gauss_filt_pva,'nanout','edge'); 
                    data_v_b_smooth = nanconv(data_v_b,gauss_filt_pva,'nanout','edge'); 
                    
                    for y = 2:(ylen_wrf-1)
                        for x = 2:(xlen_wrf-1)
                            if(domain == 2)
                                dx_pva = 6000;
                                dy_pva = 6000;
                            else
                                dx_pva = 18000;
                                dy_pva = 18000;
                            end
                            data_pva_a(y,x) = -1*(data_u_a_smooth(y,x)*(data_vort_a_smooth(y,x+1)-data_vort_a_smooth(y,x-1))/dx_pva + data_v_a_smooth(y,x)*(data_vort_a_smooth(y-1,x)-data_vort_a_smooth(y+1,x))/dy_pva);
                            data_pva_b(y,x) = -1*(data_u_b_smooth(y,x)*(data_vort_b_smooth(y,x+1)-data_vort_b_smooth(y,x-1))/dx_pva + data_v_b_smooth(y,x)*(data_vort_b_smooth(y-1,x)-data_vort_b_smooth(y+1,x))/dy_pva);
                        end
                    end
                    
                %end
                
                % Redirect to main data holder vars
                if(upper(data_type) == "T")
                    data_a = data_t_a;
                    data_b = data_t_b;
                    aligned = true;
                    flattened = true;
                elseif(upper(data_type) == "MSLP")
                    data_a = mslp_a;
                    data_b = mslp_b;
                    aligned = false;
                    flattened = true;
                elseif(upper(data_type) == "GPH")
                    data_a = data_gph_a; 
                    data_b = data_gph_b;
                    aligned = true;
                    flattened = true;
                elseif(upper(data_type) == "VORT")
                    data_a = data_vort_a;
                    data_b = data_vort_b;
                    aligned = true;
                    flattened = true;
                elseif(upper(data_type) == "WIND") % WE DON'T ACTUALLY USE THIS VALUE, JUST NEED A PLACEHOLDER
                    data_a = data_u_a;
                    data_b = data_u_b; 
                    aligned = true;
                    flattened = true;
                elseif(upper(data_type) == "Q")
                    data_a = data_q_a;
                    data_b = data_q_b;
                elseif(upper(data_type) == "OMEGA")
                    data_a = data_omega_a;
                    data_b = data_omega_b;
                    aligned = true;
                    flattened = true;
                elseif(upper(data_type) == "T2M")
                    data_a = ncread(full_filepath_a,data_name); 
                    data_b = ncread(full_filepath_b,data_name);
                    data_a = data_a - 273.15;
                    data_b = data_b - 273.15;
                % If using base reflectivity
                elseif(use_base_refl)
                    load(sprintf(input_format_base_refl,input_path_base_refl_a,data_src_a,'base',domain,member_string,upper(data_type),timestamp_file),'wrf_base_refl');
                    data_a = wrf_base_refl;
                    clearvars wrf_base_refl;
                    load(sprintf(input_format_base_refl,input_path_base_refl_b,data_src_b,'base',domain,member_string,upper(data_type),timestamp_file),'wrf_base_refl');
                    data_b = wrf_base_refl;
                    clearvars wrf_base_refl;
                else     
                    data_a = ncread(full_filepath_a,data_name); % Retrieve data
                    data_b = ncread(full_filepath_b,data_name); % Retrieve data
                end
                
                dims = [ylen_wrf xlen_wrf zlen_wrf];
                
                % All data must be 2D to be processed and plotted. If not
                % composite or interpolated already:
                if(~use_base_refl && ~aligned)
                    data_a = align_data(data_a,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                    data_b = align_data(data_b,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
                    if(is_3d && ~z_composite_data && ~flattened) % Interpolate to desired pressure level
                        data_a = interpolate_to_pl(data_a,dims,pressure_a,pressure_target);
                        data_b = interpolate_to_pl(data_b,dims,pressure_b,pressure_target);
                        %function data_interp = interpolate_to_pl(data_raw,dims_interp,pressure_3d,pressure_target)
                    end
                end
    
                % Store values for next loop
                ens_mem_a(:,:,member_idx) = data_a;
                ens_mem_b(:,:,member_idx) = data_b;
                ens_mem_u_a(:,:,member_idx) = data_u_a;
                ens_mem_u_b(:,:,member_idx) = data_u_b;
                ens_mem_v_a(:,:,member_idx) = data_v_a;
                ens_mem_v_b(:,:,member_idx) = data_v_b;
                ens_mem_vort_a(:,:,member_idx) = data_vort_a;
                ens_mem_vort_b(:,:,member_idx) = data_vort_b;
                ens_mem_pva_a(:,:,member_idx) = data_pva_a;
                ens_mem_pva_b(:,:,member_idx) = data_pva_b;
                ens_mem_gph_a(:,:,member_idx) = data_gph_a;
                ens_mem_gph_b(:,:,member_idx) = data_gph_b;
                ens_mem_omega_a(:,:,member_idx) = data_omega_a;
                ens_mem_omega_b(:,:,member_idx) = data_omega_b;
            end
            
            % Compute thin indices for wind barbs
            thin_y_idxs = 1:thinning_factor:ylen_wrf;
            thin_x_idxs = 1:thinning_factor:xlen_wrf;

            %% Unit conversion
            if(data_type == "SNOWH")
                data_a = data_a.*100;
                data_b = data_b.*100;
                if(use_snow_inches)
                    data_a = data_a.*cm_to_in;
                    data_b = data_b.*cm_to_in;
                end
            end
            
            %% 3D: Apply area mask
            
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
    
            %% 3E: Make dif, compute RMSE, bias, ETS
            
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
    
            % Save error scores in matrix for output at end of script
            if(member_string == "mean")
                idx = num_members+1;
            else
                idx = str2num(member_string);
            end
    
            RMSE_list(idx) = rmse_err;
            bias_list(idx) = bias;
            
            if(upper(data_type) == "WIND")
                data_dif_u = data_u_a - data_u_b;
                data_dif_v = data_v_a - data_v_b;
            end
    
            %% Store for specific plots
            data_storage_a(:,:,member_idx,height_idx) = data_a_compare;
            data_storage_b(:,:,member_idx,height_idx) = data_b_compare;
    
            %% 3G: Plot
            
            if(plot_mean_only && member_idx ~= num_members+1)
                if(show_progress)
                    toc
                    fprintf('Skipping member plots...\n'); % debug statement
                end
                continue;
            elseif(plot_first_and_mean_only && member_idx ~= 1 && member_idx ~= num_members+1)
                if(show_progress)
                    toc
                    fprintf('Skipping all but 1st member plots...\n'); % debug statement
                end
                continue;
            end
    
            if(show_progress)
                toc
                fprintf('Plotting...\n'); % debug statement
            end
            
            gauss_filt = fspecial('gaussian',gauss_width,gauss_sd);
    
            pressure_string = string(pressure_target);
            if(~is_3d && pressure_target == 1000)
                pressure_string = "sfc";
            elseif(~is_3d && pressure_target == 0)
                pressure_string = "VC"; % Vertical composite
            elseif(pressure_target == 1)
                pressure_string = 'X';
            end
            
            for plot_idx = 1:3
                if(plot_idx == 1)
                    if(~plot_raw) continue; end
                    plot_type = 'raw';  % Specify fit, dif, trim, raw
                    data_to_plot = data_a_compare;
                    data_vort = data_vort_a.*10^5;
                    data_pva = data_pva_a.*10^5;
                    data_u = data_u_a;
                    data_v = data_v_a;
                    data_gph = data_gph_a;
                    exp_name_plot = exp_name_a_plot;
                    if(data_type == "REFL")
                        cmap = 'reflmap';
                    elseif(data_type == "SNOWH")
                        cmap = jet_modded;
                    else
                        cmap = 'jet';
                    end
                    plot_filename = sprintf(output_format_var,output_path_large,lower(exp_name_a_plot),lower(exp_name_a_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),pressure_target,timestamp_file);
                    plot_filename_small = sprintf(output_format_var,output_path_small,lower(exp_name_a_plot),lower(exp_name_a_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),pressure_target,timestamp_file);
                elseif(plot_idx == 2)
                    if(~plot_raw) continue; end
                    plot_type = 'raw';  % Specify fit, dif, trim, raw
                    data_to_plot = data_b_compare;
                    data_vort = data_vort_b.*10^5;
                    data_pva = data_pva_b.*10^5;
                    data_u = data_u_b;
                    data_v = data_v_b;
                    data_gph = data_gph_b;
                    exp_name_plot = exp_name_b_plot;
                    if(data_type == "REFL")
                        cmap = 'reflmap';
                    elseif(data_type == "SNOWH")
                        cmap = jet_modded;
                    else
                        cmap = 'jet';
                    end
                    plot_filename = sprintf(output_format_var,output_path_large,lower(exp_name_b_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),pressure_target,timestamp_file);
                    plot_filename_small = sprintf(output_format_var,output_path_small,lower(exp_name_b_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),pressure_target,timestamp_file);
                elseif(plot_idx == 3)
                    if(~plot_dif) continue; end
                    if(compare_to_background)
                        plot_type = 'inc';
                    else
                        plot_type = 'dif';
                    end
                    data_to_plot = data_dif;
                    data_gph = data_gph_a;
                    cmap = 'redblue';
                    plot_filename = sprintf(output_format_var,output_path_large,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),pressure_target,timestamp_file);
                    plot_filename_small = sprintf(output_format_var,output_path_small,lower(exp_name_a_plot),lower(exp_name_b_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type),pressure_target,timestamp_file);
                else
                    fprintf('ERROR: You should not ever see this message. You have a bug in the plot code.\n');
                end
    
                if(isfile(plot_filename) && ~remake_plots) % If shouldn't override existing plots, don't
                    continue;
                end
                
                if(upper(data_type) == "OMEGA")
                    cmap = 'redblue';
                end
    
                % Plot LARGE
                f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                
                if(show_progress && announce_plots)
                    toc
                    fprintf('%s...\n',plot_type); % debug statement
                end
                
                if(plot_type == "raw" && use_contours)
                    if(ismember(upper(data_type),["GPH","WIND","VORT"]) && combo_gph_full)
                        data_vort_smooth = nanconv(data_vort,gauss_filt,'nanout','edge'); % To stop the scale from being messy. Labled on colorbar.
                        h = pcolor(lon_wrf,lat_wrf,data_vort_smooth); % Plot the data
                        set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                        shading interp; % Smooth out plot from grid-boxes
                        c = colorbar('FontSize',axes_font_size); % Make colorbar
                        colormap(flip(hot)); % Set colors
                        if(limit_raw_colorbar)
                            caxis([clim_lower_vort clim_upper_vort]);
                        end
                        ylabel(c,"Cyclonic Vorticity (10^{-5} s^{-1})");
                        hold on;
                        h = quiver(lon_wrf(thin_y_idxs,thin_x_idxs),lat_wrf(thin_y_idxs,thin_x_idxs),data_u(thin_y_idxs,thin_x_idxs),data_v(thin_y_idxs,thin_x_idxs)); % replace with other fileexchange flags?
                    end
                    hold on;
                    data_smooth = nanconv(data_to_plot,gauss_filt,'nanout','edge'); % Keep contours nice and clean
                    [M,h] = contour(lon_wrf,lat_wrf,data_smooth,'EdgeColor',"#0031d9","LineWidth",3,"LevelStep",contour_steps); % Plot the data
                    if(label_contours)
                        clabel(M,h,'FontSize',contour_font_size,'LabelSpacing',250);
                    end
                elseif(upper(data_type) == "WIND")
                    rmse_err = NaN;
                    bias = NaN;
                    if(include_gph)
                        data_smooth = nanconv(data_gph,gauss_filt,'nanout','edge'); % Keep contours nice and clean
                        [M,h] = contour(lon_wrf,lat_wrf,data_smooth,'EdgeColor',"#0031d9","LineWidth",3,"LevelStep",contour_steps); % Plot the data
                        hold on;
                        if(label_contours)
                            clabel(M,h,'FontSize',contour_font_size,'LabelSpacing',250);
                        end
                    end
                    %wind_magnitudes = sqrt(data_u.^2 + data_v.^2);
                    if(plot_idx == 3)
                        data_v = data_dif_v;
                        data_u = data_dif_u;
                        scale_factor = 0.04;
                    else
                        scale_factor = 0.008;
                    end
                    h = quiverwcolorbar(lon_wrf(thin_y_idxs,thin_x_idxs),lat_wrf(thin_y_idxs,thin_x_idxs),data_u(thin_y_idxs,thin_x_idxs),data_v(thin_y_idxs,thin_x_idxs),scale_factor);
                    %quiver(lon_ wrf(thin_y_idxs,thin_x_idxs),lat_wrf(thin_y_idxs,thin_x_idxs),data_u(thin_y_idxs,thin_x_idxs),data_v(thin_y_idxs,thin_x_idxs),'Color',wind_magnitudes); % replace with other fileexchange flags?  
                else
                    if(smooth_data)
                        data_to_plot = nanconv(data_to_plot,gauss_filt,'nanout','edge');
                        clim_mod = 1/gauss_sd;
                    else
                        clim_mod = 1;
                    end
                    h = pcolor(lon_wrf,lat_wrf,data_to_plot); % Plot the data
                    set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                    shading interp; % Smooth out plot from grid-boxes
                    colorbar('FontSize',axes_font_size); % Make colorbar
                    colormap(cmap); % Set colors
    
                    if(plot_type == "dif" || plot_type == "inc")
                        caxis([clim_lower_dif*clim_mod clim_upper_dif*clim_mod]);% Plotting dif
                    else
                        if(limit_raw_colorbar)
                            caxis([clim_lower*clim_mod clim_upper*clim_mod]);% Plotting dif
                        end
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
                    title(sprintf(title_format_raw,exp_name_plot,domain,member_string,upper(data_type),units,pressure_string,timestamp_title),'FontSize',title_font_size); 
                elseif(plot_type == "dif")
                    title(sprintf(title_format_dif,exp_name_a_plot,exp_name_b_plot,domain,member_string,plot_type,upper(data_type),units,pressure_string,timestamp_title,rmse_err,bias),'FontSize',title_font_size); 
                elseif(plot_type == "inc")
                    title(sprintf(title_format_inc,exp_name_a,domain,member_string,upper(data_type),units,pressure_string,timestamp_title),'FontSize',title_font_size); 
                else
                    title(sprintf(title_format_trim,exp_name_plot,plot_type,upper(data_type),units,timestamp_title),'FontSize',title_font_size); 
                end
    
                xlabel('Longitude (deg)','FontSize',label_font_size);
                ylabel('Latitude (deg)','FontSize',label_font_size);
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
                title(sprintf(title_format_small,upper(sprintf('%s-%s',data_src_a,data_src_b)),plot_type,upper(data_type),timestamp_file),'FontSize',title_font_size); % Replace title
                if(show_progress && announce_plots)
                    toc
                    fprintf('Saving small...\n'); % debug statement
                end
                saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
                close('all');
                
        %% 3H: Plot Additional Supporting Fields
                
                % Positive Vorticity Advection (PVA)
                if(include_pva && ismember(upper(data_type),["GPH","WIND","VORT"]) && plot_idx ~=3)
                    plot_type = "PVA";
                    data_type_temp = 'PVA';
                    f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                    h = pcolor(lon_wrf,lat_wrf,data_pva); % Plot the data
                    set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                    shading interp; % Smooth out plot from grid-boxes
                    colorbar('FontSize',axes_font_size); % Make colorbar
                    colormap(flip(hot)); % Set colors
                    caxis([clim_lower_pva clim_upper_pva]);% Plotting dif
                    
                    hold on;
                    borders('continental us','black','linewidth',1); 
                    hold off;
                    if(limit_borders)
                        xlim([w_lim e_lim]);
                        ylim([s_lim n_lim]);
                    end
                    plot_filename = sprintf(output_format_var,output_path_large,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                    plot_filename_small = sprintf(output_format_var,output_path_small,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                    title(sprintf('[%s|d0%d] PVA (10^{-5} s^{-2})| %dmb | %s',exp_name_plot,domain,pressure_target,timestamp_file),'FontSize',title_font_size); % Replace title
                    xlabel('Longitude (deg)','FontSize',label_font_size);
                    ylabel('Latitude (deg)','FontSize',label_font_size);
                    set(gca,'Fontsize',axes_font_size);
                    saveas(gcf,plot_filename); % Save as .png
                    f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
                    saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
                    close('all');
                end
                
                % Plot ECMWF MSLP (pressure) obs
                if(upper(data_type) == "MSLP" && member_idx == (num_members+1) && plot_idx ~= 3)
                    
                    plot_type = "MSLP";
                    data_type_temp = 'MSLP';
                    f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                    [M,h] = contour(obs_lon_full,obs_lat_full,obs_mslp_full,'EdgeColor',"#0031d9","LineWidth",3,"LevelStep",contour_steps);
                    clabel(M,h,'FontSize',contour_font_size,'LabelSpacing',250);
                    hold on;
                    borders('continental us','black','linewidth',1); 
                    if(limit_borders)
                        xlim([w_lim e_lim]);
                        ylim([s_lim n_lim]);
                    end
                    hold off;
                    plot_filename = sprintf(output_format_var,output_path_large,'obs','ecmwf',domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                    plot_filename_small = sprintf(output_format_var,output_path_small,'obs','ecmwf',domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                    title(sprintf('[obs|d%02.f] MSLP(mb) | %s',domain,timestamp_title),'FontSize',title_font_size); 
                    xlabel('Longitude (deg)','FontSize',label_font_size);
                    ylabel('Latitude (deg)','FontSize',label_font_size);
                    set(gca,'Fontsize',axes_font_size);
                    saveas(gcf,plot_filename); % Save as .png
                    f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
                    saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
                    close('all');
                    
                    % Plot 2D L centerpoint scatter
                    plot_type = "LLOC";
                    data_type_temp = 'MSLP';
                    f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                    %h = scatter(L_loc(:,1,plot_idx),L_loc(:,2,plot_idx),25,'MarkerFaceColor',L_dist(:,plot_idx),'MarkerEdgeColor','black');
                    h = scatter(L_loc(:,1,plot_idx),L_loc(:,2,plot_idx),30,L_dist(:,plot_idx),'filled');
                    hold on;
                    h = scatter(L_lon_o,L_lat_o,100,"red","filled","pentagram");
                    labelpoints(L_loc(:,1,plot_idx),L_loc(:,2,plot_idx),[string(1:num_members) "Mean"],'FontSize',axes_font_size-5,'FontWeight','bold');
                    c = colorbar('FontSize',axes_font_size); % Make colorbar
                    colormap(winter); % Set colors
                    ylabel(c,"Distance from member L to obs L (m)");
                    %caxis([clim_lower_pva clim_upper_pva]);% Plotting dif
                    % h = scatter(lon_thin,lat_thin,10,pt_clusters,'filled'); % Plot the data
                    borders('continental us','black','linewidth',1); 
                    if(tight_L)
                        x_limits = [min(min(L_loc(:,1,plot_idx)),L_lon_o)*1.001 max(max(L_loc(:,1,plot_idx)),L_lon_o)*0.999];
                        y_limits = [min(min(L_loc(:,2,plot_idx)),L_lat_o)*0.999 max(max(L_loc(:,2,plot_idx)),L_lat_o)*1.001];
                        x_dif = x_limits(2) - x_limits(1);
                        y_dif = y_limits(2) - y_limits(1);
                        axes_ratio = (x_dif/y_dif);
                        if(axes_ratio < 1)
                            %x_limits(1) = x_limits(1)*(1/axes_ratio);
                        elseif(axes_ratio > 1)
                            %
                        end
                        xlim(x_limits);
                        ylim(y_limits);
                    elseif(limit_borders)
                        xlim([w_lim e_lim]);
                        ylim([s_lim n_lim]);
                    end
                    hold off;
                    plot_filename = sprintf(output_format_var,output_path_large,exp_name_plot,exp_name_plot,domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                    plot_filename_small = sprintf(output_format_var,output_path_small,exp_name_plot,exp_name_plot,domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                    title(sprintf('[%s|d%02.f] L Center Position | %s',exp_name_plot,domain,timestamp_title),'FontSize',title_font_size); 
                    xlabel('Longitude (deg)','FontSize',label_font_size);
                    ylabel('Latitude (deg)','FontSize',label_font_size);
                    set(gca,'Fontsize',axes_font_size);
                    
                    saveas(gcf,plot_filename); % Save as .png
                    f.Position = [fig_x fig_y fig_width_small*1.2 fig_height_small]; % Shrink figure
                    saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
                    close('all');
                
                % Plot 2D freezing line contours
                elseif(upper(data_type) == "T" && plot_fl_np && member_idx == (num_members+1))
                    
                    % Only do it once at the end
                    if(multi_height && height_idx ~= num_heights)
                        continue;
                    end
                    
                    plot_type = 'T0';
                    data_type_temp = 'T0';
                    if(plot_idx == 1)
                        exp_name_plot = exp_name_a_plot;
                    elseif(plot_idx == 2)
                        exp_name_plot = exp_name_b_plot;
                    else
                        close('all');
                        clearvars f h data_to_plot;
                        continue;
                    end
                    
                    output_path_large_temp = sprintf('%s/%s/%s/%s/%s',output_path_base,run_name,upper(data_type_temp),'large');
                    output_path_small_temp = sprintf('%s/%s/%s/%s/%s',output_path_base,run_name,upper(data_type_temp),'small');
    
                    % If output folders do not exist, create them
                    if(~isfolder(output_path_large_temp) || ~isfolder(output_path_small_temp))
                        mkdir(output_path_large_temp);
                        mkdir(output_path_small_temp);
                    end
                    
                    % For a single-height freezing line
                    if(~multi_height)
                        cmap = 'jet';
                        plot_filename = sprintf(output_format_var,output_path_large_temp,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                        plot_filename_small = sprintf(output_format_var,output_path_small_temp,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                          
                        %colors = jet(num_members);
                        colors = reds;
                        
                        % Plot LARGE
                        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                        hold on;
                        borders('continental us','black','linewidth',1); 
                        
                        for mem_sub_idx = 1:num_members
                            if(plot_idx == 1)
                                data_to_plot = data_storage_a(:,:,mem_sub_idx);
                            else
                                data_to_plot = data_storage_b(:,:,mem_sub_idx);
                            end
                            data_smooth = nanconv(data_to_plot,gauss_filt,'nanout','edge'); % Keep contours nice and clean
                            [M,h] = contour(lon_wrf,lat_wrf,data_smooth,[0 0],"LineWidth",1.1,'EdgeColor',colors(mem_sub_idx,:)); % Plot the data
                        end
                        hold off;
    
                        if(limit_borders)
                            xlim([w_lim e_lim]);
                            ylim([s_lim n_lim]);
                        end
    
                        title(sprintf('[%s|d0%d] Ens. Freezing Line (%dmb) | %s',exp_name_plot,domain,pressure_target,timestamp_file),'FontSize',title_font_size); % Replace title
                    else % if multiple heights
                        output_format_fm = '%s/%s_%s_d0%d_sbd%d_%s_%s_%s_%s_%smb_%s.png';
                        
                        plot_filename = sprintf(output_format_fm,output_path_large_temp,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),"multi",timestamp_file);
                        plot_filename_small = sprintf(output_format_fm,output_path_small_temp,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),"multi",timestamp_file);
    
                        contour_storage = struct("h1",[],"h2",[],"h3",[],"h4",[]);
                        
                        % Plot LARGE
                        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                        hold on;
                        borders('continental us','black','linewidth',1); 
                        
                        
                        for height_sub_idx = 1:num_heights
                            
                            switch height_sub_idx
                                case 1
                                    colors = reds;
                                case 2
                                    colors = greens;
                                case 3
                                    colors = blues;
                                case 4
                                    colors = purples;
                                otherwise
                                    error("Too many heights");
                            end
    
                            for mem_sub_idx = 1:num_members
                                if(plot_idx == 1)
                                    data_to_plot = data_storage_a(:,:,mem_sub_idx,height_sub_idx);
                                else
                                    data_to_plot = data_storage_b(:,:,mem_sub_idx,height_sub_idx);
                                end
                                data_smooth = nanconv(data_to_plot,gauss_filt,'nanout','edge'); % Keep contours nice and clean
                                    % DEBUG
                                    fprintf('RGB n%d: [%03.f %03.f %03.f]\n',mem_sub_idx, colors(mem_sub_idx,:)*255);
                                [M,h] = contour(lon_wrf,lat_wrf,data_smooth,[0 0],"LineWidth",1.1,'EdgeColor',colors(mem_sub_idx,:),'DisplayName',sprintf('%dmb',pressure_target_list(height_sub_idx))); % Plot the data
                                if(mem_sub_idx == num_members/2)
                                    switch height_sub_idx
                                        case 1
                                            contour_storage.h1 = h;
                                        case 2
                                            contour_storage.h2 = h;
                                        case 3
                                            contour_storage.h3 = h;
                                        case 4
                                            contour_storage.h4 = h;
                                        otherwise
                                            error("Too many distinct heights");
                                    end
                                end
                            end
                        end
                        hold off;
    
                        if(limit_borders)
                            xlim([w_lim e_lim]);
                            ylim([s_lim n_lim]);
                        end
    
                        title(sprintf('[%s|d0%d] Ens. Freezing Line: Height Range | %s',exp_name_plot,domain,timestamp_file),'FontSize',title_font_size); % Replace title
                        legend([contour_storage.h1 contour_storage.h2 contour_storage.h3 contour_storage.h4]);
                    end
                    
                    xlabel('Longitude (deg)','FontSize',label_font_size);
                    ylabel('Latitude (deg)','FontSize',label_font_size);
                    set(gca,'Fontsize',axes_font_size);
                    hold off;
                    saveas(gcf,plot_filename); % Save as .png
                    f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
                    saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
                    
                end

                close('all');
                clearvars f h data_to_plot;
            
            end
    
            %%
            if(plot_eqt)
                cmap = 'jet';
                plot_type = "era5";
                data_type_temp = "Q";
                exp_name_plot = "era5";
    
                f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                h = pcolor(eqt_lon,eqt_lat,eqt_q(:,:,eqt_p_choice,eqt_time_choice));
                set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                shading interp; % Smooth out plot from grid-boxes
                colorbar('FontSize',axes_font_size); % Make colorbar
    
                colormap(cmap); % Set colors
                %caxis([clim_lower_pva clim_upper_pva]);% Plotting dif
                
                hold on;
                borders('continental us','black','linewidth',1); 
                hold off;
                if(limit_borders)
                    xlim([w_lim e_lim]);
                    ylim([s_lim n_lim]);
                end
                plot_filename = sprintf(output_format_var,output_path_large,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                plot_filename_small = sprintf(output_format_var,output_path_small,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                title(sprintf('[Era5] Q (kg/kg) | %dmb | %s',pressure_target,timestamp_file),'FontSize',title_font_size); % Replace title
                xlabel('Longitude (deg)','FontSize',label_font_size);
                ylabel('Latitude (deg)','FontSize',label_font_size);
                set(gca,'Fontsize',axes_font_size);
                saveas(gcf,plot_filename); % Save as .png
                f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
                saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
                close('all');
    
                data_type_temp = "T";
    
                f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
                h = pcolor(eqt_lon,eqt_lat,eqt_t(:,:,eqt_p_choice,  eqt_time_choice)-273.15);
                set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
                shading interp; % Smooth out plot from grid-boxes
                colorbar('FontSize',axes_font_size); % Make colorbar
    
                colormap(cmap); % Set colors
                %caxis([clim_lower_pva clim_upper_pva]);% Plotting dif
                
                hold on;
                borders('continental us','black','linewidth',1); 
                hold off;
                if(limit_borders)
                    xlim([w_lim e_lim]);
                    ylim([s_lim n_lim]);
                end
                unit_t = sprintf('%sC',char(176));
                plot_filename = sprintf(output_format_var,output_path_large,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                plot_filename_small = sprintf(output_format_var,output_path_small,lower(exp_name_plot),lower(exp_name_plot),domain,subdomain,member_string,bao_short,plot_type,upper(data_type_temp),pressure_target,timestamp_file);
                title(sprintf('[Era5] T (%s) | %dmb | %s',unit_t,pressure_target,timestamp_file),'FontSize',title_font_size); % Replace title
                xlabel('Longitude (deg)','FontSize',label_font_size);
                ylabel('Latitude (deg)','FontSize',label_font_size);
                set(gca,'Fontsize',axes_font_size);
                saveas(gcf,plot_filename); % Save as .png
                f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
                saveas(gcf,sprintf('%s_small.png',extractBefore(plot_filename_small,".png"))); % Save as .png
                close('all');
            end
            
            %% Clear big variables between members
            clearvars rmse_err bias data* -except data_storage* data_type data_name* data_gis_trim data_src* data_wrf_mean_ref;
        end
        end
    
        %% 3I. Write out error value table
    
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
    
            % Make the numbers pretty
            RMSE_list = round(RMSE_list,error_dec);
            bias_list = round(bias_list,error_dec);
    
            % Write
            error_q_table = table(RMSE_list',bias_list','RowNames',err_rowNames,'VariableNames',["RMSE","Bias"]);
            writetable(error_q_table,sprintf(error_table_format,output_path_large,run_name,pressure_target,timestamp_file),'Delimiter','\t','WriteRowNames',true);
        end
        
        clearvars ens_mem* obs_* error_q_table RMSE_list bias_list member_sum* *_colnames data_wrf_mean_ref;
        clearvars rmse_err bias data* -except data_type data_name* data_gis_trim data_src* data_wrf_mean_ref;
        
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
end

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


                