% ptp_analysis.m
%
% Purpose: Generate point-to-point (PTP) diagnostic comparison plots of 
% outputs from IMPACTS WRF to EXRAD and NEXRAD observations
% 
% Author(s): Jon Seibert
% Last updated: 19 Feb 2023
% 
% Inputs: [WRF Outputs].nc, [Observations].png/.nc
% Outputs: [WRF Raw].png, [WRF Fit].png, [Obs Trim].png, [WRF-Obs Dif].png
% Dependencies: reflmap.m, borders.m, find_nn_idx_irregular_exp.m,
%               bilinear_interpolate_irregular_to_grid.m
%
% Process: (For each comparable timestep)
%   1) Read in WRF experiment and/or obs data for a given timestamp
%   2) Calculate basic ensemble mean, use mean for comparisons
%   3) Regrid WRF mean to each obs for direct comparison, if needed
%   4) Calculate dif, RMSE and Bias (mean of dif) for each
%   5) Plot all
%
% TODO:
%   - ADD BACKGROUND INPUT TO COMPARE IT TO ANALYSIS AND/OR OBS (previous
%   timestep) Format: wrfinput_d02_2022-01-29_13:00:00_001, in 1200 folder,
%   for 1300
%   - Reconcile "trim_e" limit from snowband_tracking with current
%   dimensions (NOT NEEDED?)
%
% NOTES:
%  - TIME_COUNT is not robust when crossing month boundaries- replace with
%    manual time-tick count if using it to do so.
%  - NOT ROBUST TO USE IN THE SOUTHERN HEMISPHERE LATITUDES!

% TODO:
%   - Add Equitable Threat Scores (ETS) based on exceedance of a 20-30 dbz
%   threshold

%% 0A. Script Controls

script_name = 'ptp_analysis_ens_full.m';

use_auto_run_name = false;
run_name = 'sd_color_test';

% Local or server?
%local = false; % true = running on local device, false = running on server
local = true;

% Date & time
% Range: 2020-02-07-1300 : 2020-02-07-1800
% But 1300 is the same as 1400? 1400 is same across all experiments
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 7;
start_hour = 14;
end_hour = 14;

exp_choice_a = 1;
exp_choice_b = 0;
use_base_ref = true;
process_all_members = true; % true = all 40 members + mean; false = just ensemble mean

show_plots = true;
show_progress = true; % If true, prints out progress bar/debug messages

%% 0B. General Settings

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020';
    input_path_gis = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    intermediate_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/testing';
    path_to_borders = './borders';
    % Experiment paths
    exp_path_1 = ''; % TESTING
    exp_path_2 = 'CONV/fc/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/';
    input_path_gis = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/GIS';
    intermediate_path = '/storage/home/jjs5895/projects/IMPACTS/intermediate';
    output_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/ptp';
    path_to_borders = '/storage/home/jjs5895/projects/IMPACTS/code/borders'; % Specify path to borders.m
    % Experiment paths
    exp_path_1 = 'AIRCFT/fc';
    exp_path_2 = 'CONV/fc/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
end

addpath(path_to_borders);

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

% Spatial limits of analysis and plot (degrees lat/lon)
w_lim = -79;
e_lim = -69.75;
s_lim = 36;
n_lim = 46;
limit_borders = true; % Whether to apply spatial limits to plot
trim_e = false; % Whether to trim the eastern edge
plot_raw_wrf = false;
plot_sd = true;

% Figure font sizes
title_font_size = 18;
label_font_size = 18;
axes_font_size = 16;

data_type = "refl"; % Variable being analyzed
units = "dBZ";

cmap = 'reflmap'; % Choice of colorbar colormap
clim_lower = -30; % Colorbar limits
clim_upper = 75;
clim_lower_dif = -40;
clim_upper_dif = 40;
sd_max = 10000;

dbz_nodata = -30;
dbz_nodata_gis = -66;

edge_buffer = 42; % Minimum number of points to search from theoretical min in each direction (WARNING: Artifacting may occur at values lower than 36.)

ets_threshold_list = [20,35]; % dBZ threshold for Equitable Threat Scores (ETS)
ets_radius = 9000; % Neighborhood radius for "hits" in threat score calculations

error_dec = 3; % decimal places in error table

%rps_min = 0;
%rps_max = 50;

%% 1A. Settings- WRF Data (2022 Case)

% Basic details
domain = 2;
%num_members = 40;
num_members = 2; % TESTING
data_src_wrf = "wrf"; % Dataset label

% NetCDF retrieval
lon_name_wrf = "XLONG"; 
lat_name_wrf = "XLAT"; 
data_name_wrf = "REFL_10CM"; % Variable being analyzed
time_name_wrf = "Times";

% Note whether transposition or compositing is necessary
transpose_wrf = true;
flip_wrf = true;
composite_wrf = true;

filename_format_wrf = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #
filename_format_alt = 'wrfout_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00'; % domain, year, month, day, hour

%% 1B. Settings- NEXRAD Image Composites via GIS, Iowa Env. Mesonet

data_src_gis = "nex";

file_format_gis = '%s/n0q_%04.f%02.f%02.f%02.f%02.f.png';

lat_gis_raw = 49.9975:-0.005:23.0025;
lon_gis_raw = -125.9975:0.005:-65.0025;
[lon_gis,lat_gis] = meshgrid(lon_gis_raw,lat_gis_raw);

% NOTES: Comes in image format, does not contain own latlon data
% Conversion formula: data_gis = double(imread())*0.5 - 32.5; 
        
%% 2. Preliminary Processing & Setup

if(~show_plots)
    set(groot,'DefaultFigureVisible','off') % Turn off figure popups for local
end

% Announce operating mode
fprintf('Starting %s in %s mode.\n',script_name,mode);

tic; % Start script timer

% Lay out title and filename formats
title_format_fit =   '[%s-%s|d0%d|%s] %s: %s(%s)|%s'; % exp name 1, exp name 2, domain, member #/mean, background/analysis/obs (bao), bao, plot_type, data type, units, timestamp
title_format_dif =   '[%s-%s|d0%d|%s] %s: %s(%s)|%s\nRMSE = %0.3f, Bias = %0.3f'; % exp name 1, exp name 2, domain, member #/mean, background/analysis/obs (bao), bao, plot_type, data type, units, timestamp, err, bias
title_format_trim =  '[%s] %s: %s(%s)|%s'; % data source, plot_type, data type, units, timestamp
title_format_raw =   '[%s|d0%d|%s] Raw: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
title_format_sd =   '[%s|d0%d] SD: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
title_format_small = '[%s-%s] %s-%s-%s'; % exp name 1, exp name 2,plot_type,data type,timestamp
output_format =         '%s/%s_%s_d0%d_%s_%s_%s_%s_%s.png'; % output path, exp name 1, exp name 2, domain, mem/mean, bao_short, plot_type, data type, timestamp
output_format_small =   '%s/%s_%s_d0%d_%s_%s_%s_%s_%s_small.png'; % output path, exp name 1, exp name 2, domain, mem/mean, bao_short, plot_type, data type, timestamp
input_format_base_refl = '%s/%s_%s_d%02.f_%s_%s_%s.mat'; % output path, data source, plot type, domain, member #, data type, datetime
datetime_format_file =  '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

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

bao_2 = 'Analysis';
switch exp_choice_b
    case 0
        compare_to_obs = true;
        bao_2 = 'Obs';
        bao_short = strcat(bao_short,'o');
        input_path_b = input_path_gis;
        exp_name_b = 'OBS';
    case 1
        bao_short = strcat(bao_short,'a');
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_1);
        exp_name_b = exp_name_1;
    case 2
        bao_short = strcat(bao_short,'a');
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_2);
        exp_name_b = exp_name_2;
    case 3
        bao_short = strcat(bao_short,'a');
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_3);
        exp_name_b = exp_name_3;
    case 4
        bao_short = strcat(bao_short,'a');
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_4);
        exp_name_b = exp_name_4;
    case 5
        bao_short = strcat(bao_short,'a');
        input_path_b = sprintf('%s/%s',input_path_base,exp_path_5);
        exp_name_b = exp_name_5;
    otherwise
        error('Unrecognized experiment choice');
end

if(use_auto_run_name)
    run_name = sprintf('%s_%s_%s',lower(exp_name_a),lower(exp_name_b),bao_short);
end
if(show_progress)
    toc
    fprintf('Run designation: %s.\n',run_name);
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
filename_wrf = sprintf(filename_format_wrf,domain,member_idx);
filename_gis = sprintf(file_format_gis,year,month,day,hour,minute);

load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'lon_wrf_trim','lat_wrf_trim'); % Load in preestablished latlon values
lon_wrf = lon_wrf_trim;
lat_wrf = lat_wrf_trim;
clearvars lon_wrf_trim lat_wrf_trim;

dimensions_wrf = size(lon_wrf);
ylen_wrf = dimensions_wrf(1);
xlen_wrf = dimensions_wrf(2);
clearvars dimensions_wrf;

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

% Define table to hold all the error numbers and later write to a file
err_rowNames = 1:(num_members+3);
err_rowNames = string(err_rowNames);
err_rowNames(num_members+1) = "EMean";
err_rowNames(num_members+2) = "Median";
err_rowNames(num_members+3) = "MeanEV";

% Split into two sections: compare to obs, compare to another WRF output

%% 3A. Main Loop: If comparing WRF to OBS
%if(compare_to_obs)

input_path_wrf = input_path_a;    
exp_choice = exp_choice_a;
exp_name = exp_name_a;

if(show_progress)
    toc
    fprintf('Comparison: %s-obs.\n',exp_name_a);
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

%% 3B. For each timestamp being analyzed:
for time_idx = 1:time_count
   
    timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
    timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
    data_src = lower(exp_name);
    
    RMSE_list = zeros(size(err_rowNames));
    bias_list = zeros(size(err_rowNames));
    SD_list = zeros(size(err_rowNames));
    CSI_list_a = zeros(size(err_rowNames));
    CSI_list_b = zeros(size(err_rowNames));
    ETS_list_a = zeros(size(err_rowNames));
    ETS_list_b = zeros(size(err_rowNames));
    brier_list_a = zeros(size(err_rowNames));
    brier_list_b = zeros(size(err_rowNames));
    sd_grid = zeros(ylen_gis,xlen_gis);
    
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

    % Read in each member and compute the mean
    for member_idx = 1:(num_members+1)
        
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

            if(use_base_ref)
                input_path_base_refl = sprintf('%s/BASE_REF/%s',input_path_base,exp_name);
                load(sprintf(input_format_base_refl,input_path_base_refl,data_src,'base',domain,member_string,data_type,timestamp_file),'wrf_base_refl');
                data_wrf = wrf_base_refl;
                clearvars wrf_base_refl;
            elseif(ismember(exp_choice,alt_input_exps))
                filename = sprintf(filename_format_alt,domain,year,month,day,hour);
                data_wrf = ncread(sprintf('%s/%s/%s',input_path_wrf,member_string,filename),data_name_wrf); % Retrieve data
            else
                filename = sprintf(filename_format_wrf,domain,member_idx);       
                data_wrf = ncread(sprintf('%s/%s/%s',input_path_wrf,timestamp_file,filename),data_name_wrf); % Retrieve data
            end

            if(~use_base_ref)
                % If data needs to be composited down to 2D, do so
                if(composite_wrf)
                    data_wrf = max(data_wrf,[],3);
                end

                % If data needs to be transposed to align with lat/lon grid, do so
                if(transpose_wrf) 
                    data_wrf = data_wrf';
                end

                if(flip_wrf)
                    data_wrf = flip(data_wrf);
                end
            end

            member_sum = member_sum + data_wrf;
        
        end
        
        % Determine whether individual member will be fully processed
        process_this_member = (process_all_members || member_string == "mean");
        
        data_wrf_regridded_filename = sprintf('%s/data_%s_%s_d0%d_%s_%s_%s_%s.mat',intermediate_path,exp_name_a,exp_name_b,domain,member_string,bao_short,data_type,timestamp_file);
         
        % Don't recompute if we don't need to
        if(isfile(data_wrf_regridded_filename))
            load(data_wrf_regridded_filename,'data_wrf_regridded');
        else
            data_wrf_regridded = zeros(ylen_gis,xlen_gis); % Composite anyway, no need for vertical

            %% Loop through each GIS point and find corresponding WRF value 
            if(show_progress)
                toc
                fprintf('Regridding...\n'); % debug statement
            end
            for y = 1:ylen_gis
                for x = 1:xlen_gis
                    [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE] = find_nn_idx_irregular(lon_gis_trim(y,x),lat_gis_trim(y,x),lon_wrf,lat_wrf,edge_buffer); % Get nn index values
                    nn_points = [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE];        
                    if(any(nn_points == 0) || any(nn_points(1:4) > xlen_wrf) || any(nn_points(5:8) > ylen_wrf)) % If indices could not be found / reference points are NaN/ goes off edge
                        data_wrf_regridded(y,x) = dbz_nodata; % Lowest dBZ value used on color scale
                    else
                        data_wrf_regridded(y,x) = bilinear_interpolate_irregular_to_grid(y,x,lon_gis_trim,lat_gis_trim,data_wrf,lon_wrf,lat_wrf,nn_points);
                    end
                end
            end

            % Calculate the bounds of the real data for proper comparison with
            % GIS / count how many NaN rows and columns there are
            [first_data_lat,last_data_lat,first_data_lon,last_data_lon] = find_data_edges_nas(dimensions_gis_trim(1),dimensions_gis_trim(2),data_gis_trim);

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
        end
        
        %% Apply mask
        
        if(show_progress)
            toc
            fprintf('Masking...\n'); % debug statement
        end

        load('mask_wrf_on_gis','wrf_gis_mask');

        data_wrf_masked = data_wrf_regridded;
        data_wrf_masked(data_wrf_masked == 0) = 42069;
        data_wrf_masked = data_wrf_masked.*wrf_gis_mask;
        data_wrf_masked(data_wrf_masked == 0) = dbz_nodata_gis;
        data_wrf_masked(data_wrf_masked == 42069) = 0;
        data_wrf_compare = data_wrf_masked;

        data_gis_masked = data_gis_trim;
        data_gis_masked(data_gis_masked == 0) = 42069;
        data_gis_masked = data_gis_masked.*wrf_gis_mask;
        data_gis_masked(data_gis_masked == 0) = dbz_nodata_gis;
        data_gis_masked(data_gis_masked == 42069) = 0;
        data_gis_compare = data_gis_masked;
        clearvars wrf_gis_mask;

        %% Make dif, compute RMSE, bias, ETS
        
        if(show_progress)
            toc
            fprintf('Computing error values...\n'); % debug statement
        end

        data_dif = data_wrf_compare - data_gis_compare;
        data_dif_clean = data_dif(~isnan(data_dif));
        bias = mean(data_dif_clean,'all');
        rmse_err = sqrt(sum(data_dif_clean.^2)/(numel(data_dif_clean)));
        sde = 0;
        sd_grid = sd_grid + (data_dif.^2);
        
        %sd_grid = zeros(ylen_gis,xlen_gis);
        
        % Brier Score
        % For each gridpoint, squared dif of (ensemble probability of >=
        % threshold) to (0/1 obs >= threshold)
        % BS = mean of the above
        % BS = (1/n)*sum((EPi - Oi)^2)
        % Need to save binary of each ensemble member... member_sum_bin?
        % Divide mem_sum_bin by num_members (k) [n = num gridpoints]
        % Then take sqdif of mem_sum_bin to obs at each grid point, take mean

        for ets_idx = 1:2

            ets_threshold = ets_threshold_list(ets_idx);

            % Equitable Threat Score (ETS)
            data_wrf_thresh_list = find(data_wrf_compare>=ets_threshold);
            data_gis_thresh_list = find(data_gis_compare>=ets_threshold);
            data_wrf_bin = zeros(size(data_wrf_compare));
            data_gis_bin = zeros(size(data_gis_compare));
            data_wrf_bin(data_wrf_thresh_list) = 1;
            data_gis_bin(data_gis_thresh_list) = 1;

            if(process_this_member)
            
                a = 0; % Hit (event in both obs and model)
                b = 0; % Miss (event in obs, not in model)
                c = 0; % False alarm (event in model, not in obs)
                d = 0; % Correct miss (neither)

                for y = 1:ylen_gis
                    for x = 1:xlen_gis
                        if(data_gis_bin(y,x) == 1) % If event occurs in obs:
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
            brier_list_a(idx) = mean(power((ens_prob_a - data_gis_bin),2),'all');
            brier_list_b(idx) = mean(power((ens_prob_b - data_gis_bin),2),'all');
            sd_grid = sqrt(sd_grid./((ylen_gis*xlen_gis)-1));
            sde = mean(sd_grid,'all');
        else
            idx = str2num(member_string);
        end

        RMSE_list(idx) = rmse_err;
        bias_list(idx) = bias;
        SD_list(idx) = sde;
        CSI_list_a(idx) = CSI_a;
        CSI_list_b(idx) = CSI_b;
        ETS_list_a(idx) = ETS_a;
        ETS_list_b(idx) = ETS_b;

        %% WRF-GIS.2: Plot

        if(show_progress)
            toc
            fprintf('Plotting...\n'); % debug statement
        end
        
        for plot_idx = 1:3
            if(plot_idx == 1)
                plot_type = 'fit';  % Specify fit, dif, trim, raw
                data_to_plot = data_wrf_regridded;
                cmap = 'reflmap';
            elseif((plot_idx == 2) && (member_idx == 1))
                plot_type = 'trim';
                data_to_plot = data_gis_trim;
                cmap = 'reflmap';
            else
                plot_type = 'dif';
                data_to_plot = data_dif;
                cmap = 'redblue';
            end

            % Plot LARGE
            f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
            % Focus on desired area, remove whitespace
            
            toc
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

            toc
            % Plot state borders
            hold on;
            borders('continental us','black','linewidth',1); 
            hold off;
            
            toc
            if(limit_borders)
                xlim([w_lim e_lim]);
                ylim([s_lim n_lim]);
            end

            toc
            % Apply labels
            if(plot_type == "fit")
                title(sprintf(title_format_fit,exp_name_a,exp_name_b,domain,member_string,plot_type,data_type,units,timestamp_title),'FontSize',title_font_size); 
            elseif(plot_type == "dif")
                title(sprintf(title_format_dif,exp_name_a,exp_name_b,domain,member_string,plot_type,data_type,units,timestamp_title,rmse_err,bias),'FontSize',title_font_size); 
            else
                title(sprintf(title_format_trim,exp_name_b,plot_type,data_type,units,timestamp_title),'FontSize',title_font_size); 
            end
            xlabel('Longitude','FontSize',label_font_size);
            ylabel('Latitude','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);
            toc
            saveas(h,sprintf(output_format,output_path_large,lower(exp_name_a),lower(exp_name_b),domain,member_string,bao_short,plot_type,data_type,timestamp_file)); % Save as .png

            % Plot SMALL
            f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
            title(sprintf(title_format_small,upper(data_src),plot_type,data_type,timestamp_file),'FontSize',title_font_size); % Replace title
            saveas(h,sprintf(output_format_small,output_path_small,lower(exp_name_a),lower(exp_name_b),domain,member_string,bao_short,plot_type,data_type,timestamp_file)); % Save as .png
            %delete(h);
            %delete(f);
            close('all');
            clearvars data_to_plot f h;
            toc
        end

        %% Plot WRF on WRF grid for reference

        if(plot_raw_wrf)
            cmap = 'reflmap';
            plot_type = 'raw';

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
            title(sprintf(title_format_raw,exp_name_a,domain,member_string,data_type,units,timestamp_title),'FontSize',title_font_size); 
            xlabel('Longitude','FontSize',label_font_size);
            ylabel('Latitude','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);
            saveas(h,sprintf(output_format,output_path_large,'wrf',lower(exp_name_a),domain,member_string,bao_short,plot_type,data_type,timestamp_file)); % Save as .png
            %delete(h);
            %delete(f);
            close('all');
        end 
        
        if(plot_sd)
            data_to_plot = sd_grid;
            cmap = 'jet';
            plot_type = 'sd';

            % Plot LARGE
            f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
            % Focus on desired area, remove whitespace
            
            h = pcolor(lon_gis_trim,lat_gis_trim,data_to_plot); % Plot the data
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
            title(sprintf(title_format_sd,exp_name_a,domain,data_type,units,timestamp_title),'FontSize',title_font_size); 
            xlabel('Longitude','FontSize',label_font_size);
            ylabel('Latitude','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);
            saveas(h,sprintf(output_format,output_path_large,lower(exp_name_a),lower(exp_name_b),domain,member_string,bao_short,plot_type,data_type,timestamp_file)); % Save as .png

            % Plot SMALL
            f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
            title(sprintf(title_format_small,upper(data_src),upper(plot_type),data_type,timestamp_file),'FontSize',title_font_size); % Replace title
            saveas(h,sprintf(output_format_small,output_path_small,lower(exp_name_a),lower(exp_name_b),domain,member_string,bao_short,plot_type,data_type,timestamp_file)); % Save as .png
            %delete(h);
            %delete(f);
            close('all');
        end
        
        clearvars f h data_to_plot;
        
        %% Clear variables between members
        clearvars data* -except data_type data_name* data_gis_trim data_sr* sd_grid;
        
    end

    %% Write out error value table

    if(show_progress)
        fprintf('Saving error table...\n'); % debug statement
        toc;
    end
    
    % Compute mean and median values across ensemble
    RMSE_list(num_members+2) = median(RMSE_list(1:num_members));
    RMSE_list(num_members+3) = mean(RMSE_list(1:num_members));
    bias_list(num_members+2) = median(bias_list(1:num_members));
    bias_list(num_members+3) = mean(bias_list(1:num_members));
    SD_list(num_members+2) = median(SD_list(1:num_members));
    SD_list(num_members+3) = mean(SD_list(1:num_members));
    CSI_list_a(num_members+2) = median(CSI_list_a(1:num_members));
    CSI_list_a(num_members+3) = mean(CSI_list_a(1:num_members));
    CSI_list_b(num_members+2) = median(CSI_list_b(1:num_members));
    CSI_list_b(num_members+3) = mean(CSI_list_b(1:num_members));
    ETS_list_a(num_members+2) = median(ETS_list_a(1:num_members));
    ETS_list_a(num_members+3) = mean(ETS_list_a(1:num_members));
    ETS_list_b(num_members+2) = median(ETS_list_b(1:num_members));
    ETS_list_b(num_members+3) = mean(ETS_list_b(1:num_members));

    % Make the numbers pretty
    RMSE_list = round(RMSE_list,error_dec);
    bias_list = round(bias_list,error_dec);
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
    error_q_table = table(RMSE_list',bias_list',ETS_list_a',ETS_list_b',CSI_list_a',CSI_list_b',SD_list',brier_list_a',brier_list_b','RowNames',err_rowNames,'VariableNames',["RMSE","Bias",ETS_colnames(1),ETS_colnames(2),CSI_colnames(1),CSI_colnames(2),"SD",Brier_colnames(1),Brier_colnames(2)]);
    writetable(error_q_table,sprintf('%s/ev_table_%s_%s.txt',output_path_large,run_name,timestamp_file),'Delimiter','\t','WriteRowNames',true);
    clearvars error_q_table RMSE_list bias_list SDE_list CSI_list* ETS_list* brier_list* RPS_list* member_sum* *_colnames;
    
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