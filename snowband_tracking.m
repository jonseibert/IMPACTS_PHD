% snowband_tracking.m
%
% Purpose: Generate physical movement tracks of core snowbands in a winter
% cyclone event
% 
% Author(s): Jon Seibert
% Last updated: 6 Dec 2024
% 
% Inputs: *.nc [Various NetCDF files]
% Outputs: *.png
%
% Usage: Edit sections "0. General Settings" and  "1. Settings" to specify the inputs, run entire script
%
% TODO: Implement more snowband verification checks, generalize file names
% and titles; ALLOW DEMONS MATCHING TO OBS
%  
% Dependencies: borders.m, reflmap.m, haversine_distance.m, fit_ellipse.m,
% image processing toolbox 
%
% NOTES:
%   - Currently designed to work only for WRF output members/mean
%   - WILL BREAK IF CROSSING MONTH BOUNDARIES WITH DATA (find variable
%     'time_count' to adjust)
%   - Assumes lat/lon consistent across all times and members

script_name = 'snowband_tracking.m';

%% 0A. Script Controls
% Meant to be changed often

run_name = 'burst_testing'; % Name of output folder

% Which experiments to use (0 is obs, only for B)
% To compare background and analysis, use the same experiment # for both
% To compare to observations, use 0 for experiment B
exp_choice_a = 3; 
exp_choice_b = 4;

data_type = "REFL"; % Working name of variable being analyzed (REFL, T, MSLP, GPH, VORT, WIND, Q, OMEGA) [+T2M,SNOWH(depth),SNOWNC(grid total per timestep)]
use_base_refl = 1; % Whether to override REFL inputs with precalculated Base Reflectivity values

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

domain = 1; % 2: d02, E coast region; 1: d01, eastern CONUS
subdomain = 0; %1 for E Coast Only subdomain, 2 for N only, 0 for all d02 included

% Vertical levels to use if is 3D
is_3d = 0; % Is the input data 3-dimensional?
pressure_target = 1; %mb [0 = vc, 1 = no height, 1000 = sfc]

z_composite_data = 1; % Unused if use_base_ref is true
multi_height = 0; % Whether to loop over multiple target pressure values
pressure_target_list = [1000,925]; % [925,850,700]; Overrides above pressure_target value if multi_height true

% Plot suite selections
smooth_data = 0; % If true AND the variable isn't already being smoothed, smooth it before plotting (gaussian)

remake_plots = 1; % If true, overrides existing outputs. If false, only plots 'new' files.
plot_mean_only = 0; % Suppresses individual member plots
plot_first_and_mean_only = 1;

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
%z_composite_data = 1;

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

% Filler values
nodata_val = NaN;
nodata_val_gis = NaN;

edge_buffer = 42; % Minimum number of points to search from theoretical min in each direction (WARNING: Artifacting may occur at values lower than 36.)

process_all_members = true; % Whether to compute error table values for each ensemble member
use_auto_run_name = false; % Whether to override manual run name

cm_to_in = 0.3937008;
use_snow_inches = true;

% Constants
Re = 6.3781e6; % Radius pf Earth in m
g = 9.80665; % m/s^2 standard gravity at sea level
R = 287.0600676; % J/kg/K, specific gas constant

% Colorbar settings: overridden by some choices of data variables
limit_raw_colorbar = 1; % Whether to apply specified min and max colorbar values
jet_blue_percent = 0.2; % Cut off coldest (Darkest blues) _% of jet colorbar for ease of analysis

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

filename_format_analysis = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #
filename_format_background = 'wrfinput_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00_%03.f'; % domain, year, month, day, hour, member #
filename_format_background_alt = 'wrfinput_d%02.f_%04.f-%02.f-%02.f_%02.f-00-00_%03.f'; % domain, year, month, day, hour, member #
filename_format_free = 'wrfout_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00'; % domain, year, month, day, hour

%% End modern settings paste

%% 0. General Settings

% Print out processing notes?
print_warnings = false;
progress_bar = false;
show_plots = false;

% Lay out title and filename formats
% Input filenames, plot titles, output filenames, datetime strings
% Compressed = for small plots
filename_format = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #
filename_format_alt = 'wrfinput_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00_%03.f'; % domain, year, month, day, hour, member #
%title_format = '[%s-d%02.f-%s] Neigh. Prob. of %s > %d %s, ROI %d km (%s)'; % Data source (caps), domain, member #, data type, cutoff value, units, ROI/1000, datetime
%title_format_compressed = 'NP>%d|%s|%s'; % cutoff, member #, datetime
output_format = '%s/%s_%s_d%02.f_%s_%s_t%s_%s_c%dk_%s_%s_%s.png'; % output path, data source, plot type, domain, member #, data type, threshold, bounding shape, cutoff dist, testing values string, datetime, demons flag
            % for obs, domain  = d00, mem = 000; bounding shapes: bx, el, nb; 
            % if not using point radius and req, n0k_0; demons flag = d, nd (timestamp should be that of unmodified version)
demons_format = '%s/demons_%s_d%02.f_%s_%s_t%s_c%dk_n%02.fk_%d_%s.mat'; % output path, data source, domain, member #, data type, threshold, cutoff dist, nbh radius, nbh req, datetime
datetime_format_file = '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

%% 1. Settings- Data (WRF 2020 Case)

% Denotes output folder name
run_name = 'air_bugtesting_3';
exp_name = 'AIR';

% Basic details
domain = 2;
%num_members = 40;
num_members = 1; %TESTING
data_src = "air"; % Dataset label
data_type = "refl"; % Variable being analyzed
units = "dBZ";

% NetCDF retrieval
lon_name = "XLONG"; 
lat_name = "XLAT"; 
data_name = "REFL_10CM"; % Variable being analyzed

% Snowband determination vars
snowband_length_req = 250000; % 250 km
snowband_wl_ratio = 1/3;
snowband_portion_cold = 0.5; % Proportion of snowband that has to be in "cold" zone (below freezing)

% Threshold settings
object_threshold = 0; %dBZ; minimum value to be considered a "reflectivity object"
band_threshold = 50; % dBZ: minimum value to be considered a "snowband"
sd_mult = 1.25; % standard deviation value to take for embedded inner snowband
use_sd_threshold = true;

% SLINK Parameters
distance_cutoff = 10000; % (meters) Maximum distance allowed to be bridged when adding a point to a cluster
min_clusters = 1; % Minimum number of clusters allowed
point_req = 100; % Minimum number of points to consider a cluster for snowband status.
nbh_radius = 10000; % Radius  of neighborhood circle [m] (6 km; current clustering cutoff radius is 10 km... try that? link the two?)
nbh_req = 26; % Number of neighbors in radius required to be considered part of a shared band? % Will I be using this?
dt_cutoff_value = 1.2; % Euclidean distance from edge threshold for merging SLINK clusters
weight_centroids = true; % True = centroids weighted by dBZ % of object total

% Display options  
plot_major_axis = true;
plot_centroids = true;
centroid_size = 60;
extend_major_axis = true;

% System version controls
use_v6 = false; % False = use v5 SLINK with neighborhood cutoff; True = V6 with distance transform cutoff
use_demons = true; % Make demons transforms and analyze them too?

% Sensistivity testing controls
radius_list = [10000];
req_list = [26];
%distance_list = [14000];
distance_list = [10000];
dtc_list = [2.4];
%dtc_list = [1.2];

%% 1-B. Setup

if(~show_plots)
    set(groot,'DefaultFigureVisible','off') % Turn off figure popups for local
end

if(use_v6)
    slink_version = 6;
else
    slink_version = 5;
end

% Announce operating mode
fprintf('Starting %s in %s mode using SLINK V%d.\n',script_name,mode,slink_version);
tic

%% 2. Initial Processing

% Directory structure:
% output/np/(run_name)/(variable_name)/(plot_size)
output_path_large = sprintf('%s/%s/%s/%s',output_path_base,run_name,data_type,'large');
output_path_small = sprintf('%s/%s/%s/%s',output_path_base,run_name,data_type,'small');

% Time counters [NOTE: NOT ROBUST TO CROSSING MONTH BOUNDARIES]
max_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 30, 31, 31]; % Number of days in each month
time_count = (end_day - start_day + 1)*24 - (start_hour) - (24 - end_hour) + 1;

% If output folders do not exist, create them
if(~isfolder(output_path_large) || ~isfolder(output_path_small))
    mkdir(output_path_large);
    mkdir(output_path_small);
end

% Specify threshold strings
if(use_sd_threshold)
    thresh_string = sprintf('sd%01d-%02.f',floor(sd_mult),100*mod(sd_mult,1));
else
    thresh_string = sprintf('bt%02.f',band_threshold);
end

% Trim WRF data down to specified domain limits    

%lon_wrf_trim = lon_wrf(n_idx:s_idx,w_idx:e_idx);
%lat_wrf_trim = lat_wrf(n_idx:s_idx,w_idx:e_idx);

load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'lon_wrf_trim','lat_wrf_trim','n_idx','s_idx','e_idx','w_idx'); % Load in preestablished latlon values
    

for nbh_radius = radius_list
for nbh_req = req_list
for distance_cutoff = distance_list
for dt_cutoff_value = dtc_list
    
    if(use_v6)
        testing_value_string = sprintf('dt%01d-%02.f',floor(dt_cutoff_value),100*mod(dt_cutoff_value,1));
        fprintf('Initializing: C = %dkm, DT >= %01.1f\n',distance_cutoff,dt_cutoff_value);
    else
        testing_value_string = sprintf('n%dk_%01.1f',nbh_radius/1000,nbh_req);
        fprintf('Initializing: C = %dkm, NR %dkm, PR %d\n',distance_cutoff,nbh_radius/1000,nbh_req);
    end 

%% 3. MAIN LOOP

% For each member:
for member = 1:num_members 
    
    year = start_year;
    month = start_month;
    day = start_day;
    hour = start_hour;
    
    % For each timestamp being analyzed:
    for time_idx = 1:time_count
        
        
        %% 3A. Read in data 
        
        % Define current loop's time strings
        timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
        timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
        
        fprintf('Processing: Member %d, %s\n',member,timestamp_title);

        % Specify input filename 
        %if(time_idx == 1) %TESTING
            %filename = sprintf(filename_format_alt,domain,year,month,day,hour+1,member); % Special exception to name schema because WRF
        %else
            filename = sprintf(filename_format,domain,member);
        %end
        member_string = sprintf('%03.f',member);

        data = ncread(sprintf('%s/%s/%s',input_path_base,timestamp_file,filename),data_name); % Retrieve data

        % If data needs to be composited down to 2D, do so
        data_trim = align_data(data,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
        dims_trim = size(data_trim);

        % Save original data
        %data_original = data;
        %data = data_trim;
     
        %% 3B. Threshold and set up data
        
        % Base threshold: cut off below [0] dBZ
        flat_list_obj = find(data_trim>object_threshold);
        data_obj = data_trim(flat_list_obj);
        num_points_obj = length(flat_list_obj);
        
        sd = std(data_obj); % Compute standard deviation
        mn = mean(data_obj); % Compute mean
        
        if(use_sd_threshold) % Compute new band threshold dynamically
            band_threshold = mn + (sd*sd_mult); 
        end % Otherwise uses preset value in settings
        
        data_a = data_trim;
        data_a(data_a<=band_threshold) = 0;
        
        %% DEMONS Algorithm Loop
        for demons_id = 1:2 
            
        if((demons_id == 2) && (~use_demons)) % If not using Demons this time, skip demons half of the loop
            continue;
        end
        
        endtime_skip = false;
        demons_filename = sprintf(demons_format,intermediate_path,data_src,domain,member_string,data_type,thresh_string,distance_cutoff/1000,nbh_radius/1000,nbh_req,timestamp_file);
        
        % If the demons transformed version of data_thresh already exists,
        % load it. If not, create it.
        if((exist(demons_filename,'file') == 2))
            load(demons_filename,'D','a_trans');
        elseif(use_demons)
            % If not at the last timestep
            if(time_idx < time_count)
                % Load in next timestep
                filename = sprintf(filename_format,domain,member);
                data_b = ncread(sprintf('%s/%s/%s',input_path_base,sprintf(datetime_format_file,year,month,day,hour+1),filename),data_name); % Retrieve data
                if(z_composite_data)
                    data_b = max(data_b,[],3);
                end
                data_b = data_b(:,:,1);
                if(transpose_data) 
                    data_b = data_b';
                end
                if(flip_data)
                    data_b = flip(data_b);
                end
                
                data_b = data_b(n_idx:s_idx,w_idx:e_idx);
                data_b(data_b<=band_threshold) = 0;
                
                % Perform image registration (phase correlation & demons)
                tformEstimate = imregcorr(data_a,data_b,"translation");
                Rfixed = imref2d(size(data_b));
                a_corr = imwarp(data_a,tformEstimate,"OutputView",Rfixed);
                [D,a_trans] = imregdemons(a_corr,data_b);
                save(demons_filename,'D','a_trans'); % Save as .png
            else
                endtime_skip = true;
            end
        end
        
        % Set flags
        if(demons_id == 1)
            demons_flag = false;
            demons_string = 'nd';
        elseif(endtime_skip)
            continue;
        else
            data_trim = a_trans;
            demons_flag = true;
            demons_string = 'd';
        end
        
        % Get indices of above-threshold points
        %[row_list,col_list] = find(data_trim>band_threshold); % FIND A WAY TO CUT OFF THE SE CORNER REFLECTIVITY??
        flat_list = find(data_trim>band_threshold);
        num_points = length(flat_list); 
        
        % Extract above-threshold points only
        data_thin = data_trim(flat_list);
        lon_thin = lon_wrf_trim(flat_list);
        lat_thin = lat_wrf_trim(flat_list);
        
        % DISTANCE TRANSFORM
        mask = single(zeros(dims_trim(1),dims_trim(2)));
        
        mask(flat_list) = 1;
        dist_trans = bwdist(~mask);
        dist_trans_thin = dist_trans(flat_list);
        %D_norm = D./(max(D,[],"all"));
        
        %% Single Linkage Clustering algorithm
        
        % (Consider storing matrices in .m file)
        % Assumption: Cluster count is initially n, not 0

        % GENERATE DISTANCE MATRIX
            % FOR EACH POINT: LOOP THROUGH EVERY OTHER POINT
                % COMPUTE DISTANCE: SQUARED EUCLIDEAN? HAVERSINE?
                % INPUT INTO MATRIX

        % LOOP: WHILE (NUMBER OF CLUSTERS > THRESHOLD_A)
            % GET INDICES OF SMALLEST DISTANCE
            % ASSIGN ALL MEMBERS OF SMALLER CLUSTER TO INDEX OF LARGER CLUSTER
            % UPDATE DISTANCE MATRIX
                % DISTANCE OF NEW CLUSTER TO ALL OLD LINKS = MIN(DISTANCE OF ADDED POINT(S),DISTANCE OF OLD CLUSTER)
            % DECREASE CLUSTER COUNT BY 1
  
        
        %% 3C. IMPLEMENTATION 6 - Distance Transform instead of Neighbor Count
        % SLINK is the optimal implementation of single-linkage clustering:
        % O(n^2) instead of O(n^3)
        % Track nearest neighbor distances and ids instead of checking
        % matrix each loop
        
        % Note: relevant data are 1D lists of lat and lon values
        % corresponding to the thresholded high-reflectivity data points
        
        if(use_v6)
        
        num_clusters = num_points;
        pt_clusters = 1:num_points; % Tracker for which cluster a point belongs to
        dist_matrix = zeros(num_points,num_points);
        dist_matrix(:,:) = NaN; % n-n distance is 0, but we don't want the algorithm to find "0" as the min distance
        
        nn_idx_list = zeros(num_points,1); % Tracker for nearest neighbor of each point/cluster
        nn_dist_list = zeros(num_points,1); % Tracker for distance between nearest neighbors
        
        %dt_cutoff_value = 0.5; % Don't merge if the point isn't at least this percentile towards the interior
        
        % GENERATE DISTANCE MATRIX
        if(progress_bar)
            tic;
            fprintf('Generating distance matrix...\n')
        end
        % FOR EACH POINT: LOOP THROUGH EVERY OTHER POINT
        for idx_a = 1:num_points
            min_dist = 9999999;
            for idx_b = 1:num_points
                % COMPUTE DISTANCE: SQUARED EUCLIDEAN, INPUT INTO MATRIX
                if(idx_a ~= idx_b) % Don't set diagonals to 0: stay as NaN
                    d = haversine_distance(lon_thin(idx_a),lat_thin(idx_a),lon_thin(idx_b),lat_thin(idx_b));
                    dist_matrix(idx_a,idx_b) = d; dist_matrix(idx_b,idx_a) = d; % Matrix is, by nature, mirrored
                    if(d < min_dist)
                        min_dist = d;
                        min_idx_b = idx_b;
                    end
                end
            end
            nn_idx_list(idx_a) = min_idx_b;
            nn_dist_list(idx_a) = min_dist;
            %if(progress_bar)
            %    if(idx_a == floor(num_points)/4)
            %        fprintf('25%% complete.\n')
            %    elseif(idx_a == floor(num_points/2))
            %        fprintf('50%% complete.\n')
            %    elseif(idx_a == floor(num_points*(3/4)))
            %        fprintf('75%% complete.\n')
            %    end
            %end
        end
        
        if(progress_bar)
            fprintf('Done.\n')
            toc
            tic;
            fprintf('Clustering...\n')
        end
        
        % LOOP: WHILE (NUMBER OF CLUSTERS > MINIMIM CLUSTERS && MINIMUM DISTANCE < DISTANCE CUTOFF)
        while(num_clusters > min_clusters)
            
            % GET INDICES OF SMALLEST DISTANCE
            [min_value,idx_a] = min(nn_dist_list,[],"all","linear");
            idx_b = nn_idx_list(idx_a);
            dt_a = dist_trans_thin(idx_a); dt_b = dist_trans_thin(idx_b);
            
            %sprintf('Min value = %d',min_value) % DEBUG
            if(min_value > distance_cutoff) % If the closest point is too far away, stop early.
                break;
            %elseif(nbh_a < nbh_req || nbh_b < nbh_req) % If at an edge,
            %don't merge? Does this work? would it stop the bits at the
            %edges of bands from ever getting merged? [TRY && INSTEAD OF
            %||- THAT WAY EDGES ONLY DON'T GET MERGED DIRECTLY WITH OTHER
            %EDGES? DOES THAT SOLVE THE PROBLEM?]
            elseif(dt_a < dt_cutoff_value && dt_b < dt_cutoff_value)
                dist_matrix(idx_a,idx_b) = NaN;
                dist_matrix(idx_b,idx_a) = NaN;
                % By definition, b is the nn of a: adjust
                [next_min,next_min_idx] = min(dist_matrix(idx_a,:));
                nn_dist_list(idx_a) = next_min;
                nn_idx_list(idx_a) = next_min_idx;
                % If a is also the nn of b:
                if(nn_idx_list(idx_b) == idx_a)
                    [next_min,next_min_idx] = min(dist_matrix(idx_b,:));
                    nn_dist_list(idx_b) = next_min;
                    nn_idx_list(idx_b) = next_min_idx;
                end
                continue;
            end
            
            % Retrieve cluster numbers
            % Merge (assign all pts in src to dest's cluster number)
            src_pts = find(pt_clusters == pt_clusters(idx_b));
            pt_clusters(src_pts) = pt_clusters(idx_a);
            
            % UPDATE DISTANCE MATRIX:
            % DISTANCE OF NEW CLUSTER TO ALL OLD LINKS = MIN(DISTANCE OF ADDED POINT(S),DISTANCE OF OLD CLUSTER)

            % Edit the row and column corresponding to the absorbing
            % cluster to use the new min distances 
            dist_matrix(idx_a,:) = min(dist_matrix(idx_a,:),dist_matrix(idx_b,:));
            dist_matrix(:,idx_a) = min(dist_matrix(:,idx_a),dist_matrix(:,idx_b));
            dist_matrix(idx_a,idx_a) = NaN; % Distance to itself is 0, but we don't want it merging with itself

            % "Delete" the row and column corresponding to the point
            % being absorbed/the smaller or higher-numbered cluster
            % NOTE: THIS VERSION SETS ALL TO NAN INSTEAD OF DELETING
            dist_matrix(idx_b,:) = NaN;
            dist_matrix(:,idx_b) = NaN; 
            
            %UPDATE NN_DIST_LIST AND NN_IDX_LIST (Av and Ad)
            nn_dist_list(idx_b) = NaN; % Eventually nn_idx_list will be a map of what cluster that cluster merged with. huh.
            [next_min,next_min_idx] = min(dist_matrix(idx_a,:));
            nn_dist_list(idx_a) = next_min;
            nn_idx_list(idx_a) = next_min_idx;
            
            % Find all points that had idx_b as their nearest neighbor and
            % redirect to idx_a
            redir_idx = find(nn_idx_list == idx_b);
            nn_idx_list(redir_idx) = idx_a;
            
            % PLOT?
            
            num_clusters = num_clusters - 1; % Decrement counter- there will always be one fewer independent "cluster"
            
            if(progress_bar)
                if(num_clusters == floor(num_points/4))
                    fprintf('75%% complete.\n')
                elseif(num_clusters == floor(num_points/2))
                    fprintf('50%% complete.\n')
                elseif(num_clusters == floor(num_points*(3/4)))
                    fprintf('25%% complete.\n')
                end
            end
        end
        if(progress_bar) 
            fprintf('Done.\n');
            toc 
            fprintf('Analyzing snowband criteria...\n');
            tic;
        end
        % Output: pt_clusters now contains cluster assignments for all
        % points
        
        else
        
        %% 3C-2. IMPLEMENTATION 5 - SLINK with neighbor-count detection
        % SLINK is the optimal implementation of single-linkage clustering:
        % O(n^2) instead of O(n^3)
        % Track nearest neighbor distances and ids instead of checking
        % matrix each loop
        % 
        % Proposal to better split clusters: track number of
        % above-threshold points within a certain radius of current point
        % being considered; if that value is at a local minimum, that is a
        % split location.
        % Try computing neighbor count all at once during distance matrix
        % computation, checking chains of neighbor counts? i.e. if nn_count
        % down to 4 at current point A, and also 4 at current point B, if
        % either already belongs to a cluster, they mark the border?
        % Diamond shape: i+j displacement no greater than the index radius?
        % Cruder than circle but no extra dist calculations needed.... but
        % we already have all the distances computed, don't we... maybe
        % true circle isn't hard then. Try that.
        
        % Note: relevant data are 1D lists of lat and lon values
        % corresponding to the thresholded high-reflectivity data points 
        
        num_clusters = num_points;
        pt_clusters = 1:num_points; % Tracker for which cluster a point belongs to
        dist_matrix = zeros(num_points,num_points);
        dist_matrix(:,:) = NaN; % n-n distance is 0, but we don't want the algorithm to find "0" as the min distance
        
        nn_idx_list = zeros(num_points,1); % Tracker for nearest neighbor of each point/cluster
        nn_dist_list = zeros(num_points,1); % Tracker for distance between nearest neighbors
        
        %nbh_radius = 20000; % Radius  of neighborhood circle [m] (6 km; current clustering cutoff radius is 10 km... try that? link the two?)
        %nbh_req = 40; % Number of neighbors in radius required to be considered part of a shared band? % Will I be using this?
        nbh_counts = zeros(num_points,1); % Counting variable
        
        % GENERATE DISTANCE MATRIX
        if(progress_bar)
            tic;
            fprintf('Generating distance matrix...\n')
        end
        % FOR EACH POINT: LOOP THROUGH EVERY OTHER POINT
        for idx_a = 1:num_points
            min_dist = 9999999;
            for idx_b = 1:num_points
                % COMPUTE DISTANCE: SQUARED EUCLIDEAN, INPUT INTO MATRIX
                if(idx_a ~= idx_b) % Don't set diagonals to 0: stay as NaN
                    d = haversine_distance(lon_thin(idx_a),lat_thin(idx_a),lon_thin(idx_b),lat_thin(idx_b));
                    dist_matrix(idx_a,idx_b) = d; dist_matrix(idx_b,idx_a) = d; % Matrix is, by nature, mirrored
                    if(d < min_dist)
                        min_dist = d;
                        min_idx_b = idx_b;
                    end
                    if(d <= nbh_radius) % If within neighbor radius, increment counter
                        nbh_counts(idx_a) = nbh_counts(idx_a) + 1;
                    end
                end
            end
            nn_idx_list(idx_a) = min_idx_b;
            nn_dist_list(idx_a) = min_dist;
            %if(progress_bar)
            %    if(idx_a == floor(num_points)/4)
            %        fprintf('25%% complete.\n')
            %    elseif(idx_a == floor(num_points/2))
            %        fprintf('50%% complete.\n')
            %    elseif(idx_a == floor(num_points*(3/4)))
            %        fprintf('75%% complete.\n')
            %    end
            %end
        end
        
        if(progress_bar)
            fprintf('Done.\n')
            toc
            tic;
            fprintf('Clustering...\n')
        end
        
        % LOOP: WHILE (NUMBER OF CLUSTERS > MINIMIM CLUSTERS && MINIMUM DISTANCE < DISTANCE CUTOFF)
        while(num_clusters > min_clusters)
            
            % GET INDICES OF SMALLEST DISTANCE
            [min_value,idx_a] = min(nn_dist_list,[],"all","linear");
            idx_b = nn_idx_list(idx_a);
            nbh_a = nbh_counts(idx_a); nbh_b = nbh_counts(idx_b);
            
            %sprintf('Min value = %d',min_value) % DEBUG
            if(min_value > distance_cutoff) % If the closest point is too far away, stop early.
                break;
            %elseif(nbh_a < nbh_req || nbh_b < nbh_req) % If at an edge,
            %don't merge? Does this work? would it stop the bits at the
            %edges of bands from ever getting merged? [TRY && INSTEAD OF
            %||- THAT WAY EDGES ONLY DON'T GET MERGED DIRECTLY WITH OTHER
            %EDGES? DOES THAT SOLVE THE PROBLEM?]
            elseif(nbh_a < nbh_req && nbh_b < nbh_req)
                dist_matrix(idx_a,idx_b) = NaN;
                dist_matrix(idx_b,idx_a) = NaN;
                % By definition, b is the nn of a: adjust
                [next_min,next_min_idx] = min(dist_matrix(idx_a,:));
                nn_dist_list(idx_a) = next_min;
                nn_idx_list(idx_a) = next_min_idx;
                % If a is also the nn of b:
                if(nn_idx_list(idx_b) == idx_a)
                    [next_min,next_min_idx] = min(dist_matrix(idx_b,:));
                    nn_dist_list(idx_b) = next_min;
                    nn_idx_list(idx_b) = next_min_idx;
                end
                continue;
            end
            
            % Retrieve cluster numbers
            % Merge (assign all pts in src to dest's cluster number)
            src_pts = find(pt_clusters == pt_clusters(idx_b));
            pt_clusters(src_pts) = pt_clusters(idx_a);
            
            % UPDATE DISTANCE MATRIX:
            % DISTANCE OF NEW CLUSTER TO ALL OLD LINKS = MIN(DISTANCE OF ADDED POINT(S),DISTANCE OF OLD CLUSTER)

            % Edit the row and column corresponding to the absorbing
            % cluster to use the new min distances 
            dist_matrix(idx_a,:) = min(dist_matrix(idx_a,:),dist_matrix(idx_b,:));
            dist_matrix(:,idx_a) = min(dist_matrix(:,idx_a),dist_matrix(:,idx_b));
            dist_matrix(idx_a,idx_a) = NaN; % Distance to itself is 0, but we don't want it merging with itself

            % "Delete" the row and column corresponding to the point
            % being absorbed/the smaller or higher-numbered cluster
            % NOTE: THIS VERSION SETS ALL TO NAN INSTEAD OF DELETING
            dist_matrix(idx_b,:) = NaN;
            dist_matrix(:,idx_b) = NaN; 
            
            %UPDATE NN_DIST_LIST AND NN_IDX_LIST (Av and Ad)
            nn_dist_list(idx_b) = NaN; % Eventually nn_idx_list will be a map of what cluster that cluster merged with. huh.
            [next_min,next_min_idx] = min(dist_matrix(idx_a,:));
            nn_dist_list(idx_a) = next_min;
            nn_idx_list(idx_a) = next_min_idx;
            
            % Find all points that had idx_b as their nearest neighbor and
            % redirect to idx_a
            redir_idx = find(nn_idx_list == idx_b);
            nn_idx_list(redir_idx) = idx_a;
            
            % PLOT?
            
            num_clusters = num_clusters - 1; % Decrement counter- there will always be one fewer independent "cluster"
            
            if(progress_bar)
                if(num_clusters == floor(num_points/4))
                    fprintf('75%% complete.\n')
                elseif(num_clusters == floor(num_points/2))
                    fprintf('50%% complete.\n')
                elseif(num_clusters == floor(num_points*(3/4)))
                    fprintf('25%% complete.\n')
                end
            end
        end
        if(progress_bar) 
            fprintf('Done.\n');
            toc 
            fprintf('Analyzing snowband criteria...\n');
            tic;
        end
        % Output: pt_clusters now contains cluster assignments for all
        % points
        
        end
        
        %% 3D. Identify bands among clustered objects
        % 1) Are the bands long enough?
        
        % Holding variables
        cluster_ids = unique(pt_clusters);
        length_flags = zeros(num_clusters,1);
        ratio_flags = zeros(num_clusters,1);
        midpoints = zeros(num_clusters,2);
        major_axis_pts = zeros(num_clusters,4);
        lon_corners = zeros(num_clusters,4);
        lat_corners = zeros(num_clusters,4);
        centroids = zeros(num_clusters,2);
        
        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        colormap(jet);

        % Plot state borders
        hold on;
        borders('continental us','black','linewidth',1); 
  
        % Apply labels
        title(sprintf('[%s|d%02.f|%s] Refl Band Locations: %s',upper(data_src),domain,upper(member_string),timestamp_title),'FontSize',title_font_size); 
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);

        % For each cluster, identify and plot
        for id_num = 1:length(cluster_ids)
            
            %if(id_num == 22)
            %    5; % DEBUG
            %end
            
            %% 3D-1. Retrieve cluster data
            % Retrieve points in cluster
            id = cluster_ids(id_num);
            band_idx = find(pt_clusters == id);
            num_band_pts = length(band_idx);
            
            % Make sure there's at least 10 points in the "band"
            if(num_band_pts <= point_req)
                if(print_warnings)
                    fprintf('WARNING: Band #%02.f has only %02.f member points. Skipping.\n',id_num,num_band_pts);
                end
                continue;
            end
            
            % Get cluster points' latlon
            band_lon = double(lon_thin(band_idx)); % Don't use singles??
            band_lat = double(lat_thin(band_idx));
            band_val = double(data_thin(band_idx));
            
            if(weight_centroids)
                val_sum = sum(band_val);
                band_weights = band_val./val_sum;
                
                centroid_lat = sum(band_weights.*band_lat);
                centroid_lon = sum(band_weights.*band_lon);
            else
                centroid_lat = mean(band_lat);
                centroid_lon = mean(band_lon);
            end
            
            centroids(id_num,:) = [centroid_lon,centroid_lat];
            nonzero_idx = find(centroids(:,1) ~= 0);
            centroids_trim = centroids(nonzero_idx,:);
            
            % Find 2 furthest apart points in cluster, act as major axis ends
            max_len = 0;
            max_dist_ids = [0 0];
            for idx_a = 1:num_band_pts
                for idx_b = 1:num_band_pts
                    d = haversine_distance(band_lon(idx_a),band_lat(idx_a),band_lon(idx_b),band_lat(idx_b));
                    if(d > max_len)
                        max_len = d;
                        max_dist_ids = [idx_a idx_b];
                    end
                end
            end
            
            % Check distance between them
            % Is it >= 250km? 
            if(max_len >= snowband_length_req)
                % If yes, check width, mark as true
                length_flags(id_num) = true;
            end
            
            idx_a = max_dist_ids(1);
            idx_b = max_dist_ids(2);
            lat_a = band_lat(idx_a);
            lon_a = band_lon(idx_a);
            lat_b = band_lat(idx_b);
            lon_b = band_lon(idx_b);
            
            major_axis_pts(id_num,:) = [lon_a,lon_b,lat_a,lat_b];
            %plot([lon_a,lon_b],[lat_a,lat_b]) % Debug
            
            % Determine values of major and minor axes
            % y = mx + b; lat = m*lon + b;
            major_slope = (lat_b-lat_a)/(lon_b-lon_a);
            major_int = lat_a - major_slope*lon_a;
            
            minor_slope = -1/major_slope;
            minor_int_a = lat_a - minor_slope*lon_a;
            minor_int_b = lat_b - minor_slope*lon_b;
            
            %% 3D-2. Find minor axis length
            % Find the two points furthest from the major axis line in each
            % direction: C and D
            pt_line_dist = zeros(num_band_pts,1);
            
            x1 = band_lon(idx_a); x2 = band_lon(idx_b);
            y1 = band_lat(idx_a); y2 = band_lat(idx_b);
                      
            for pt_idx = 1:num_band_pts
                if(pt_idx ~= idx_a && pt_idx ~= idx_b)
                    x0 = band_lon(pt_idx);
                    y0 = band_lat(pt_idx);
                    d = abs(((x2-x1)*(y1-y0))-((x1-x0)*(y2-y1)))/sqrt(((x2-x1)^2)+((y2-y1)^2));
                    pt_line_dist(pt_idx) = d;
                end
            end
            
            dist_temp = pt_line_dist; % Make copy to mess with
            
            % Loop down the highest distances until we find two on opposing sides of the major axis line
            found_edge1 = false;
            found_edge2 = false;
            counter = 1;
            while(~(found_edge1 && found_edge2)) 
                [max_pl_dist,max_pl_idx] = max(dist_temp);
                
                % If the max remaining distance is 0, didn't find D.
                % Distance to A and B will always be 0.
                if(max_pl_dist == 0)
                    break;
                end
                
                pt_lat = band_lat(max_pl_idx);
                pt_lon = band_lon(max_pl_idx);
                % y = mx+b; lat = major_slope*lon + major_int;
                % x = (y-b)*(1/m); lon = (lat - major_int)*(1/major_slope)
                maj_lat_1 = pt_lat; % Pt 1: same latitude, find W/E lon
                maj_lat_2 = major_slope*pt_lon + major_int; % Pt 2: same lon, find N/S lat
                maj_lon_1 = (pt_lat - major_int)*(1/major_slope);
                maj_lon_2 = pt_lon;
                
                if(pt_lon > maj_lon_1) % Is pt E of line?
                    pt_x = 1;
                else % Is pt W of line?
                    pt_x = -1;
                end
                
                if(pt_lat > maj_lat_2) % Is point N of line?
                    pt_y = 1;
                else % Is point S of line?
                    pt_y = -1;
                end
                
                if(found_edge1)
                    if(((pt_y == -1*pt_y_1) && (pt_x == -1*pt_x_1))) % If new point is not on the other side, skip to next
                        edge2_idx = max_pl_idx;
                        found_edge2 = true;
                    end
                else
                    edge1_idx = max_pl_idx;
                    pt_y_1 = pt_y;
                    pt_x_1 = pt_x;
                    found_edge1 = true;
                end
                
                dist_temp(max_pl_idx) = 0; % Set max dist in temp copy to 0 to find next largest dist
                counter = counter+1;
                % If checked every point and not found one on the other side, there are none.
                % In that event, set that point to the midpoint of the
                % major axis.
                if(counter >= num_band_pts) 
                    break;
                end
            end
            
            idx_c = edge1_idx; % These both need to use the major axis slope, intercept defined by their location
            lon_c = band_lon(idx_c);
            lat_c = band_lat(idx_c);
            
            if(found_edge2) % If found points on both sides of the major axis:
                idx_d = edge2_idx;
                lon_d = band_lon(idx_d);
                lat_d = band_lat(idx_d);
            else % If all points are on the same side of the major axis: create fake point D on major axis
                idx_d = 0;
                lon_d = (lon_a + lon_b)/2;
                lat_d = (lat_a + lat_b)/2;
            end
            
            major_int_c = lat_c - major_slope*lon_c;
            major_int_d = lat_d - major_slope*lon_d;
            
            %% 3D-3. Correct major axis length: find outer extreme points
            % Find the two points furthest from the minor axis line in each
            % direction: A2 and B2
            
            if(extend_major_axis)
                %fprintf('Extending...'); % DEBUG
            
            
            % Define true minor axis
            % minor_slope
            midpoint_lat = (lat_a + lat_b)/2;
            midpoint_lon = (lon_a + lon_b)/2;
            
            % Equation for line of minor axis:
            % lat = m*lon + b; b = lat - m*lon;
            midpoint_int = midpoint_lat - minor_slope*midpoint_lon;
            
            % Minor axis:         lat = minor_slope*lon + midpoint_int;
            % Equation of side C: lat = major_slope*lon + major_int_c;
            % Equation of side D: lat = major_slope*lon + major_int_d;
            lon_mc = (midpoint_int - major_int_c)/(major_slope-minor_slope);
            lat_mc = major_slope*lon_mc + major_int_c;
            lon_md = (midpoint_int - major_int_d)/(major_slope-minor_slope);
            lat_md = major_slope*lon_md + major_int_d;
            
            pt_line_dist = zeros(num_band_pts,1);
            
            x1 = lon_mc; x2 = lon_md;
            y1 = lat_mc; y2 = lat_md;
                      
            for pt_idx = 1:num_band_pts
                %if(pt_idx ~= idx_a && pt_idx ~= idx_b)
                    x0 = band_lon(pt_idx);
                    y0 = band_lat(pt_idx);
                    d = abs(((x2-x1)*(y1-y0))-((x1-x0)*(y2-y1)))/sqrt(((x2-x1)^2)+((y2-y1)^2)); % Distance from a point to a line
                    pt_line_dist(pt_idx) = d;
                %end
            end
            
            dist_temp = pt_line_dist; % Make copy to mess with
            
            % Loop down the highest distances until we find two on opposing sides of the minor axis line
            found_edge1 = false;
            found_edge2 = false;
            counter = 1;
            while(~(found_edge1 && found_edge2)) 
                [max_pl_dist,max_pl_idx] = max(dist_temp);
                
                % If the max remaining distance is 0, didn't find D.
                % Distance to axis points will always be 0.
                if(max_pl_dist == 0)
                    break;
                end
                
                pt_lat = band_lat(max_pl_idx);
                pt_lon = band_lon(max_pl_idx);
                % y = mx+b; lat = major_slope*lon + major_int;
                % x = (y-b)*(1/m); lon = (lat - major_int)*(1/major_slope)
                      
                min_lat_1 = pt_lat; % Pt 1: same latitude, find W/E lon
                min_lat_2 = minor_slope*pt_lon + midpoint_int; % Pt 2: same lon, find N/S lat
                min_lon_1 = (pt_lat - midpoint_int)*(1/minor_slope);
                min_lon_2 = pt_lon;
                
                if(pt_lon > min_lon_1) % Is pt E of line?
                    pt_x = 1;
                else % Is pt W of line?
                    pt_x = -1;
                end
                
                if(pt_lat > min_lat_2) % Is point N of line?
                    pt_y = 1;
                else % Is point S of line?
                    pt_y = -1;
                end
                
                if(found_edge1)
                    if(((pt_y == -1*pt_y_1) && (pt_x == -1*pt_x_1))) % If new point is not on the other side, skip to next
                        edge2_idx = max_pl_idx;
                        found_edge2 = true;
                    end
                else
                    edge1_idx = max_pl_idx;
                    pt_y_1 = pt_y;
                    pt_x_1 = pt_x;
                    found_edge1 = true;
                end
                
                dist_temp(max_pl_idx) = 0; % Set max dist in temp copy to 0 to find next largest dist
                counter = counter+1;
                % If checked every point and not found one on the other side, there are none.
                % In that event, set that point to the midpoint of the
                % major axis.
                if(counter >= num_band_pts) 
                    break;
                end
            end
            
            idx_a2 = edge1_idx; % These both need to use the major axis slope, intercept defined by their location
            lon_a2 = band_lon(idx_a2);
            lat_a2 = band_lat(idx_a2);
            
            if(found_edge2) % If found points on both sides of the major axis:
                idx_b2 = edge2_idx;
                lon_b2 = band_lon(idx_b2);
                lat_b2 = band_lat(idx_b2);
            else % If all points are on the same side of the major axis: create fake point D on minor axis
                idx_b2 = 0;
                lon_b2 = (lon_c + lon_d)/2;
                lat_b2 = minor_slope*lon_b2 + midpoint_int;
            end
            
            % REPLACE A AND B INTERCEPTS WITH ADJUSTED BORDERS
            minor_int_a = lat_a2 - minor_slope*lon_a2;
            minor_int_b = lat_b2 - minor_slope*lon_b2;
            
            end
            
            %% 3D-4. Compute intersections and bounding box, fit model
            
            % y = mx + b; b = (y-mx);
            % Major axis = intercepts A and B
            % Minor axis = undefined
            % N edge = minor slope through B; y = minor_slope*x + minor_int_b;
            % S edge = minor slope through A; y = minor_slope*x + minor_int_a;
            % E edge = major slope through C
            % W edge = major slope through D (or vv with C)
            
            
            % Find intersection points
            % y = major_slope*x + major_int_c; y = major_slope*x + major_int_d;
            % major_slope*lon + major_int_c = minor_slope*lon + minor_int_a;
            % lon(major-minor) = minor_ia - major_ic;
            % lon = (minor_ia - major_ic)/(major-minor);
            lon_ac = (minor_int_a - major_int_c)/(major_slope-minor_slope);
            lat_ac = major_slope*lon_ac + major_int_c;
            lon_ad = (minor_int_a - major_int_d)/(major_slope-minor_slope);
            lat_ad = major_slope*lon_ad + major_int_d;
            lon_bc = (minor_int_b - major_int_c)/(major_slope-minor_slope);
            lat_bc = major_slope*lon_bc + major_int_c;
            lon_bd = (minor_int_b - major_int_d)/(major_slope-minor_slope);
            lat_bd = major_slope*lon_bd + major_int_d;
            % Box = lines between AC, AD, BC, BD
            
            max_width = haversine_distance(lon_ac,lat_ac,lon_ad,lat_ad);
            
            % Is the ratio of length:width 3:1 or greater?
            if(max_width <= (max_len*snowband_wl_ratio))
                ratio_flags(id_num) = true;
            end
            
            lon_corners(id_num,:) = [lon_ac,lon_ad,lon_bd,lon_bc];
            lat_corners(id_num,:) = [lat_ac,lat_ad,lat_bd,lat_bc];
            
            % Fit quadratic line to band
            mdl = fitlm(band_lon,band_lat,'quadratic');
            a = mdl.Coefficients.Estimate(2);
            b = mdl.Coefficients.Estimate(3);
            c = mdl.Coefficients.Estimate(1);
            
            test_fit = a*band_lon + b*power(band_lon,2) + c;
            %plot(mdl);
            h = plot(band_lon,test_fit);
            
            % Find midpoint to use as "centroid"
            %midpoint_lat = (lat_a+lat_b)/2;
            mid_lon_fit = (lon_a+lon_b)/2;
            mid_lat_fit = a*mid_lon_fit + b*power(mid_lon_fit,2) + c;
            midpoints(id_num,:) = [mid_lon_fit, mid_lat_fit];
        end
        
        % Save tracking stats for snowband_dif
        stats_filename = sprintf(output_format,intermediate_path,data_src,'stats',domain,member_string,data_type,thresh_string,'nb',distance_cutoff/1000,testing_value_string,timestamp_file,demons_string);
        if(demons_string == 'd')
            D_save = D;
        else
            D_save = NaN;
        end
        save(stats_filename,'pt_clusters','cluster_ids','centroids_trim','data_trim','D_save');
        
        %% 4A. Plotting - Cluster fits
        
        % Focus on desired area, remove whitespace
        if(limit_borders)
            xlim([w_lim e_lim]);
            ylim([s_lim n_lim]);
        else
            xlim([min(lon_thin) max(lon_thin)]);
            ylim([min(lat_thin) max(lat_thin)]);
        end
        
        nonzero_idx = find(midpoints(:,1) ~= 0);
        midpoints_trim = midpoints(nonzero_idx,:);
        
        % Plot band midpoints
        h = scatter(midpoints_trim(:,1),midpoints_trim(:,2),'b',"filled");
        
        plot_type = 'cluster_fit';
        bounding_shape = 'nb';

        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png

        
            % output_format = '%s/%s_%s_d%02.f_%s_%s_t%s_%s_c%dk_n%02.fk_%d_%s_%s.png'; % output path, data source, plot type, domain, member #, data type, threshold, bounding shape, cutoff dist, nbh radius, nbh req, datetime, demons flag
            % for obs, domain  = d00, mem = 000; bounding shapes: bx, el, nb; 
            % if not using point radius and req, n0k_0; demons flag = d, nd (timestamp should be that of unmodified version)
        
        % Plot bounding boxes
        for id_num = 1:num_clusters
            box_lon = [lon_corners(id_num,1),lon_corners(id_num,2),lon_corners(id_num,3),lon_corners(id_num,4),lon_corners(id_num,1)];
            box_lat = [lat_corners(id_num,1),lat_corners(id_num,2),lat_corners(id_num,3),lat_corners(id_num,4),lat_corners(id_num,1)];
            h = plot(box_lon,box_lat,'r','LineWidth',3);
            if(plot_major_axis)
                h = plot([major_axis_pts(1),major_axis_pts(2)],[major_axis_pts(3),major_axis_pts(4)]);
            end
        end
        
        bounding_shape = 'bx';
        
        title(sprintf('[%s|d%02.f|%s] Refl Band Locations (Box): %s',upper(data_src),domain,upper(member_string),timestamp_title),'FontSize',title_font_size); 
        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png

            
            
        hold off;
        close("all");
        
        is_snowband = (length_flags & true);  % WIP
        
        %% 4B. Plotting - Cluster groups
        
        % PLOT CLUSTERS INDIVIDUALLY AND DRAW OVALS AROUND THEM?
        
        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        h = scatter(lon_thin,lat_thin,10,pt_clusters,'filled'); % Plot the data
        ax = gca;
        c = colorbar('FontSize',axes_font_size); % Make colorbar
        %colormap(lines); % Set colors
        colormap(jet);

        % Plot state borders
        hold on;
        borders('continental us','black','linewidth',1); 

        % Focus on desired area, remove whitespace
        if(limit_borders)
            xlim([w_lim e_lim]);
            ylim([s_lim n_lim]);
        else
            xlim([min(lon_thin) max(lon_thin)]);
            ylim([min(lat_thin) max(lat_thin)]);
        end
        
        plot_type = 'cluster_shp';
        bounding_shape = 'nb';

        % Apply labels
        title(sprintf('[%s|d%02.f|%s] Refl Band Clusters: %s',upper(data_src),domain,upper(member_string),timestamp_title),'FontSize',title_font_size); 
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);

        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png

        % Plot mean ellipses
        for id_num = 1:num_clusters
            if(is_snowband(id_num))
                band_idx = find(pt_clusters == cluster_ids(id_num));
                fit_ellipse(lon_thin(band_idx),lat_thin(band_idx),ax);
            end
        end
        
        hold off; 
        
        bounding_shape = 'el';
        
        title(sprintf('[%s|d%02.f|%s] Refl Band Clusters (Mean Ellipse): %s',upper(data_src),domain,upper(member_string),timestamp_title),'FontSize',title_font_size); 
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);
            
        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png

        close("all");
        
        % Replot base clusters with boxes
        
        % Save cluster files for later use
        %save(sprintf('intermediate/clusters_%s_%s',member_string,timestamp_file),"data_thresh","lon_wrf_trim","lat_wrf_trim","pt_clusters");
        
        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        h = scatter(lon_thin,lat_thin,10,pt_clusters,'filled'); % Plot the data
        ax = gca;
        c = colorbar('FontSize',axes_font_size); % Make colorbar
        %colormap(lines); % Set colors
        colormap(jet);

        % Plot state borders
        hold on;
        borders('continental us','black','linewidth',1); 

        % Focus on desired area, remove whitespace
        if(limit_borders)
            xlim([w_lim e_lim]);
            ylim([s_lim n_lim]);
        else
            xlim([min(lon_thin) max(lon_thin)]);
            ylim([min(lat_thin) max(lat_thin)]);
        end

        % Apply labels
        %title(sprintf('[%s|d%02.f|%s] Refl Band Clusters: %s',upper(data_src),domain,upper(member_string),timestamp_title),'FontSize',title_font_size); 
        %xlabel('Longitude','FontSize',label_font_size);
        %ylabel('Latitude','FontSize',label_font_size);
        %set(gca,'Fontsize',axes_font_size);
        %if(use_v6)
        %    saveas(h,sprintf('%s/clusters_%s_d%02.f_%s_%s_%s.png',output_path_large,data_src,domain,member_string,data_type,timestamp_file)); % Save as .png
        %else
        %    %saveas(h,sprintf('%s/clusters_%s_d%02.f_%s_%s_t%s_%s_c%dk_n%02.fk_%d.png',output_path_large,data_src,domain,member_string,data_type,thresh_string,timestamp_file,distance_cutoff/1000,nbh_radius/1000,nbh_req)); % Save as .png
        %    saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,nbh_radius/1000,nbh_req,timestamp_file,demons_string)); % Save as .png
        %end
        
        % Plot bounding boxes
        for id_num = 1:num_clusters
            box_lon = [lon_corners(id_num,1),lon_corners(id_num,2),lon_corners(id_num,3),lon_corners(id_num,4),lon_corners(id_num,1)];
            box_lat = [lat_corners(id_num,1),lat_corners(id_num,2),lat_corners(id_num,3),lat_corners(id_num,4),lat_corners(id_num,1)];
            h = plot(box_lon,box_lat,'r','LineWidth',3);
            if(plot_major_axis)
                h = plot([major_axis_pts(1),major_axis_pts(2)],[major_axis_pts(3),major_axis_pts(4)]);
            end
        end
        
        nonzero_idx = find(centroids(:,1) ~= 0);
        centroids_trim = centroids(nonzero_idx,:);
        
        if(plot_centroids)
            h = scatter(centroids_trim(:,1),centroids_trim(:,2),centroid_size,'black',"filled","pentagram");
        end
        
        bounding_shape = 'bx';
        
        title(sprintf('[%s|d%02.f|%s] Refl Band Clusters (Box): %s',upper(data_src),domain,upper(member_string),timestamp_title),'FontSize',title_font_size);
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);
        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png

        hold off;
        close("all");

        % Plot SMALL
        %f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
        %title(sprintf(title_format_compressed,cutoff,member_string,timestamp_title),'FontSize',title_font_size);
        %saveas(h,sprintf(output_format,output_path_small,data_src,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png
        %close('all');
        
        save(sprintf('intermediate/snowband_tracking_stats_%s_%d00_%s',lower(exp_name),hour,demons_string),'pt_clusters','cluster_ids','centroids_trim');
        
        
        %% DISTANCE TRANSFORM
        
        mask = single(zeros(dims_trim(1),dims_trim(2)));
        mask(flat_list) = 1;
        %imshow(mask);
        dist_trans = bwdist(~mask);
        %D_norm = D./(max(D,[],"all"));

        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        h = pcolor(lon_wrf_trim,lat_wrf_trim,dist_trans);
        set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
        shading interp; % Smooth out plot from grid-boxes
        c = colorbar('FontSize',axes_font_size); % Make colorbar
        colormap('gray'); % Set colors
        
        if(limit_borders)
            xlim([w_lim e_lim]);
            ylim([s_lim n_lim]);
        else
            xlim([min(lon_thin) max(lon_thin)]);
            ylim([min(lat_thin) max(lat_thin)]);
        end
        
         % Apply labels
        title(sprintf('[%s|d%02.f|%s] Refl Distance Transform: %s',upper(data_src),domain,upper(member_string),timestamp_title),'FontSize',title_font_size); 
        xlabel('Longitude','FontSize',label_font_size);
        ylabel('Latitude','FontSize',label_font_size);
        ylabel(c,'Normalized Euclidean Distance to Reflectivity Object Edge','FontSize',label_font_size)
        set(gca,'Fontsize',axes_font_size);
        
        plot_type = 'dist_trans_nb';
        bounding_shape = 'nb';
        
        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png

        hold on;
        borders('continental us','color','#20891e','linewidth',1); 
        hold off;
        
        plot_type = 'dist_trans';
        
        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png

        % Also plot base refl for reference
        
        if(demons_flag) % Make sure demons transformed output without the cyan background
            data_trim(data_trim == 0) = -30;
        end
        
        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        h = pcolor(lon_wrf_trim,lat_wrf_trim,data_trim); % Plot the data
        set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
        shading interp; % Smooth out plot from grid-boxes
        c = colorbar('FontSize',axes_font_size); % Make colorbar
        colormap(cmap); % Set colors
        caxis([clim_lower clim_upper]);

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
        
        plot_type = 'raw';
        bounding_shape = 'nb';

        title_format_raw = '[%s|d0%d|%s] Raw: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
        title(sprintf(title_format_raw,exp_name,domain,member_string,data_type,units,timestamp_title),'FontSize',title_font_size);
        
        saveas(h,sprintf(output_format,output_path_large,data_src,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png

        if(progress_bar) 
            fprintf('Done.\n');
            toc 
        end
        
        end % End demons loop
        
        %% 5. Loop increment controls
        
        % Time increment
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

end
end
end
end

set(groot,'DefaultFigureVisible','on');
fprintf('%s run complete.\n',script_name);
toc