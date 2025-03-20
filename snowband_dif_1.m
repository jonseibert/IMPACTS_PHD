% snowband_dif_obs.m
%
% Purpose: Compute object centroid distances between snowband objects
% 
% Author(s): Jon Seibert
% Last updated: 4 Mar 2025
% 
% Inputs: *.nc [Various NetCDF files]
% Outputs: *.png
%
% Usage: Edit sections "0. General Settings" and  "1. Settings" to specify the inputs, run entire script
%
% TODO: 
%   - Make version that uses D_save to compute correspondances between d
%   and nd versions of the same file, then take distances between WRF ND
%   and OBS ND!
%   - Add point weighting
%  
% Dependencies: borders.m, reflmap.m, haversine_distance.m, fit_ellipse.m,
% image processing toolbox 
%
% NOTES:
%   - ASSUMES "A" IS OBS, AND IS _ND, WITH "B" AS EXP, _D
%   - If "use_weighted_distances" is true, output units are distance per
%   point (in above-threshold cluster meeting minimum size requirements, as
%   a proportion of all others meeting said criteria)


script_name = 'snowband_dif_obs.m';

%% 0A. Script Controls

% Denotes output folder name
run_name = 'thresh_dif_5_pred';

exp_choice = 1; %  

use_weighted_distances = 1; % Weight by num points in each object (1) or by number of objects (0)

use_transformed_distances = 0; % Use distances from WRF Demons Transformation -> OBS (1), or backtrack to retrieve distances from WRF->OBS w/o Demons, using output D matrix (0)

% -------------------
use_base_refl = 1; % Whether to use vertical composite (0) or simulated base (1) reflectivity

compare_against = 0; % obs

%reg_against = 1; % 0 = obs; if not obs, compares member to member, could be at a dif time or dif experiment [NOT IMPLEMENTED]
%time_forward = 0; % How many timesteps forward to compare to [NOT IMPLEMENTED]

% Date & time
start_year = 2020;
end_year = 2020;
start_month = 2;
end_month = 2;
start_day = 7;
end_day = 7;
start_hour = 14;
end_hour = 18;
minute = 0;

num_members = 40;
use_member_list = 1;
custom_member_list = 1:2;

thresh_list = [20]; % 20, 35, 50

%%

apply_ocean_mask = 1; % Cut off the bits SE of the ocean exclusion line

% Threshold settings
object_threshold = 0; %dBZ; minimum value to be considered a "reflectivity object"
%band_threshold = 40; % dBZ; minimum value to be considered a "snowband" (Overridden)
%band_threshold_b = 35; % dBZ threshold for Demons target
band_threshold_b_adjustment = 0; % modifier for threshold for demons target (obs)
sd_mult = 1.25; % Standard deviation value to take for embedded inner snowband
use_sd_threshold = 0; % Flag: Override with adaptive threshold = mean + sd_mult*SD?

% SLINK Parameters
%distance_cutoff = 10000; % (meters) Maximum distance allowed to be bridged when adding a point to a cluster
min_clusters = 1; % Minimum number of clusters allowed
point_req = 20; % Minimum number of points to consider a cluster for snowband status.
%nbh_radius = 10000; % Radius  of neighborhood circle [m] (6 km; current clustering cutoff radius is 10 km... try that? link the two?)
%nbh_req = 26; % Number of neighbors in radius required to be considered part of a shared band? % Will I be using this?
%dt_cutoff_value = 1.2; % Euclidean distance from edge threshold for merging SLINK clusters

% Sensistivity testing controls

radius_list = [10000]; % List of nbh radii
%req_list = [20, 26, 32]; % List of nbh point count requirements [20 26 32]
req_list = [26]; % List of nbh point count requirements [20 26 32]
%distance_list = [14000]; % List of distance cutoffs
distance_list = [10000];

%override_diff = 1;
make_boxplot = 1;
make_condensed_boxplots = 1;

override_saved_stats = 0;
remake_plots = 1;
make_empty_plots = 0;

% Print out processing notes?
show_progress = 1;
show_loading_msg = 0;
print_warnings = 0;
show_plots = 1;
show_debug = 0;

y_neg_scaling = 1/8.5; % Ratio to make sure the boxplots go far enough below 0 to separate the whiskers from the ticks

%% 0B. Environment Detection

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

%% 0C. General Settings

if(local)
    % Figure font sizes
    title_font_size = 18;
    label_font_size = 18;
    axes_font_size = 14;
    contour_font_size = 14;
else
    % Figure font sizes
    title_font_size = 24;
    label_font_size = 24;
    axes_font_size = 20;
    contour_font_size = 20;
end

% Lay out title and filename formats
% Input filenames, plot titles, output filenames, datetime strings
% Compressed = for small plots
output_format_plot = '%s/%s_%s_plot_d%02.f_%s_%s_t%s_%s_c%dkm_%s_%s_%s.png'; % output path, data_tag,"dif",domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,dist_label,timestamp_file
            % for obs, domain  = d00, mem = 000; bounding shapes: bx, el, nb; 
            % if not using point radius and req, n0k_0; demons flag = d, nd (timestamp should be that of unmodified version)
input_format_stats = '%s/%s_%s_d%02.f_%s_%s_t%s_%s_c%dkm_%s_%s_%s.mat'; % output path, data source, plot type, domain, member #, data type, threshold, bounding shape, cutoff dist, testing values string, datetime, demons flag
input_format_stats_local = '%s_%s_d%02.f_%s_%s_t%s_%s_c%dkm_%s_%s_%s.mat'; % data source, plot type, domain, member #, data type, threshold, bounding shape, cutoff dist, testing values string, datetime, demons flag
datetime_format_file = '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]
output_format_table = '%s/%s_%s_table_d%02.f_%s_%s_t%s_%s_c%dkm_%s_%s_%s.txt'; % output path, data_tag,dif/stats,domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,dist_label,timestamp_file
output_format_mat = '%s/%s_%s_mat_d%02.f_%s_%s_t%s_%s_c%dkm_%s_%s_%s.mat'; % output path, data_tag,dif/stats,domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,dist_label,timestamp_file
%writetable(error_q_table,sprintf(output_table_format,output_path,data_tag,"dif",domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,timestamp_file),'Delimiter','\t','WriteRowNames',false);
        

% Reference filenames: wrf_enkf_output_d02_001_2020-02-07_14-00-00.nc

% Filepaths
if(local)
    mode = 'local';
    input_path_base = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2020';
    input_path_obs = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    intermediate_path = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code/output/tracking';
    path_to_code = "C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Code";
    path_to_extra_code = './downloaded_code';
    input_path_ecmwf = 'C:/Users/JonSe/Documents/Actual Documents/PSU/IMPACTS/Data/ECMWF';
    % Experiment paths
    exp_path_1 = 'AIR';
    exp_path_2 = 'CONV/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
else
    mode = 'server';
    input_path_base = '/storage/home/jjs5895/projects/IMPACTS/data/2020/';
    input_path_obs = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/GIS';
    intermediate_path = '/storage/home/jjs5895/projects/IMPACTS/intermediate';
    output_path_base = '/storage/home/jjs5895/projects/IMPACTS/output/tracking';
    path_to_code = "/storage/home/jjs5895/projects/IMPACTS/code";
    path_to_extra_code = '/storage/home/jjs5895/projects/IMPACTS/code/downloaded_code'; % Specify path to borders.m
    input_path_ecmwf = '/storage/home/jjs5895/projects/IMPACTS/data/obs_2020/ECMWF';
    % Experiment paths
    exp_path_1 = 'AIRCFT/fc';
    exp_path_2 = 'CONV/fc/';
    exp_path_3 = 'AIRCFT-F-18Z';
    exp_path_4 = 'CONV-F-18Z';
    exp_path_5 = 'NODA-14Z/';
end

% Experiment names
exp_name_1 = 'AIR';
exp_name_2 = 'CONV';
exp_name_3 = 'AIR-F';
exp_name_4 = 'CONV-F';
exp_name_5 = 'NODA';

%% 0D. Settings- Data (WRF 2020 Case)

% Basic details
domain = 2;
data_type = "refl"; % Variable being analyzed
units = "dBZ";

% NetCDF retrieval
lon_name = "XLONG"; 
lat_name = "XLAT"; 
data_name = "REFL_10CM"; % Variable being analyzed

% Note whether transposition or compositing is necessary to properly align
% the data with desired grid
transpose_data = true;
flip_data = true;
z_composite_data = true;

% Hour increment size
hour_step = 1; 

% Figure specs
fig_x = 50;
fig_y = 50;
fig_width = 700;
fig_height = 850;
fig_width_large = 1500;
fig_height_large = 850;
fig_width_summary = 1200;
fig_height_summary = 900;

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

cmap = 'reflmap'; % Choice of colorbar colormap
clim_lower = -30; % Colorbar limits
clim_upper = 75;

% Snowband determination vars
snowband_length_req = 250000; % 250 km
snowband_wl_ratio = 1/3;
snowband_portion_cold = 0.5; % Proportion of snowband that has to be in "cold" zone (below freezing)

% Display options  
plot_major_axis = true;
plot_centroids = true;
centroid_size = 60;
extend_major_axis = true;
label_clusters = true;
weight_centroids = true; % True = centroids weighted by dBZ % of object total

% System version controls
use_v6 = false; % False = use v5 SLINK with neighborhood cutoff; True = V6 with distance transform cutoff

edge_buffer = 42;
nodata_val = NaN;

clean_demons = true;

% Only one of these, at most, should ever be true
regrid_obs_to_wrf = 1;
regrid_wrf_to_obs = 0;

shape_string = 'nb';

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
%if(lon_gis_trim(y,x) < ((lat_gis_trim(y,x)-b)/m))

%% 0E. Obs Settings 
data_src_obs = "nex";
file_format_obs = '%s/n0q_%04.f%02.f%02.f%02.f%02.f.png';

lat_obs_raw = 49.9975:-0.005:23.0025;
lon_obs_raw = -125.9975:0.005:-65.0025;
[lon_obs,lat_obs] = meshgrid(lon_obs_raw,lat_obs_raw);

%% 1-B. Setup

if(use_weighted_distances)
    weight_flag = "ON";
    weight_string = "w1"; % Weighted
    weight_string_plot = "Weighted";
else
    weight_flag = "OFF";
    weight_string = "w0"; % Unweighted
    weight_string_plot = "Unweighted";
end

if(use_transformed_distances)
    trans_flag = "ON";
    trans_string = "dm1"; % Weighted
else
    trans_flag = "OFF";
    trans_string = "dm0"; % Unweighted
end

dist_label = sprintf('%s_%s',weight_string,trans_string);

% Announce operating mode
fprintf('Starting %s in %s mode.\n',upper(script_name),upper(mode));
fprintf('Weighted Distances: %s\n',weight_flag);
fprintf('Transformed Distances: %s\n',trans_flag);
tic;


if(~show_plots)
    set(groot,'DefaultFigureVisible','off') % Turn off figure popups for local
end

switch exp_choice
    case 0
        input_path = input_path_obs;
        exp_name = "OBS";
        num_members = 1;
    case 1
        input_path = sprintf('%s/%s',input_path_base,exp_path_1);
        exp_name = exp_name_1;
    case 2
        input_path = sprintf('%s/%s',input_path_base,exp_path_2);
        exp_name = exp_name_2;
    case 3
        input_path = sprintf('%s/%s',input_path_base,exp_path_3);
        exp_name = exp_name_3;
    case 4
        input_path = sprintf('%s/%s',input_path_base,exp_path_4);
        exp_name = exp_name_4;
    case 5
        input_path = sprintf('%s/%s',input_path_base,exp_path_5);
        exp_name = exp_name_5;
    otherwise
        error('Unrecognized experiment choice');
end
data_src = lower(exp_name);


if(regrid_obs_to_wrf)
    data_tag = sprintf('%s_%s',"OBS",exp_name);
else
    data_tag = sprintf('%s_%s',exp_name,"OBS");
end

if(compare_against == 0)
    data_tag_nd = sprintf('%s_%s',"OBS","OBS");
else
    data_tag_nd = data_tag;
end

if(use_base_refl)
    vcb_string = "ba";
else
    vcb_string = "vc";
end


%% 2. Initial Processing

% Directory structure:
% output/np/(run_name)/(variable_name)/(plot_size)
output_path = sprintf('%s/%s/%s',output_path_base,run_name,data_type);

addpath(path_to_extra_code);

% Time counters [NOTE: NOT ROBUST TO CROSSING MONTH BOUNDARIES]
max_day = [31, 28, 31, 30, 31, 30, 31, 31, 30, 30, 31, 31]; % Number of days in each month
time_count = (end_day - start_day + 1)*24 - (start_hour) - (24 - end_hour) + 1;

% If output folders do not exist, create them
if(~isfolder(output_path))
    mkdir(output_path);
end

load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'n_idx','s_idx','e_idx','w_idx','lon_wrf_trim','lat_wrf_trim');
lon_wrf = lon_wrf_trim;
lat_wrf = lat_wrf_trim;
ylen_wrf = size(lon_wrf,1);
xlen_wrf = size(lon_wrf,2);

if(use_member_list)
    member_list = custom_member_list;
    num_members = length(custom_member_list);
else
    member_list = 1:num_members;
end

num_directions = 2; % A->B and B->A
num_thresholds = length(thresh_list);
num_stats = 8; % Min, Q25, Median, Q75, Max, Mean, Sum (for weighted), Count

% Set up stats stats container (stats of the stats)
box_stats = zeros(num_members,num_directions,50);
big_box_stats = zeros(time_count,num_members,num_directions,num_thresholds,num_stats);
big_box_stats(:,:,:,:,:) = NaN;
%big_box_stats_weighted = zeros(time_count,num_members,num_directions,num_thresholds,num_stats);

centroid_count_collection = zeros(time_count,num_directions,num_thresholds);

if(num_members > 20)
    fig_width = fig_width_large;
end

%% 3. MAIN LOOP

thresh_idx = 1;

for band_threshold = thresh_list
for nbh_radius = radius_list
for nbh_req = req_list
for distance_cutoff = distance_list

% Specify threshold strings
if(use_sd_threshold)
    thresh_string = sprintf('sd%01d-%02.f',floor(sd_mult),100*mod(sd_mult,1));
else
    thresh_string = sprintf('bt%02.f',band_threshold);
end


testing_value_string = sprintf('%s_c%dkm_n%dkm_nr%01.f_pr%01.f',vcb_string,distance_cutoff/1000,nbh_radius/1000,nbh_req,point_req);
fprintf('Initializing: T %ddBZ, C %dkm, NR %dkm, NPR %d, MinPR %d\n',band_threshold,distance_cutoff/1000,nbh_radius/1000,nbh_req,point_req);
    

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
    box_stats(:,:,:) = NaN;

    % For each member:
    for member_idx = member_list

        member_string = sprintf('%03.f',member_idx);

        if(show_progress)
            fprintf('Processing: %s Member %s | %s\n',exp_name,member_string,timestamp_file);
        end
        
        %% Load

        if(compare_against == 0)
            member_string_a = "001";
        else
            member_string_a = member_string;
        end
        
        % Non-Demons file
        demons_string = 'nd';
        stats_filename = sprintf(input_format_stats,intermediate_path,data_tag_nd,'stats',domain,member_string_a,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string);
        stats_filename_local = sprintf(input_format_stats_local,data_tag_nd,'stats',domain,member_string_a,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string);
        
        %air_stats_d02_001_refl_tsd1-25_nb_c10k_n10k_26.0_202002071500_d
        %air_stats_d02_001_refl_tsd1-25_nb_c10k_n10k_26.0_202002071500_nd
        %snowband_tracking_stats_air_1400_d
        %stats_filename = 'snowband_tracking_stats_air_1400_d.mat';

        if(show_loading_msg)
            fprintf('Attempting to load: %s\n',stats_filename_local);
        end
        load(stats_filename,'pt_clusters','cluster_ids','centroids_trim','data_main','D_save');
            
        pt_clusters_a = pt_clusters;
        cluster_ids_a = cluster_ids;
        centroids_a = centroids_trim;
        data_a = data_main;
        D_a = D_save;

        % Compute cluster weights
        % NOTE: THIS ASSUMES THE CENTROIDS (CHOSEN CLUSTERS W/ POINTS >=
        % POINT_REQ) ARE STILL IN ASCENDING ID ORDER
        num_centroids_a = size(centroids_a,1);
        num_clusters_a = numel(cluster_ids_a);
        num_pts_a = numel(pt_clusters_a);

        if(show_debug)
            fprintf('A Centroids: %d\n',num_centroids_a);
            fprintf('A Running Total: %d\n',centroid_count_collection(time_idx,1,thresh_idx)+num_centroids_a);
        end
        centroid_count_collection(time_idx,1,thresh_idx) = num_centroids_a + centroid_count_collection(time_idx,1,thresh_idx);

        weights_a = zeros(num_centroids_a,1);
        use_cluster_a = ones(num_clusters_a,1);

        centroid_idx = 1;
        % Remove small clusters (and ocean clusters) from the count
        for cluster_idx = 1:num_clusters_a
            cluster_id_num = cluster_ids_a(cluster_idx);
            cluster_members = find(pt_clusters_a == cluster_id_num);
            if(numel(cluster_members) < point_req)
                num_pts_a = num_pts_a - numel(cluster_members);
                use_cluster_a(cluster_idx) = 0;
            elseif((apply_ocean_mask && (centroids_a(centroid_idx,1) >= ((centroids_a(centroid_idx,2)-b)/m))))
                centroid_idx = centroid_idx + 1;
                num_pts_a = num_pts_a - numel(cluster_members);
                use_cluster_a(cluster_idx) = 0;
            end
        end

        % Remove ocean objects from the count and centroids list (apply mask)
        if(apply_ocean_mask)
            use_centroids = ones(num_centroids_a,1);
            new_num_centroids = num_centroids_a;
            for centroid_idx = 1:num_centroids_a
                %if(lon_gis_trim(y,x) >= ((lat_gis_trim(y,x)-b)/m))
                if(centroids_a(centroid_idx,1) >= ((centroids_a(centroid_idx,2)-b)/m))
                    new_num_centroids = new_num_centroids - 1;
                    use_centroids(centroid_idx) = 0;
                end
            end
            % Adjust centroid list (WIP)
            centroids_a = centroids_a(use_centroids==1,:);
            num_centroids_a = new_num_centroids;
        end


        % Compute weights as proportion of remaining points
        centroid_idx = 1;
        for cluster_idx = 1:num_clusters_a
            cluster_id_num = cluster_ids_a(cluster_idx);
            cluster_members = find(pt_clusters_a == cluster_id_num);
            if(use_cluster_a(cluster_idx))
                weights_a(centroid_idx) = numel(cluster_members)/num_pts_a;
                centroid_idx = centroid_idx + 1;
            end
        end
        
        %% File B: Experiment (Demons-adjusted, or backtrack to pre-demons)
        
        demons_string = 'd';
        stats_filename = sprintf(input_format_stats,intermediate_path,data_tag,'stats',domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string);
        stats_filename_local = sprintf(input_format_stats_local,data_tag,'stats',domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string);

        if(show_loading_msg)
            fprintf('Attempting to load: %s\n',stats_filename_local);
        end
        load(stats_filename,'pt_clusters','cluster_ids','centroids_trim','data_main','D_save');
        
        pt_clusters_b_d = pt_clusters;
        cluster_ids_b_d = cluster_ids;
        centroids_b_d = centroids_trim;
        data_b_d = data_main;
        D_b = D_save;

        demons_string = 'nd';
        stats_filename = sprintf(input_format_stats,intermediate_path,data_tag,'stats',domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string);
        stats_filename_local = sprintf(input_format_stats_local,data_tag,'stats',domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string);

        if(show_loading_msg)
            fprintf('Attempting to load: %s\n',stats_filename_local);
        end
        load(stats_filename,'pt_clusters','cluster_ids','centroids_trim','data_main','D_save');
        pt_clusters_b_nd = pt_clusters;
        cluster_ids_b_nd = cluster_ids;
        centroids_b_nd = centroids_trim;
        data_b_nd = data_main;
    
        % Compute cluster weights
        % NOTE: THIS ASSUMES THE CENTROIDS (CHOSEN CLUSTERS W/ POINTS >=
        % POINT_REQ) ARE STILL IN ASCENDING ID ORDER
        num_centroids_b_d = size(centroids_b_d,1);
        num_clusters_b_d = numel(cluster_ids_b_d);
        num_pts_b_d = numel(pt_clusters_b_d);

        if(show_debug)
            fprintf('B Centroids: %d\n',num_centroids_b_d);
            fprintf('B Running Total: %d\n',centroid_count_collection(time_idx,2,thresh_idx)+num_centroids_b_d);
        end

        

        if(num_centroids_a == 0 || num_centroids_b_d == 0) % Skip if no data, to avoid crashes
            %;
        else
    
            weights_b_d = zeros(num_centroids_b_d,1);
            use_cluster_b_d = ones(num_clusters_b_d,1);
    
            centroid_idx = 1;
            % Remove small clusters (and ocean clusters) from the count
            for cluster_idx = 1:num_clusters_b_d
                cluster_id_num = cluster_ids_b_d(cluster_idx);
                cluster_members = find(pt_clusters_b_d == cluster_id_num);
                if(numel(cluster_members) < point_req)
                    num_pts_b_d = num_pts_b_d - numel(cluster_members);
                    use_cluster_b_d(cluster_idx) = 0;
                elseif((apply_ocean_mask && (centroids_b_d(centroid_idx,1) >= ((centroids_b_d(centroid_idx,2)-b)/m))))
                    centroid_idx = centroid_idx + 1;
                    num_pts_b_d = num_pts_b_d - numel(cluster_members);
                    use_cluster_b_d(cluster_idx) = 0;
                end
            end
    
    
            % Remove ocean objects from the count and centroids list (apply mask)
            if(apply_ocean_mask)
                use_centroids = ones(num_centroids_b_d,1);
                new_num_centroids = num_centroids_b_d;
                for centroid_idx = 1:num_centroids_b_d
                    %if(lon_gis_trim(y,x) >= ((lat_gis_trim(y,x)-b)/m))
                    if(centroids_b_d(centroid_idx,1) >= ((centroids_b_d(centroid_idx,2)-b)/m))
                        new_num_centroids = new_num_centroids - 1;
                        use_centroids(centroid_idx) = 0;
                    end
                end
                % Adjust centroid list
                centroids_b_d = centroids_b_d(use_centroids==1,:);
                num_centroids_b_d = new_num_centroids;
            end
    
            % Compute weights as proportion of remaining points
            centroid_idx = 1;
            for cluster_idx = 1:num_clusters_b_d
                cluster_id_num = cluster_ids_b_d(cluster_idx);
                cluster_members = find(pt_clusters_b_d == cluster_id_num);
                if(use_cluster_b_d(cluster_idx))
                    weights_b_d(centroid_idx) = numel(cluster_members)/num_pts_b_d;
                    centroid_idx = centroid_idx + 1;
                end
            end
        end

        if(use_transformed_distances)
            centroid_count_collection(time_idx,2,thresh_idx) = num_centroids_b_d + centroid_count_collection(time_idx,2,thresh_idx);
        end

        %% Compute distances

        output_name_save = sprintf(output_format_mat,intermediate_path,data_tag,"dif",domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,dist_label,timestamp_file);
        max_len = max(num_centroids_a,num_centroids_b_d);

        distances = zeros(num_centroids_a,num_centroids_b_d);
        ab_distances = zeros(max_len,1);
        ab_distances(:) = NaN;
        ab_links = zeros(max_len,1);
        ba_distances  = zeros(max_len,1);
        ba_distances(:) = NaN;
        ba_links = zeros(max_len,1);

        if(num_centroids_a == 0 || num_centroids_b_d == 0)
            num_centroids_a = 0;
            num_centroids_b_d = 0;
            fprintf('NO DATA. SKIPPING.\n')
            save(output_name_save,'num_centroids_a','ab_links','ab_distances','num_centroids_b_d','ba_links','ba_distances');
            continue;
        end

        for centroid_idx_a = 1:num_centroids_a
            for centroid_idx_b_d = 1:num_centroids_b_d
                dist = haversine_distance(centroids_a(centroid_idx_a,1),centroids_a(centroid_idx_a,2),centroids_b_d(centroid_idx_b_d,1),centroids_b_d(centroid_idx_b_d,2));
                distances(centroid_idx_a,centroid_idx_b_d) = round(dist,2);
            end
        end
        
        for centroid_idx_a = 1:num_centroids_a
            if(use_weighted_distances)
                dist_list = weights_a(centroid_idx_a).*distances(centroid_idx_a,:);
                [min_val,min_idx] = min(dist_list);
            else
                [min_val,min_idx] = min(distances(centroid_idx_a,:));
            end
            ab_distances(centroid_idx_a) = min_val;
            ab_links(centroid_idx_a) = min_idx;
        end
        
        for centroid_idx_b_d = 1:num_centroids_b_d
            if(use_weighted_distances)
                dist_list = weights_b_d(centroid_idx_b_d).*distances(:,centroid_idx_b_d);
                [min_val,min_idx] = min(dist_list);
            else
                [min_val,min_idx] = min(distances(:,centroid_idx_b_d));
            end
            ba_distances(centroid_idx_b_d) = min_val;
            ba_links(centroid_idx_b_d) = min_idx;
        end

        distances_d = distances;
        ab_distances_d = ab_distances;
        ab_links_d = ab_links;
        ba_distances_d  = ba_distances;
        ba_links_d = ba_links;

        %% Now backtrack with displacement matrix : first compute all initial stats again for ND version
        
        if(~use_transformed_distances)
            pt_clusters_b = pt_clusters_b_nd;
            cluster_ids_b = cluster_ids_b_nd;
            centroids_b = centroids_b_nd;
            data_b = data_b_nd;
        
            % Compute cluster weights
            % NOTE: THIS ASSUMES THE CENTROIDS (CHOSEN CLUSTERS W/ POINTS >=
            % POINT_REQ) ARE STILL IN ASCENDING ID ORDER
            num_centroids_b_nd = size(centroids_b_nd,1);
            num_clusters_b_nd = numel(cluster_ids_b_nd);
            num_pts_b_nd = numel(pt_clusters_b_nd);
    
            if(show_debug)
                fprintf('B Centroids: %d\n',num_centroids_b_nd);
                fprintf('B Running Total: %d\n',centroid_count_collection(time_idx,2,thresh_idx)+num_centroids_b_nd);
            end
    
            if(num_centroids_a == 0 || num_centroids_b_nd == 0) % Skip if no data, to avoid crashes
                %;
            else
        
                weights_b_nd = zeros(num_centroids_b_nd,1);
                use_cluster_b_nd = ones(num_clusters_b_nd,1);
        
                centroid_idx = 1;
                % Remove small clusters (and ocean clusters) from the count
                for cluster_idx = 1:num_clusters_b_nd
                    cluster_id_num = cluster_ids_b_nd(cluster_idx);
                    cluster_members = find(pt_clusters_b_nd == cluster_id_num);
                    if(numel(cluster_members) < point_req)
                        num_pts_b_nd = num_pts_b_nd - numel(cluster_members);
                        use_cluster_b_nd(cluster_idx) = 0;
                    elseif((apply_ocean_mask && (centroids_b_nd(centroid_idx,1) >= ((centroids_b_nd(centroid_idx,2)-b)/m))))
                        centroid_idx = centroid_idx + 1;
                        num_pts_b_nd = num_pts_b_nd - numel(cluster_members);
                        use_cluster_b_nd(cluster_idx) = 0;
                    end
                end
        
        
                % Remove ocean objects from the count and centroids list (apply mask)
                if(apply_ocean_mask)
                    use_centroids = ones(num_centroids_b_nd,1);
                    new_num_centroids = num_centroids_b_nd;
                    for centroid_idx = 1:num_centroids_b_nd
                        %if(lon_gis_trim(y,x) >= ((lat_gis_trim(y,x)-b)/m))
                        if(centroids_b_nd(centroid_idx,1) >= ((centroids_b_nd(centroid_idx,2)-b)/m))
                            new_num_centroids = new_num_centroids - 1;
                            use_centroids(centroid_idx) = 0;
                        end
                    end
                    % Adjust centroid list
                    centroids_b_nd = centroids_b_nd(use_centroids==1,:);
                    num_centroids_b_nd = new_num_centroids;
                end
        
                % Compute weights as proportion of remaining points
                centroid_idx = 1;
                for cluster_idx = 1:num_clusters_b_nd
                    cluster_id_num = cluster_ids_b_nd(cluster_idx);
                    cluster_members = find(pt_clusters_b_nd == cluster_id_num);
                    if(use_cluster_b_nd(cluster_idx))
                        weights_b_nd(centroid_idx) = numel(cluster_members)/num_pts_b_nd;
                        centroid_idx = centroid_idx + 1;
                    end
                end
            end

            if(~use_transformed_distances)
                centroid_count_collection(time_idx,2,thresh_idx) = num_centroids_b_nd + centroid_count_collection(time_idx,2,thresh_idx);
            end
    
            %% Match newly trimmed ND centroids to D centroids via displacement matrix
            % This will have to be any amount of ND -> whichever D are closest
            % to interpolated output point
    
            max_len = max(num_centroids_a,num_centroids_b_nd);
            distances = zeros(num_centroids_a,num_centroids_b_nd);
            ab_distances = zeros(max_len,1);
            ab_distances(:) = NaN;
            ab_links = zeros(max_len,1);
            ba_distances  = zeros(max_len,1);
            ba_distances(:) = NaN;
            ba_links = zeros(max_len,1);

            if(num_centroids_a == 0 || num_centroids_b_nd == 0)
                num_centroids_a = 0;
                num_centroids_b_d = 0;
                fprintf('NO DATA. SKIPPING.\n')
                save(output_name_save,'num_centroids_a','ab_links','ab_distances','num_centroids_b_d','ba_links','ba_distances');
                continue;
            end
    
            %D_b; % (:,:,x-displacements),(:,:,y-displacements)
            % D_b is the same size as lon_wrf and lat_wrf, as it should be.
            % Therefore- 
            % 1) find closest grid point to each centroid's listed lat/lon, 
            % 2) take displacement, 
            % 3) find new closest D centroid to that location, 
            % 4) steal that linkage, compute new linkage distance
    
            for centroid_idx = 1:num_centroids_b_nd
                clon = centroids_b_nd(centroid_idx,1);
                clat = centroids_b_nd(centroid_idx,2);
                [y_idx,x_idx] = find_closest(lat_wrf,lon_wrf,clat,clon);
                y_dis = D_b(y_idx,x_idx,2);
                x_dis = D_b(y_idx,x_idx,1);
                y_new = round(y_idx + y_dis);
                x_new = round(x_idx + x_dis);
                y_new = max(min(y_new,ylen_wrf),1);
                x_new = max(min(x_new,xlen_wrf),1);
    
                min_dist = 999999999;
                closest_centroid = 0;
                for sub_centroid_idx = 1:num_centroids_b_d
                    dist = haversine_distance(lon_wrf(y_new,x_new),lat_wrf(y_new,x_new),centroids_b_d(sub_centroid_idx,1),centroids_b_d(sub_centroid_idx,2));
                    if(dist < min_dist)
                        min_dist = dist;
                        closest_centroid = sub_centroid_idx;
                    end
                end
    
                matched_centroid_a = ba_links_d(closest_centroid);
    
                new_dist = haversine_distance(centroids_a(matched_centroid_a,1),centroids_a(matched_centroid_a,2),centroids_b_nd(centroid_idx,1),centroids_b_nd(centroid_idx,2));
               
                if(use_weighted_distances)
                    mult_a = weights_a(matched_centroid_a);
                    mult_b = weights_b_nd(centroid_idx);
                else
                    mult_a = 1; mult_b = 1;
                end

                ba_distances(centroid_idx) = new_dist*mult_b;
                ba_links(centroid_idx) = matched_centroid_a;
                if((ab_links(matched_centroid_a) == 0) || (new_dist*mult_a < ab_distances(matched_centroid_a)))
                    ab_distances(matched_centroid_a) = new_dist*mult_a;
                    ab_links(matched_centroid_a) = centroid_idx;
                end
    
            end

            % Compute displaced latlon grid so the OA can be calculated
            lon_wrf_dis = zeros(size(lon_wrf));
            lat_wrf_dis = zeros(size(lat_wrf));
            %lon_wrf_dis_idx = zeros(size(lon_wrf));
            %lat_wrf_dis_idx = zeros(size(lat_wrf));
            for y = 1:ylen_wrf
                for x = 1:xlen_wrf
                    y_dis = D_b(y,x,2);
                    x_dis = D_b(y,x,1);
                    y_new = max(min(round(y + y_dis),ylen_wrf),1);
                    x_new = max(min(round(x + x_dis),xlen_wrf),1);
                    lon_wrf_dis(y,x) = lon_wrf(y_new,x_new);
                    lat_wrf_dis(y,x) = lat_wrf(y_new,x_new);
                    %lon_wrf_dis_idx(y,x) = x_new;
                    %lat_wrf_dis_idx(y,x) = y_new;
                end
            end

            % For each Observations centroid
            % 1) take the D centroid link it had, 
            % 2) get its coordinates, 
            % 3) then reference the fucked up grid computed above to see what 
            % gridpoint was displaced into (or closest to) that one
            % 4) find ND centroid closest to that new gridpoint, set as
            % link and take distance
            for centroid_idx = 1:num_centroids_a
                linked_b_d = ab_links_d(centroid_idx);
                clon = centroids_b_d(linked_b_d,1);
                clat = centroids_b_d(linked_b_d,2);
                [y_idx,x_idx] = find_closest(lat_wrf_dis,lon_wrf_dis,clat,clon);
                og_lon = lon_wrf(y_idx,x_idx);
                og_lat = lat_wrf(y_idx,x_idx);
    
                min_dist = 999999999;
                closest_centroid = 0;
                for sub_centroid_idx = 1:num_centroids_b_nd
                    dist = haversine_distance(og_lon,og_lat,centroids_b_nd(sub_centroid_idx,1),centroids_b_nd(sub_centroid_idx,2));
                    if(dist < min_dist)
                        min_dist = dist;
                        closest_centroid = sub_centroid_idx;
                    end
                end
    
                matched_centroid_b = closest_centroid;
    
                if(use_weighted_distances)
                    mult_a = weights_a(matched_centroid_a);
                    mult_b = weights_b_nd(centroid_idx);
                else
                    mult_a = 1; mult_b = 1;
                end

                ab_distances(centroid_idx) = min_dist;
                ab_links(centroid_idx) = matched_centroid_b;
    
            end
    
            num_centroids_b = num_centroids_b_nd;
            pt_clusters_b = pt_clusters_b_nd;
            cluster_ids_b = cluster_ids_b_nd;
            centroids_b = centroids_b_nd;
            data_b = data_b_nd;
        end
        

%%


        if(~exist(output_name_save) || override_saved_stats)
            save(output_name_save,'num_centroids_a','ab_links','ab_distances','num_centroids_b','ba_links','ba_distances');
        end

        %% Print

        dist_strings_a = string(ab_distances);
        dist_strings_b = string(ba_distances);
        link_strings_a = string(ab_links);
        link_strings_b = string(ba_links);

        %max_len = num_centroids_a;
        if(num_centroids_a > num_centroids_b)
            dist_strings_b(num_centroids_b+1:num_centroids_a) = "N/A";
            link_strings_b(num_centroids_b+1:num_centroids_a) = "N/A";
        elseif(num_centroids_b > num_centroids_a)
            dist_strings_a(num_centroids_a+1:num_centroids_b) = "N/A";
            link_strings_a(num_centroids_a+1:num_centroids_b) = "N/A";
            %max_len = num_centroids_b;
        end

        link_strings_a(max_len+1) = "Mean";
        dist_strings_a(max_len+1) = string(mean(ab_distances(1:num_centroids_a)));
        link_strings_b(max_len+1) = "N/A";
        dist_strings_b(max_len+1) = string(mean(ba_distances(1:num_centroids_b)));

        link_strings_a(max_len+2) = "Median";
        dist_strings_a(max_len+2) = string(median(ab_distances(1:num_centroids_a)));
        link_strings_b(max_len+2) = "N/A";
        dist_strings_b(max_len+2) = string(median(ba_distances(1:num_centroids_b)));

        % If output folders do not exist, create them
        if(~isfolder(output_path))
            mkdir(output_path);
        end

        error_q_table = table(link_strings_a,dist_strings_a,link_strings_b,dist_strings_b,'VariableNames',["AB-Links(#)","AB-Dist(m)","BA-Links(#)","BA-Dist(m)"]);
        writetable(error_q_table,sprintf(output_format_table,output_path,data_tag,"dif",domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,dist_label,timestamp_file),'Delimiter','\t','WriteRowNames',false);

        %% Add to container
        % box_stats = zeros(num_members,num_directions,num_entries)
        num_entries = length(ab_distances);
        box_stats(member_idx,1,1:num_entries) = ab_distances;
        box_stats(member_idx,2,1:num_entries) = ba_distances;
        % big_box_stats = zeros(time_count,num_members,num_directions,num_thresholds,num_stats);
        big_box_stats(time_idx,member_idx,1,thresh_idx,1:num_stats) = [min(ab_distances),quantile(ab_distances,0.25),quantile(ab_distances,0.5),quantile(ab_distances,0.75),max(ab_distances),nanmean(ab_distances),nansum(ab_distances),numel(ab_distances(~isnan(ab_distances)))];
        big_box_stats(time_idx,member_idx,2,thresh_idx,1:num_stats) = [min(ba_distances),quantile(ba_distances,0.25),quantile(ba_distances,0.5),quantile(ba_distances,0.75),max(ba_distances),nanmean(ba_distances),nansum(ba_distances),numel(ba_distances(~isnan(ba_distances)))];
    end
        
    %% Boxplot

    % Make box & whisker plot of distances for each ensemble member in this
    % variation
    if(make_boxplot)
        box_stats_km = box_stats./1000; % m -> km
        boxlabels = 1:num_members;

        box_stats_plot = box_stats_km(:,1,:);
        box_stats_plot = squeeze(box_stats_plot);
    
        f = figure('Position',[fig_x fig_y fig_width fig_height]);
        subplot(2,1,1);
        h = boxplot(box_stats_plot',boxlabels);
        hold on;
        grid minor;
        title(sprintf("[OBS->%s] Matched Obj. Dist. Stats (T%s)",exp_name,string(band_threshold)),'FontSize',title_font_size);
        xlabel("Ensemble Member",'FontSize',label_font_size);
        ylabel("Distance (km)",'FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);
        yl = ylim;
        ylim([-1*y_neg_scaling*yl(2) yl(2)]);
        
        subplot(2,1,2)
        box_stats_plot = box_stats_km(:,2,:);
        box_stats_plot = squeeze(box_stats_plot);
        %h = boxplot(box_stats_plot',boxlabels,'BoxStyle','filled','MedianStyle','target','OutlierSize',4,'Symbol','o','Jitter',0.5);
        h = boxplot(box_stats_plot',boxlabels,'PlotStyle','compact','LabelOrientation','horizontal');
        grid minor;
        title(sprintf("[%s->OBS] Matched Obj. Dist. Stats (T%s)",exp_name,string(band_threshold)),'FontSize',title_font_size);
        xlabel("Ensemble Member",'FontSize',label_font_size);
        ylabel("Distance (km)",'FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);
        yl = ylim;
        ylim([-1*y_neg_scaling*yl(2) yl(2)]);

        %set(gca,'xtick',[]);
    
        saveas(gcf,sprintf(output_format_plot,output_path,data_tag,"boxwhisk",domain,"000",data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,dist_label,timestamp_file));
        close("all");
    end

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
end
end
end
thresh_idx = thresh_idx+1;
end
%%
% Make a summary boxplot for the entire ensemble as one box, showing all
% timestamps and thresholds on one plot
output_name_plot = sprintf(output_format_plot,output_path,data_tag,"boxwhisk",domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,dist_label,timestamp_file);
if(make_condensed_boxplots && (remake_plots || ~exist(output_name_plot))) %Only bother if flag is set, and it either doesn't exist or should be overridden

    centroid_total_a = centroid_count_collection(1,1,1); % Sample
    centroid_total_b = centroid_count_collection(1,2,1);
    distance_collection_a = zeros(time_count*num_thresholds,centroid_total_a).*NaN;
    distance_collection_b = zeros(time_count*num_thresholds,centroid_total_b).*NaN;

    boxlabels = string(zeros(time_count*num_thresholds,1));

    year = start_year;
    month = start_month;
    day = start_day;
    hour = start_hour;

    % Collect total distance summary variable
    for time_idx = 1:time_count
        timestamp_file = sprintf(datetime_format_file,year,month,day,hour); 
        for thresh_idx = 1:num_thresholds
            permutation_idx = (time_idx-1)*num_thresholds + thresh_idx;

            band_threshold = thresh_list(thresh_idx);
            if(use_sd_threshold)
                thresh_string = sprintf('sd%01d-%02.f',floor(sd_mult),100*mod(sd_mult,1));
            else
                thresh_string = sprintf('bt%02.f',band_threshold);
            end

            boxlabels(permutation_idx) = sprintf("T%d|%02.f00",band_threshold,hour);

            centroid_total_a = centroid_count_collection(time_idx,1,thresh_idx);
            centroid_total_b = centroid_count_collection(time_idx,2,thresh_idx);
            member_collection_a = zeros(centroid_total_a,1).*NaN;
            member_collection_b = zeros(centroid_total_b,1).*NaN;
            a_idx = 1;
            b_idx = 1;
            for member_sub_idx = 1:num_members
                member_string = sprintf('%03.f',member_sub_idx);
                output_name_save = sprintf(output_format_mat,intermediate_path,data_tag,"dif",domain,member_string,data_type,thresh_string,shape_string,distance_cutoff/1000,testing_value_string,dist_label,timestamp_file);
                load(output_name_save,'num_centroids_a','ab_distances','num_centroids_b','ba_distances');
                member_collection_a(a_idx:a_idx+num_centroids_a-1) = ab_distances(1:num_centroids_a);
                member_collection_b(b_idx:b_idx+num_centroids_b-1) = ba_distances(1:num_centroids_b);
                a_idx = a_idx + num_centroids_a;
                b_idx = b_idx + num_centroids_b;
            end
            distance_collection_a(permutation_idx,1:centroid_total_a) = member_collection_a;
            distance_collection_b(permutation_idx,1:centroid_total_b) = member_collection_b;
        end
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

    dist_plot_a = distance_collection_a./1000; % m -> km
    dist_plot_b = distance_collection_b./1000;

    f = figure('Position',[fig_x fig_y fig_width_summary fig_height_summary]);
    subplot(2,1,1);
    h = boxplot(dist_plot_a',boxlabels');
    hold on;
    grid minor;
    title(sprintf("[OBS->%s] All-Ensemble Distances: Threshold Sensitivity Test (%s)",exp_name,weight_string_plot),'FontSize',title_font_size);
    xlabel("Time (Hour UTC) & Threshold (dBZ)",'FontSize',label_font_size);
    if(use_weighted_distances)
        ylabel("Distance (km*pt)",'FontSize',label_font_size)
    else
        ylabel("Distance (km)",'FontSize',label_font_size);
    end
    set(gca,'Fontsize',axes_font_size);
    yl = ylim;
    ylim([-1*y_neg_scaling*yl(2) yl(2)]);

    subplot(2,1,2)
    h = boxplot(dist_plot_b',boxlabels');
    grid minor;
    title(sprintf("[%s->OBS] All-Ensemble Distances: Threshold Sensitivity Test (%s)",exp_name,weight_string_plot),'FontSize',title_font_size);
    xlabel("Time (Hour UTC) & Threshold (dBZ)",'FontSize',label_font_size);
    if(use_weighted_distances)
        ylabel("Distance (km*pt)",'FontSize',label_font_size)
    else
        ylabel("Distance (km)",'FontSize',label_font_size);
    end
    set(gca,'Fontsize',axes_font_size);
    yl = ylim;
    ylim([-1*y_neg_scaling*yl(2) yl(2)]);

    saveas(gcf,output_name_plot);
end

%% Save big stats table


%big_box_stats = zeros(time_count,num_members,num_directions,num_thresholds,num_stats);
timestamp = string(start_hour:end_hour);
member_num = string(1:num_members);
direction = ["OA","AO"];
threshold = string(thresh_list);
stats_strings = ["Min","Q25","Median","Q75","Max","Mean","Sum","Count"];

%T1 = [repelem(time_strings,numel(member_strings))' repmat(member_strings',numel(time_strings),1)];
%T2 = [repelem(T1,numel(direction_strings))' repmat(direction_strings',size(T2,1),1)];
T = combinations(timestamp,member_num,direction,threshold);
num_rows = size(T,1);
T.Min = zeros(num_rows,1);
T.Q25 = zeros(num_rows,1);
T.Median = zeros(num_rows,1);
T.Q75 = zeros(num_rows,1);
T.Max = zeros(num_rows,1);
T.Mean = zeros(num_rows,1);
T.Sum = zeros(num_rows,1);
T.Count = zeros(num_rows,1);

% big_box_stats = zeros(time_count,num_members,num_directions,num_thresholds,num_stats);
row_idx = 1;
for time_idx = 1:time_count
    for member_idx = member_list
        for direction_idx = 1:2
            for threshold_idx = 1:num_thresholds
                T(row_idx,5:(5+num_stats-1)) = num2cell(squeeze(big_box_stats(time_idx,member_idx,direction_idx,threshold_idx,:))');
                row_idx = row_idx+1;
            end
        end
    end
end

writetable(T,sprintf(output_format_table,output_path,data_tag,"stats",domain,"XXX",data_type,"btXX",shape_string,distance_cutoff/1000,testing_value_string,dist_label,"20200207XXXX"),'Delimiter','\t','WriteRowNames',false);
        
%%
fprintf('%s run complete.\n',script_name);
toc
fprintf('DONE.\n')

%% Local functions

% Find indices of closest point in latlon grids
function [y_out,x_out] = find_closest(lat_grid,lon_grid,target_lat,target_lon)
    y_out = NaN; x_out = NaN;
    min_dist = 999999999;
    for y = 1:size(lat_grid,1)
        prev_dist = 999999999;
        for x = 1:size(lat_grid,2)
            dist = haversine_distance(target_lon,target_lat,lon_grid(y,x),lat_grid(y,x));
            if(dist < min_dist)
                min_dist = dist;
                y_out = y;
                x_out = x;
            elseif(dist > prev_dist)
                break;
            end
            prev_dist = dist;
        end
    end
end