% condense_snowband_tables.m
%
% Purpose: 
% 
% Author(s): Jon Seibert
% Last updated: 22 Jan 2025
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


script_name = 'condense_snowband_tables.m';

%% 0A. Script Controls

% Denotes output folder name
run_name = 'thresh_dif_5_pred';

exp_choice = 1; %  
compare_against = 0; % obs

use_weighted_distances = 0; % Weight by num points in each object (1) or by number of objects (0)
use_transformed_distances = 0; % Use distances from WRF Demons Transformation -> OBS (1), or backtrack to retrieve distances from WRF->OBS w/o Demons, using output D matrix (0)
                    % OPTION 0 NOT YET IMPLEMENTED
use_base_refl = 1; % Whether to use vertical composite (0) or simulated base (1) reflectivity

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
custom_member_list = 1:3;

% Threshold settings
object_threshold = 0; %dBZ; minimum value to be considered a "reflectivity object"
%band_threshold = 40; % dBZ; minimum value to be considered a "snowband" (Overridden)
%band_threshold_b = 35; % dBZ threshold for Demons target
band_threshold_b_adjustment = 0; % modifier for threshold for demons target (obs)
sd_mult = 1.25; % Standard deviation value to take for embedded inner snowband
use_sd_threshold = 0; % Flag: Override with adaptive threshold = mean + sd_mult*SD?

% SLINK Parameters
distance_cutoff = 10000; % (meters) Maximum distance allowed to be bridged when adding a point to a cluster
min_clusters = 1; % Minimum number of clusters allowed
point_req = 20; % Minimum number of points to consider a cluster for snowband status.
nbh_radius = 10000; % Radius  of neighborhood circle [m] (6 km; current clustering cutoff radius is 10 km... try that? link the two?)
nbh_req = 26; % Number of neighbors in radius required to be considered part of a shared band? % Will I be using this?
%dt_cutoff_value = 1.2; % Euclidean distance from edge threshold for merging SLINK clusters

% Sensistivity testing controls
thresh_list = [20,35,50];
%radius_list = [10000]; % List of nbh radii
%req_list = [20, 26, 32]; % List of nbh point count requirements [20 26 32]
%req_list = [26]; % List of nbh point count requirements [20 26 32]
%distance_list = [14000]; % List of distance cutoffs
%distance_list = [10000];

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

%% 0E. Obs Settings 
data_src_obs = "nex";
file_format_obs = '%s/n0q_%04.f%02.f%02.f%02.f%02.f.png';

lat_obs_raw = 49.9975:-0.005:23.0025;
lon_obs_raw = -125.9975:0.005:-65.0025;
[lon_obs,lat_obs] = meshgrid(lon_obs_raw,lat_obs_raw);

%% 1-B. Setup

switch exp_choice
    case 0
        input_path = input_path_obs;
        exp_name = "OBS";
        num_members = 1;
        use_member_list = 0;
        use_demons = 0;
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
tic;



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

%testing_value_string = sprintf('%s_n%dkm_%01.1f',vcb_string,nbh_radius/1000,nbh_req);
testing_value_string = sprintf('%s_c%dkm_n%dkm_nr%01.f_pr%01.f',vcb_string,distance_cutoff/1000,nbh_radius/1000,nbh_req,point_req);

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

if(use_member_list)
    member_list = custom_member_list;
    num_members = length(custom_member_list);
else
    member_list = 1:num_members;
end

num_directions = 2; % A->B and B->A
num_thresholds = length(thresh_list);
num_stats = 6; % Min, Q25, Median, Q75, Max, Mean


%% Import big stats tables

exp_name = "AIR";
weight_string = "w0";
data_tag = sprintf('%s_%s',"OBS",exp_name);
dist_label = sprintf('%s_%s',weight_string,trans_string);
table_name = sprintf(output_format_table,output_path,data_tag,"stats",domain,"XXX",data_type,"btXX",shape_string,distance_cutoff/1000,testing_value_string,dist_label,"20200207XXXX");
T_A = readtable(table_name);

exp_name = "CONV";
weight_string = "w0";
data_tag = sprintf('%s_%s',"OBS",exp_name);
dist_label = sprintf('%s_%s',weight_string,trans_string);
table_name = sprintf(output_format_table,output_path,data_tag,"stats",domain,"XXX",data_type,"btXX",shape_string,distance_cutoff/1000,testing_value_string,dist_label,"20200207XXXX");
T_B = readtable(table_name);

exp_name = "AIR";
weight_string = "w1";
data_tag = sprintf('%s_%s',"OBS",exp_name);
dist_label = sprintf('%s_%s',weight_string,trans_string);
table_name = sprintf(output_format_table,output_path,data_tag,"stats",domain,"XXX",data_type,"btXX",shape_string,distance_cutoff/1000,testing_value_string,dist_label,"20200207XXXX");
T_C = readtable(table_name);

exp_name = "CONV";
weight_string = "w1";
data_tag = sprintf('%s_%s',"OBS",exp_name);
dist_label = sprintf('%s_%s',weight_string,trans_string);
table_name = sprintf(output_format_table,output_path,data_tag,"stats",domain,"XXX",data_type,"btXX",shape_string,distance_cutoff/1000,testing_value_string,dist_label,"20200207XXXX");
T_D = readtable(table_name);

exp_name = "XXXX";

%% Extract values

% Variables:
% 1) AIR Ens Mean # Objects (#)
% 2) CONV Ens Mean # Objects (#)
% 3) AIR Ens Mean of Per-Member Mean Match Distance (m?km? CHECK)
% 4) CONV Ens Mean of Per-Member Mean Match Distance (m)
% 5) Ens Mean Difference in Per-Member Mean Match Distances (
% 6) No. Instances AIR Closer
% 7) No. Instances CONV Closer
% 8) % Instances AIR Closer
num_vars = 8;
num_times = time_count;
dir_list = ["OA","AO"];
hour_list = start_hour:end_hour;

stats_thresh_w0 = zeros(num_directions,num_thresholds,num_vars).*NaN;
stats_thresh_w1 = zeros(num_directions,num_thresholds,num_vars).*NaN;
stats_time_w0 = zeros(num_directions,num_times,num_vars).*NaN;
stats_time_w1 = zeros(num_directions,num_times,num_vars).*NaN;

tdirs_a = string(cell2mat(T_A.direction));
tdirs_b = string(cell2mat(T_B.direction));

% Unweighted
for dir_idx = 1:num_directions
    for thresh_idx = 1:num_thresholds
        rows_a = T_A.threshold == thresh_list(thresh_idx) & tdirs_a == dir_list(dir_idx);
        rows_b = T_B.threshold == thresh_list(thresh_idx) & tdirs_b == dir_list(dir_idx);
        stats_thresh_w0(dir_idx,thresh_idx,1) = mean(T_A{rows_a,"Count"},"omitnan");
        stats_thresh_w0(dir_idx,thresh_idx,2) = mean(T_B{rows_b,"Count"},"omitnan");
        stats_thresh_w0(dir_idx,thresh_idx,3) = mean(T_A{rows_a,"Mean"},"omitnan");
        stats_thresh_w0(dir_idx,thresh_idx,4) = mean(T_B{rows_b,"Mean"},"omitnan");
        dif = T_A{rows_a,"Mean"}-T_B{rows_b,"Mean"};
        stats_thresh_w0(dir_idx,thresh_idx,5) = mean(dif,"omitnan");
        stats_thresh_w0(dir_idx,thresh_idx,6) = length(dif(dif < 0));
        stats_thresh_w0(dir_idx,thresh_idx,7) = length(dif(dif > 0));
        stats_thresh_w0(dir_idx,thresh_idx,8) = length(dif(dif < 0))/length(dif);
    end
    for time_idx = 1:num_times
        rows_a = T_A.timestamp == hour_list(time_idx) & tdirs_a == dir_list(dir_idx);
        rows_b = T_B.timestamp == hour_list(time_idx) & tdirs_b == dir_list(dir_idx);
        stats_time_w0(dir_idx,time_idx,1) = mean(T_A{rows_a,"Count"},"omitnan");
        stats_time_w0(dir_idx,time_idx,2) = mean(T_B{rows_b,"Count"},"omitnan");
        stats_time_w0(dir_idx,time_idx,3) = mean(T_A{rows_a,"Mean"},"omitnan");
        stats_time_w0(dir_idx,time_idx,4) = mean(T_B{rows_b,"Mean"},"omitnan");
        dif = T_A{rows_a,"Mean"}-T_B{rows_b,"Mean"};
        stats_time_w0(dir_idx,time_idx,5) = mean(dif,"omitnan");
        stats_time_w0(dir_idx,time_idx,6) = length(dif(dif < 0));
        stats_time_w0(dir_idx,time_idx,7) = length(dif(dif > 0));
        stats_time_w0(dir_idx,time_idx,8) = length(dif(dif < 0))/length(dif);
    end
end

% For weighted distances, the "sum" field is effectively the weighted mean,
% so use that.
tdirs_c = string(cell2mat(T_C.direction));
tdirs_d = string(cell2mat(T_D.direction));
for dir_idx = 1:num_directions
    for thresh_idx = 1:num_thresholds
        rows_c = T_C.threshold == thresh_list(thresh_idx) & tdirs_c == dir_list(dir_idx);
        rows_d = T_D.threshold == thresh_list(thresh_idx) & tdirs_d == dir_list(dir_idx);
        stats_thresh_w1(dir_idx,thresh_idx,1) = mean(T_C{rows_c,"Count"},"omitnan");
        stats_thresh_w1(dir_idx,thresh_idx,2) = mean(T_D{rows_d,"Count"},"omitnan");
        stats_thresh_w1(dir_idx,thresh_idx,3) = mean(T_C{rows_c,"Sum"},"omitnan");
        stats_thresh_w1(dir_idx,thresh_idx,4) = mean(T_D{rows_d,"Sum"},"omitnan");
        %% CHECK %%
        %% Should this be "Sum"?
        dif = T_C{rows_c,"Mean"}-T_D{rows_d,"Mean"};
        stats_thresh_w1(dir_idx,thresh_idx,5) = mean(dif,"omitnan");
        stats_thresh_w1(dir_idx,thresh_idx,6) = length(dif(dif < 0));
        stats_thresh_w1(dir_idx,thresh_idx,7) = length(dif(dif > 0));
        stats_thresh_w1(dir_idx,thresh_idx,8) = length(dif(dif < 0))/length(dif);
    end
    for time_idx = 1:num_times
        rows_c = T_C.timestamp == hour_list(time_idx) & tdirs_c == dir_list(dir_idx);
        rows_d = T_D.timestamp == hour_list(time_idx) & tdirs_d == dir_list(dir_idx);
        stats_time_w1(dir_idx,time_idx,1) = mean(T_C{rows_c,"Count"},"omitnan");
        stats_time_w1(dir_idx,time_idx,2) = mean(T_D{rows_d,"Count"},"omitnan");
        stats_time_w1(dir_idx,time_idx,3) = mean(T_C{rows_c,"Mean"},"omitnan");
        stats_time_w1(dir_idx,time_idx,4) = mean(T_D{rows_d,"Mean"},"omitnan");
        dif = T_C{rows_c,"Mean"}-T_D{rows_d,"Mean"};
        stats_time_w1(dir_idx,time_idx,5) = mean(dif,"omitnan");
        stats_time_w1(dir_idx,time_idx,6) = length(dif(dif < 0));
        stats_time_w1(dir_idx,time_idx,7) = length(dif(dif > 0));
        stats_time_w1(dir_idx,time_idx,8) = length(dif(dif < 0))/length(dif);
    end
end



%%




%% Save big stats table

% Unweighted


hour_strings = string(start_hour:end_hour);
dir_strings = ["OA","AO"];
thresh_strings = string(thresh_list);
%colnames = ["MNO-A","MNO-C","MMOD-A","MMOD-C","M-DifMOD","NIC-A","NIC-C","PC-A"];

T_thresh = combinations(dir_strings,thresh_strings);
T_time = combinations(dir_strings,hour_strings);

T = T_thresh;

num_rows = size(T,1);
T.MNO_A = zeros(num_rows,1);
T.MNO_C = zeros(num_rows,1);
T.MMOD_A = zeros(num_rows,1);
T.MMOD_C = zeros(num_rows,1);
T.M_difMOD = zeros(num_rows,1);
T.NIC_A = zeros(num_rows,1);
T.NIC_C = zeros(num_rows,1);
T.PC_A = zeros(num_rows,1);

% Variables:
% 1) AIR Ens Mean # Objects (MNO_A)
% 2) CONV Ens Mean # Objects (MNO_C)
% 3) AIR Ens Mean of Per-Member Mean Match Distance (MMOD_A)
% 4) CONV Ens Mean of Per-Member Mean Match Distance (MMOD_C)
% 5) Ens Mean Difference in Per-Member Mean Match Distances (M_difMOD)
% 6) No. Instances AIR Closer (NIC_A)
% 7) No. Instances CONV Closer (NIC_C)
% 8) % Instances AIR Closer (PC_A)
stats = stats_thresh_w0;
subname = "finalstats_thresh_w0";
row_idx = 1;
for dir_idx = 1:2
    for thresh_idx = 1:num_thresholds
        for col_idx = 1:8
            T(row_idx,col_idx+2) = num2cell(stats(dir_idx,thresh_idx,col_idx));
        end
        row_idx = row_idx+1;
    end
end
table_name = sprintf(output_format_table,output_path,data_tag,subname,domain,"XXX",data_type,"btXX",shape_string,distance_cutoff/1000,testing_value_string,dist_label,"20200207XXXX");
writetable(T,table_name,'Delimiter','\t','WriteRowNames',false);

T = T_thresh;
num_rows = size(T,1);
T.MNO_A = zeros(num_rows,1);
T.MNO_C = zeros(num_rows,1);
T.MMOD_A = zeros(num_rows,1);
T.MMOD_C = zeros(num_rows,1);
T.M_difMOD = zeros(num_rows,1);
T.NIC_A = zeros(num_rows,1);
T.NIC_C = zeros(num_rows,1);
T.PC_A = zeros(num_rows,1);
stats = stats_thresh_w1;
subname = "finalstats_thresh_w1";
row_idx = 1;
for dir_idx = 1:2
    for thresh_idx = 1:num_thresholds
        for col_idx = 1:8
            T(row_idx,col_idx+2) = num2cell(stats(dir_idx,thresh_idx,col_idx));
        end
        row_idx = row_idx+1;
    end
end
table_name = sprintf(output_format_table,output_path,data_tag,subname,domain,"XXX",data_type,"btXX",shape_string,distance_cutoff/1000,testing_value_string,dist_label,"20200207XXXX");
writetable(T,table_name,'Delimiter','\t','WriteRowNames',false);

T = T_time;
num_rows = size(T,1);
T.MNO_A = zeros(num_rows,1);
T.MNO_C = zeros(num_rows,1);
T.MMOD_A = zeros(num_rows,1);
T.MMOD_C = zeros(num_rows,1);
T.M_difMOD = zeros(num_rows,1);
T.NIC_A = zeros(num_rows,1);
T.NIC_C = zeros(num_rows,1);
T.PC_A = zeros(num_rows,1);
stats = stats_time_w0;
subname = "finalstats_time_w0";
row_idx = 1;
for dir_idx = 1:2
    for time_idx = 1:time_count
        for col_idx = 1:8
            T(row_idx,col_idx+2) = num2cell(stats(dir_idx,time_idx,col_idx));
        end
        row_idx = row_idx+1;
    end
end
table_name = sprintf(output_format_table,output_path,data_tag,subname,domain,"XXX",data_type,"btXX",shape_string,distance_cutoff/1000,testing_value_string,dist_label,"20200207XXXX");
writetable(T,table_name,'Delimiter','\t','WriteRowNames',false);

T = T_time;
num_rows = size(T,1);
T.MNO_A = zeros(num_rows,1);
T.MNO_C = zeros(num_rows,1);
T.MMOD_A = zeros(num_rows,1);
T.MMOD_C = zeros(num_rows,1);
T.M_difMOD = zeros(num_rows,1);
T.NIC_A = zeros(num_rows,1);
T.NIC_C = zeros(num_rows,1);
T.PC_A = zeros(num_rows,1);
stats = stats_time_w1;
subname = "finalstats_time_w1";
row_idx = 1;
for dir_idx = 1:2
    for time_idx = 1:time_count
        for col_idx = 1:8
            T(row_idx,col_idx+2) = num2cell(stats(dir_idx,time_idx,col_idx));
        end
        row_idx = row_idx+1;
    end
end
table_name = sprintf(output_format_table,output_path,data_tag,subname,domain,"XXX",data_type,"btXX",shape_string,distance_cutoff/1000,testing_value_string,dist_label,"20200207XXXX");
writetable(T,table_name,'Delimiter','\t','WriteRowNames',false);

       
%%
fprintf('%s run complete.\n',script_name);
toc
fprintf('DONE.\n')