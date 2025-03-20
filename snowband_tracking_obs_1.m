% snowband_tracking_obs.m
%
% Purpose: Identify precipitation objects, group into snowband candidates.
% This version uses observations as the Demons target, regridding them to
% WRF scale.
% 
% Author(s): Jon Seibert
% Last updated: 19 Feb 2025
% 
% Inputs: *.nc [Various NetCDF files]
% Outputs: cluster files, *.png
%
% Usage: Edit sections "0. General Settings" and  "1. Settings" to specify the inputs, run entire script
%
% TODO:
%   - Implement more snowband verification checks
%   - generalize file names and titles
%   - Implement demons imreg towards choice of target, not always obs
%  
% Dependencies: borders.m, reflmap.m, haversine_distance.m, fit_ ellipse.m,
% image processing toolbox, latlon_wrf_trim.mat, latlon_obs_trim.mat
% (Premade grid values)
%
% NOTES:
%   - Distance matrix becomes too large to run with grids exceeding
%   ~1000x1000 points. Overlarge grids must be downscaled.
%   - CURRENTLY ASSUMES THAT DEMONS IMREG MUST BE TOWARD CORRESPONDING REFL
%   OBSERVATIONS

script_name = 'snowband_tracking_obs.m';

%% 0A. Script Controls

% Denotes output folder name
run_name = 'CONV_base_3'; %AIR_SLC_4
exp_choice = 2;

use_base_refl = 1; % Whether to use vertical composite (0) or simulated base (1) reflectivity

reg_against = 0; % 0 = obs; if not obs, compares member to member, could be at a dif time or dif experiment [NOT IMPLEMENTED]
time_forward = 0; % How many timesteps forward to compare to [NOT IMPLEMENTED]

% Date & time
year = 2020;
month = 2;
day = 7;
hour = 14;
minute = 0;

num_members = 40;
use_member_list = 0; 
custom_member_list = 1:40;

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
point_req = 20; % Minimum number of points to consider a cluster for object status.
%nbh_radius = 10000; % Radius  of neighborhood circle [m] (6 km; current clustering cutoff radius is 10 km... try that? link the two?)
%nbh_req = 26; % Number of neighbors in radius required to be considered part of a shared band? % Will I be using this?
%dt_cutoff_value = 1.2; % Euclidean distance from edge threshold for merging SLINK clusters

% Sensistivity testing controls
thresh_list = [20,35,50];
%thresh_list = [50];
radius_list = [10000]; % List of nbh radii
%req_list = [20, 26, 32]; % List of nbh point count requirements [20 26 32]
req_list = [26]; % List of nbh point count requirements [20 26 32]
%distance_list = [14000]; % List of distance cutoffs
distance_list = [10000];

dtc_list = [2.4]; % Binary distance transform threshold [Unused, but needed]
%dtc_list = [1.2];

use_demons = true; % Make demons transforms and analyze them too

% Which plots to save? (ONLY MAKE RAW AND THRESH ONCE PER TIMESTAMP!)
make_raw = 1; % Raw reflectivity: usually vertical composite, to have something to work with.
make_thresh = 1;
make_basic = 0;
make_arcs = 0;
make_ellipse = 0;
make_boxes = 1;
make_dt = 0;
make_dt_states = 0;

%override_output = 0;
override_regrid = 0; % Remake wrf-obs regrid?
override_demons = 0; % Remake demons ImReg transform?
override_stats = 0; % Remake stats files?

% Print out processing notes?
show_progress = 1;
print_warnings = 0;
show_plots = 0;

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
    axes_font_size = 16;
    contour_font_size = 14;
else
    % Figure font sizes
    title_font_size = 24;
    label_font_size = 24;
    axes_font_size = 22;
    contour_font_size = 20;
end

% Lay out title and filename formats
% Input filenames, plot titles, output filenames, datetime strings
% Compressed = for small plots
filename_format = 'wrf_enkf_output_d%02.f_%03.f'; % domain, member #
filename_format_alt = 'wrfinput_d%02.f_%04.f-%02.f-%02.f_%02.f:00:00_%03.f'; % domain, year, month, day, hour, member #
%title_format = '[%s-d%02.f-%s] Neigh. Prob. of %s > %d %s, ROI %d km (%s)'; % Data source (caps), domain, member #, data type, cutoff value, units, ROI/1000, datetime
%title_format_compressed = 'NP>%d|%s|%s'; % cutoff, member #, datetime
output_format = '%s/%s_%s_d%02.f_%s_%s_t%s_%s_c%dkm_%s_%s_%s.png'; % output path, data source, plot type, domain, member #, data type, threshold, bounding shape, cutoff dist, testing values string, datetime, demons flag
            % for obs, domain  = d00, mem = 000; bounding shapes: bx, el, nb; 
            % if not using point radius and req, n0k_0; demons flag = d, nd (timestamp should be that of unmodified version)
output_format_stats = '%s/%s_%s_d%02.f_%s_%s_t%s_%s_c%dkm_%s_%s_%s.mat'; % output path, data source, plot type, domain, member #, data type, threshold, bounding shape, cutoff dist, testing values string, datetime, demons flag
output_format_raw = '%s/%s_%s_d%02.f_%s_%s_%s.png'; % output path, data source, plot type, domain, member #, data type, datetime
 

demons_format = '%s/demons_%s_d%02.f_%s_%s_%s_t%s_c%dkm_n%02.fk_%d_%s.mat'; % output path, data source, domain, member #, data type, VC/B, threshold, cutoff dist, nbh radius, nbh req, datetime
datetime_format_file = '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

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

addpath(path_to_extra_code); % Make sure Matlab can see the addon code

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
fig_x = 100;
fig_y = 100;
fig_width = 925;
fig_height = 900;
fig_width_small = 525;
fig_height_small = 500;

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

%% 0E. Obs Settings 
data_src_obs = "nex";
file_format_obs = '%s/n0q_%04.f%02.f%02.f%02.f%02.f.png';

lat_obs_raw = 49.9975:-0.005:23.0025;
lon_obs_raw = -125.9975:0.005:-65.0025;
[lon_obs,lat_obs] = meshgrid(lon_obs_raw,lat_obs_raw);

%% 1. Setup

if(~show_plots)
    set(groot,'DefaultFigureVisible','off') % Turn off figure popups for local
end

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

if(use_base_refl)
    vcb_string = "ba";
else
    vcb_string = "vc";
end

if(use_v6)
    slink_version = 6;
else
    slink_version = 5;
end

% Announce operating mode
fprintf('Starting %s in %s mode using SLINK V%d.\n',script_name,mode,slink_version);
tic;

% Directory structure:
% output/np/(run_name)/(variable_name)/(plot_size)
output_path_large = sprintf('%s/%s/%s/%s',output_path_base,run_name,data_type,'large');
output_path_small = sprintf('%s/%s/%s/%s',output_path_base,run_name,data_type,'small');

addpath(path_to_extra_code);

% If output folders do not exist, create them
if(~isfolder(output_path_large) || ~isfolder(output_path_small))
    mkdir(output_path_large);
    mkdir(output_path_small);
end

for band_threshold = thresh_list
for nbh_radius = radius_list
for nbh_req = req_list
for distance_cutoff = distance_list
for dt_cutoff_value = dtc_list
    
    if(use_v6)
        error("V6 removed");
    else
        testing_value_string = sprintf('%s_c%dkm_n%dkm_nr%01.f_pr%01.f',vcb_string,distance_cutoff/1000,nbh_radius/1000,nbh_req,point_req);
        fprintf('Initializing: T %ddBZ, C %dkm, NR %dkm, NPR %d, MinPR %d\n',band_threshold,distance_cutoff/1000,nbh_radius/1000,nbh_req,point_req);
    end 

%% 2. Initial Processing       

% Specify threshold strings
if(use_sd_threshold)
    thresh_string = sprintf('sd%01d-%02.f',floor(sd_mult),100*mod(sd_mult,1));
else
    thresh_string = sprintf('bt%02.f',band_threshold);
end

load(sprintf('%s/%s',intermediate_path,'latlon_wrf_trim.mat'),'n_idx','s_idx','e_idx','w_idx','lon_wrf_trim','lat_wrf_trim');
lon_wrf = lon_wrf_trim;
lat_wrf = lat_wrf_trim;
dimensions_wrf_trim = size(lon_wrf_trim);
ylen_wrf = dimensions_wrf_trim(1);
xlen_wrf = dimensions_wrf_trim(2);

load(sprintf('%s/%s',intermediate_path,'latlon_obs_trim.mat'),'n_idx_obs','s_idx_obs','e_idx_obs','w_idx_obs','lon_obs_trim','lat_obs_trim');
lon_obs = lon_obs_trim;
lat_obs = lat_obs_trim;
dimensions_obs = size(lon_obs);
ylen_obs = dimensions_obs(1);
xlen_obs = dimensions_obs(2);

if(regrid_obs_to_wrf)
    data_tag = sprintf('%s_%s',"OBS",exp_name);
else
    data_tag = sprintf('%s_%s',exp_name,"OBS");
end

%% 3. MAIN LOOP

if(use_member_list)
    member_list = custom_member_list;
else
    member_list = 1:num_members;
end

% For each member:
for member_idx = member_list
    
    member_string = sprintf('%03.f',member_idx);

    %% 3A. Read in data 
    
    % Define current loop's time strings
    timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
    timestamp_file = sprintf(datetime_format_file,year,month,day,hour);

    demons_string = 'd';
    stats_filename = sprintf(output_format_stats,intermediate_path,data_tag,'stats',domain,member_string,data_type,thresh_string,'nb',distance_cutoff/1000,testing_value_string,timestamp_file,demons_string);
    if(~override_stats && exist(stats_filename)) % If the output already exists, don't bother with this one
        fprintf('\nSKIPPING: Member %d, %s\n',member_idx,timestamp_title);
        continue;
    else
        fprintf('\nPROCESSING: Member %d, %s\n',member_idx,timestamp_title);
    end
    

    % Read in obs data
    % Convert from greyscale to dBZ values: dBZ = grey*0.5 - 32.5; (Still a bit unclear on WHY this is the conversion.)
    data = double(imread(sprintf(file_format_obs,input_path_obs,year,month,day,hour,minute)))*0.5 - 32.5; 
    data_obs = data(s_idx_obs:1:n_idx_obs,w_idx_obs:1:e_idx_obs); % Trim down to WRF domain

    % Specify input type
    if(exp_choice == 0) % Use OBS
        data_main = data_obs;
        lon_main = lon_obs;
        lat_main = lat_obs;
    else % Read in WRF data
        filename = sprintf(filename_format,domain,member_idx);
        data_wrf = ncread(sprintf('%s/%s/%s',input_path,timestamp_file,filename),data_name); % Retrieve data
        data_main = align_data(data_wrf,transpose_data,flip_data,z_composite_data,n_idx,s_idx,e_idx,w_idx);
        lon_main = lon_wrf;
        lat_main = lat_wrf;
    end

    clearvars data dimensions*;
    
    dims_trim = size(data_main);
    % Save original data
    %data_original = data;
    %data = data_main;

    %% Plot data prior to regrid

    % Also plot base refl for reference
        
    if(make_raw)
        % Plot LARGE
        f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
        h = pcolor(lon_main,lat_main,data_main); % Plot the data
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
        xlabel('Longitude (deg)','FontSize',label_font_size);
        ylabel('Latitude (deg)','FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);
        
        plot_type = 'raw';
        bounding_shape = 'nb';
        demons_string = 'nd';
    
        title_format_raw = '[%s|d0%d|%s] Raw: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
        title(sprintf(title_format_raw,exp_name,domain,member_string,data_type,units,timestamp_title),'FontSize',title_font_size);
        
        saveas(h,sprintf(output_format_raw,output_path_large,data_tag,plot_type,domain,member_string,data_type,timestamp_file)); % Save as .png
    end

    

    %% Regrid data to match obs, if needed



    if(exp_choice ~=0 && regrid_wrf_to_obs) % This will expand the WRF grid to match OBS
        data_regridded_filename = sprintf('%s/data_%s_%s_d0%d_%s_%s_%s_%s_%s.mat',intermediate_path,exp_name,"OBS",domain,member_string,"ao",upper(data_type),vcb_string,timestamp_file);
        if(~override_regrid && isfile(data_regridded_filename))
            if(show_progress)
                toc
                fprintf("Loading A->O WRF Regrid File [%s]...\n",data_regridded_filename);
            end
            load(data_regridded_filename,'data_wrf_regridded');
        else
            data_wrf_regridded = zeros(ylen_obs,xlen_obs); % Composite anyway, no need for vertical
            
            if(show_progress)
                toc
                fprintf('Regridding...\n'); % debug statement
            end
            for y = 1:ylen_obs
                for x = 1:xlen_obs
                    [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE] = find_nn_idx_irregular(lon_obs(y,x),lat_obs(y,x),lon_wrf,lat_wrf,edge_buffer); % Get nn index values
                    nn_points = [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE];        
                    if(any(nn_points == 0) || any(nn_points(1:4) > xlen_wrf) || any(nn_points(5:8) > ylen_wrf)) % If indices could not be found / reference points are NaN/ goes off edge
                        data_wrf_regridded(y,x) = nodata_val; % Lowest dBZ value used on color scale
                    else
                        data_wrf_regridded(y,x) = bilinear_interpolate_irregular_to_grid(y,x,lon_obs,lat_obs,data_main,lon_wrf,lat_wrf,nn_points);
                    end
                end
            end
            if(show_progress)
                toc
                fprintf("Saving A->O WRF Regrid File [%s]...\n",data_regridded_filename);
            end
            save(data_regridded_filename,"data_wrf_regridded");
        end
        data_main = data_wrf_regridded;
        lon_main = lon_obs;
        lat_main = lat_obs;
    elseif(regrid_obs_to_wrf) % This will condense the OBS grid to match WRF
        data_regridded_filename = sprintf('%s/data_%s_%s_d0%d_%s_%s_%s_%s_%s.mat',intermediate_path,"OBS",exp_name,domain,member_string,"oa",upper(data_type),vcb_string,timestamp_file);
        if(~override_regrid && isfile(data_regridded_filename))
            if(show_progress)
                toc
                fprintf("Loading O->A WRF Regrid File [%s]...\n",data_regridded_filename);
            end
            load(data_regridded_filename,'data_obs_regridded');
        else
            data_obs_regridded = zeros(ylen_wrf,xlen_wrf); % Composite anyway, no need for vertical
            
            if(show_progress)
                toc
                fprintf('Regridding...\n'); % debug statement
            end
            for y = 1:ylen_wrf
                for x = 1:xlen_wrf
                    [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE] = find_nn_idx_irregular(lon_wrf(y,x),lat_wrf(y,x),lon_obs,lat_obs,edge_buffer); % Get nn index values
                    nn_points = [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE];        
                    if(any(nn_points == 0) || any(nn_points(1:4) > xlen_obs) || any(nn_points(5:8) > ylen_obs)) % If indices could not be found / reference points are NaN/ goes off edge
                        data_obs_regridded(y,x) = nodata_val; % Lowest dBZ value used on color scale
                    else
                        data_obs_regridded(y,x) = bilinear_interpolate_irregular_to_grid(y,x,lon_wrf,lat_wrf,data_obs,lon_obs,lat_obs,nn_points);
                    end
                end
            end
            if(show_progress)
                toc
                fprintf("Saving O->A WRF Regrid File [%s]...\n",data_regridded_filename);
            end
            save(data_regridded_filename,"data_obs_regridded");
        end
        data_obs = data_obs_regridded;
        if(exp_choice == 0)
            data_main = data_obs_regridded;
            lon_main = lon_wrf;
            lat_main = lat_wrf;
        end
    end


    dims_trim = size(data_main);
 
    %% 3B. Threshold and set up data
    
    % Base threshold: cut off below [0] dBZ
    flat_list_obj = find(data_main>object_threshold);
    data_obj = data_main(flat_list_obj);
    num_points_obj = length(flat_list_obj);
    
    sd = std(data_obj); % Compute standard deviation
    mn = mean(data_obj); % Compute mean
    
    if(use_sd_threshold) % Compute new band threshold dynamically
        band_threshold = mn + (sd*sd_mult); 
    end % Otherwise uses preset value in settings
    
    data_a = data_main;
    data_a(data_a<=band_threshold) = 0;
    
    %% DEMONS Algorithm Loop
    for demons_idx = 1:2 

        if(show_progress)
            if(demons_idx == 1)
                fprintf('Processing: Original data form\n');
            else
                fprintf('Processing: Demons data form\n');
            end
        end
        
        if(demons_idx == 2)
            if(~use_demons) % If not using Demons this time, skip demons half of the loop
                if(show_progress)
                    fprintf('Skipping Demons.\n');
                end
                continue;
            else
                demons_filename = sprintf(demons_format,intermediate_path,data_tag,domain,member_string,data_type,vcb_string,thresh_string,distance_cutoff/1000,nbh_radius/1000,nbh_req,timestamp_file);
                
                % If the demons transformed version of data_thresh already exists,
                % load it. If not, create it.
                if(~override_demons && (exist(demons_filename,'file') == 2))
                    if(show_progress)
                        toc;
                        fprintf('Loading Demons ImReg File [%s]...\n',demons_filename);
                    end
                    load(demons_filename,'D','a_trans');
                elseif(use_demons)
                    if(show_progress)
                        toc;
                        fprintf('Performing Demons Image Registration...\n')
                    end
                    % Read in obs data
                    % Convert from greyscale to dBZ values: dBZ = grey*0.5 - 32.5; (Still a bit unclear on WHY this is the conversion.)
                    data_b = double(imread(sprintf(file_format_obs,input_path_obs,year,month,day,hour,minute)))*0.5 - 32.5; 
                    data_b = data_b(s_idx_obs:1:n_idx_obs,w_idx_obs:1:e_idx_obs); % Trim down to WRF domain
        
                    data_b(data_b<=(band_threshold+band_threshold_b_adjustment)) = nodata_val;

                    if(regrid_obs_to_wrf)
                        data_regridded_filename = sprintf('%s/data_%s_%s_d0%d_%s_%s_%s_%s_%s.mat',intermediate_path,"OBS",exp_name,domain,member_string,"oa",upper(data_type),vcb_string,timestamp_file);
                        if(~override_regrid && isfile(data_regridded_filename))
                            if(show_progress)
                                toc
                                fprintf("Loading O->A WRF Regrid File [%s]...\n",data_regridded_filename);
                            end
                            data_obs_regridded_temp = data_obs_regridded;
                            load(data_regridded_filename,'data_obs_regridded');
                            data_b_regridded = data_obs_regridded;
                            data_obs_regridded = data_obs_regridded_temp;
                        else
                            data_b_regridded = zeros(ylen_wrf,xlen_wrf); % Composite anyway, no need for vertical
                            for y = 1:ylen_wrf
                                for x = 1:xlen_wrf
                                    [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE] = find_nn_idx_irregular(lon_wrf(y,x),lat_wrf(y,x),lon_obs,lat_obs,edge_buffer); % Get nn index values
                                    nn_points = [xNW,xNE,xSW,xSE,yNW,yNE,ySW,ySE];        
                                    if(any(nn_points == 0) || any(nn_points(1:4) > xlen_obs) || any(nn_points(5:8) > ylen_obs)) % If indices could not be found / reference points are NaN/ goes off edge
                                        data_b_regridded(y,x) = nodata_val; % Lowest dBZ value used on color scale
                                    else
                                        data_b_regridded(y,x) = bilinear_interpolate_irregular_to_grid(y,x,lon_wrf,lat_wrf,data_b,lon_obs,lat_obs,nn_points);
                                    end
                                end
                            end
                        end
                        data_b = data_b_regridded;
                    end
        
                    % Temporarily replace NaNs with filler value to avoid code
                    % errors
                    if(clean_demons)
                        data_b_clean = data_b;
                        data_b_clean(isnan(data_b)) = 0;
                        data_a_clean = data_a;
                        data_a_clean(isnan(data_a)) = 0;
                    end
                    
                    % Perform image registration (phase correlation & demons)
                    tformEstimate = imregcorr(data_a_clean,data_b_clean,"translation");
                    Rfixed = imref2d(size(data_b_clean));
                    a_corr = imwarp(data_a_clean,tformEstimate,"OutputView",Rfixed);
                    [D,a_trans] = imregdemons(a_corr,data_b_clean); % D(:,:,1) contains x-axis displacements, (...2) contains y-axis displacements 
                    if(show_progress)
                        toc;
                        fprintf('Saving Demons ImReg File [%s]...\n',demons_filename);
                    end
                    save(demons_filename,'D','a_trans'); % Save as .png
                    clearvars data_a*;
                end
            end
        end
        
        % Set flags
        if(demons_idx == 1)
            demons_flag = false;
            demons_string = 'nd';
        else
            data_main = a_trans;
            demons_flag = true;
            demons_string = 'd';
            if(regrid_wrf_to_obs)
                lon_main = lon_obs;
                lat_main = lat_obs;
            elseif(regrid_obs_to_wrf)
                lon_main = lon_wrf;
                lat_main = lat_wrf;
            end
        end
        
        % Get indices of above-threshold points
        %[row_list,col_list] = find(data_main>band_threshold); % FIND A WAY TO CUT OFF THE SE CORNER REFLECTIVITY??
        flat_list = find(data_main>band_threshold);
        num_points = length(flat_list); 
        
        % Extract above-threshold points only
        data_thin = data_main(flat_list);
        lon_thin = lon_main(flat_list);
        lat_thin = lat_main(flat_list);

        % Plot thresholded version
        if(make_thresh)
            data_thresh = data_main;
            data_thresh(:) = nodata_val;
            data_thresh(flat_list) = data_main(flat_list);

            % Plot LARGE
            f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
            h = pcolor(lon_main,lat_main,data_thresh); % Plot the data
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
            xlabel('Longitude (deg)','FontSize',label_font_size);
            ylabel('Latitude (deg)','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);
            
            plot_type = 'thresh';
            bounding_shape = 'nb';
        
            title_format_raw = '[%s|d0%d|%s] Thresholded: %s(%s)|%s'; % exp name, domain, member #/mean, data type, units, timestamp
            title(sprintf(title_format_raw,exp_name,domain,member_string,data_type,units,timestamp_title),'FontSize',title_font_size);
            
            saveas(h,sprintf(output_format,output_path_large,data_tag,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png
        end
        
        % DISTANCE TRANSFORM
        %mask = single(zeros(dims_trim(1),dims_trim(2)));
        
        %mask(flat_list) = 1;
        %dist_trans = bwdist(~mask);
        %dist_trans_thin = dist_trans(flat_list);
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
    
        if(use_v6)
            error("V6 removed");
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
        %dist_matrix = createArray(num_points,num_points,"int8");
        dist_matrix(:,:) = NaN; % n-n distance is 0, but we don't want the algorithm to find "0" as the min distance
        
        nn_idx_list = zeros(num_points,1); % Tracker for nearest neighbor of each point/cluster
        nn_dist_list = zeros(num_points,1); % Tracker for distance between nearest neighbors
        
        %nbh_radius = 20000; % Radius  of neighborhood circle [m] (6 km; current clustering cutoff radius is 10 km... try that? link the two?)
        %nbh_req = 40; % Number of neighbors in radius required to be considered part of a shared band? % Will I be using this?
        nbh_counts = zeros(num_points,1); % Counting variable
        
        % GENERATE DISTANCE MATRIX
        if(show_progress)
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
            %if(show_progress)
            %    if(idx_a == floor(num_points)/4)
            %        fprintf('25%% complete.\n')
            %    elseif(idx_a == floor(num_points/2))
            %        fprintf('50%% complete.\n')
            %    elseif(idx_a == floor(num_points*(3/4)))
            %        fprintf('75%% complete.\n')
            %    end
            %end
        end
        
        if(show_progress)
            fprintf('Done.\n')
            toc
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
            
            num_clusters = num_clusters - 1; % Decrement counter- there will always be one fewer independent "cluster"

            if(num_clusters <= 3)
                5;
            end
            
            if(show_progress)
                if(num_clusters == floor(num_points/4))
                    fprintf('75%% complete.\n')
                elseif(num_clusters == floor(num_points/2))
                    fprintf('50%% complete.\n')
                elseif(num_clusters == floor(num_points*(3/4)))
                    fprintf('25%% complete.\n')
                end
            end
        end
        clearvars dist_matrix;
        if(show_progress) 
            fprintf('Done.\n');
            toc 
            fprintf('Analyzing snowband criteria...\n');
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
        xlabel('Longitude (deg)','FontSize',label_font_size);
        ylabel('Latitude (deg)','FontSize',label_font_size);
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

        is_snowband = (length_flags & ratio_flags);  % WIP
        
        % Save tracking stats for snowband_dif
        stats_filename = sprintf(output_format_stats,intermediate_path,data_tag,'stats',domain,member_string,data_type,thresh_string,'nb',distance_cutoff/1000,testing_value_string,timestamp_file,demons_string);
        if(demons_flag)
            D_save = D;
        else
            D_save = NaN;
        end

        nonzero_idx = find(centroids(:,1) ~= 0);
        centroids_trim = centroids(nonzero_idx,:);

        if(override_stats || exist(stats_filename,'file') ~=2)
            save(stats_filename,'pt_clusters','cluster_ids','centroids_trim','data_main','D_save');
        end
         
        %% 4A. Plotting - Cluster fits

        % Initial figure made above, before cluster math loop

        % Plot clusters with fit lines
        if(make_arcs)
        
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
        
            saveas(h,sprintf(output_format,output_path_large,data_tag,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png
                    
            % output_format = '%s/%s_%s_d%02.f_%s_%s_t%s_%s_c%dkm_n%02.fk_%d_%s_%s.png'; % output path, data source, plot type, domain, member #, data type, threshold, bounding shape, cutoff dist, nbh radius, nbh req, datetime, demons flag
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
            saveas(h,sprintf(output_format,output_path_large,data_tag,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png
        
            hold off;
            close("all");
        end
        

        
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
        xlabel('Longitude (deg)','FontSize',label_font_size);
        ylabel('Latitude (deg)','FontSize',label_font_size);
        set(gca,'Fontsize',axes_font_size);
    
        if(make_basic)
            saveas(h,sprintf(output_format,output_path_large,data_tag,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png
        end

        if(make_ellipse)
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
            xlabel('Longitude (deg)','FontSize',label_font_size);
            ylabel('Latitude (deg)','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);
                
            saveas(h,sprintf(output_format,output_path_large,data_tag,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png
        end
        close("all");
        
        % Replot base clusters with boxes
        if(make_boxes)
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
                h.MarkerEdgeColor = "#ffffff";
                h.MarkerFaceColor = "#000000";
            end
        
            if(label_clusters)
                labelpoints(centroids_trim(:,1),centroids_trim(:,2),[string(1:size(centroids_trim,1))],'FontSize',axes_font_size-5,'FontWeight','bold');
            end
            
            bounding_shape = 'bx';
            
            title(sprintf('[%s|d%02.f|%s] Refl Band Clusters (Box): %s',upper(data_src),domain,upper(member_string),timestamp_title),'FontSize',title_font_size);
            xlabel('Longitude (deg)','FontSize',label_font_size);
            ylabel('Latitude (deg)','FontSize',label_font_size);
            set(gca,'Fontsize',axes_font_size);
            saveas(h,sprintf(output_format,output_path_large,data_tag,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png
        
            hold off;
            close("all");
        end
    
        % Plot SMALL
        %f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
        %title(sprintf(title_format_compressed,cutoff,member_string,timestamp_title),'FontSize',title_font_size);
        %saveas(h,sprintf(output_format,output_path_small,data_src,domain,member_string,data_type,cutoff,units,roi/1000,timestamp_file)); % Save as .png
        %close('all');
        
        %% DISTANCE TRANSFORM
        
        if(make_dt || make_dt_states)
            mask = single(zeros(dims_trim(1),dims_trim(2)));
            mask(flat_list) = 1;
            %imshow(mask);
            dist_trans = bwdist(~mask);
            %D_norm = D./(max(D,[],"all"));
        
            f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
            h = pcolor(lon_main,lat_main,dist_trans);
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
            xlabel('Longitude (deg)','FontSize',label_font_size);
            ylabel('Latitude (deg)','FontSize',label_font_size);
            ylabel(c,'Normalized Euclidean Distance to Reflectivity Object Edge','FontSize',label_font_size)
            set(gca,'Fontsize',axes_font_size);
            
            plot_type = 'dist_trans_nb';
            bounding_shape = 'nb';
            
            if(make_dt)
                saveas(h,sprintf(output_format,output_path_large,data_tag,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png
            end

            hold on;
            borders('continental us','color','#20891e','linewidth',1); 
            hold off;
            
            plot_type = 'dist_trans';
            if(make_dt_states)
                saveas(h,sprintf(output_format,output_path_large,data_tag,plot_type,domain,member_string,data_type,thresh_string,bounding_shape,distance_cutoff/1000,testing_value_string,timestamp_file,demons_string)); % Save as .png
            end
        end

        if(show_progress) 
            fprintf('Done.\n');
            toc 
        end
    
    end % End demons loop
end

if(show_progress)
    fprintf("All members completed. Moving to next iteration...\n");
end

end
end
end
end
end

set(groot,'DefaultFigureVisible','on');
% Get total runtime and print to stdout
runtime = toc;
hours = floor(runtime/3600);
mins = floor((runtime/60) - (hours*60));
secs = toc - (hours*3600) - (mins*60);
fprintf('%s run complete.\n',script_name);
fprintf('Total script runtime = %02.f:%02.f:%02.f\n',hours,mins,secs);
fprintf('#---------------------------------------------------------------------#\n');

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
