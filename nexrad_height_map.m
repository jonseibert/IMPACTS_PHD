% nexrad_height_map (v1.0.0)
%
% Purpose: Generate 2D map of the altitudes NEXRAD is unable to see (lowest
% visible)
% 
% Author(s): Jon Seibert
% Last updated: 15 Sept 2024
% 
% Inputs: NEXRAD radar station info
% Outputs: [Height map].png
% Dependencies: borders.m, linear_interpolate.m, 
%       latlon_[]_trim.mat (predefined domain boundaries)
%
% NOTES:
%   - 
%
% TODO:
%  - 
script_name = 'nexrad_height_map.m';
   
%% 0A. Script Controls

run_name = 'vslice';

% Local or server execution mode
%local = 1; % true = running on local device, false = running on server

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

num_members = 1; % TESTING

domain = 2;
subdomain = 0; %1 for E Coast Only subdomain, 2 for N only, 0 for all d02 included

limit_raw_colorbar = true;
clim_lower = 0; % Colorbar limits
clim_upper = 10000;
clim_lower_mb = 200;
clim_upper_mb = 1000;


% Plot selections
remake_plots = 1; % If true, overrides existing outputs. If false, only plots 'new' files.

show_plots = 0;
show_progress = 1; % If true, prints out progress bar/debug messages
announce_plots = 0;

gauss_sd = 2.75;
gauss_width =  2*ceil(2*gauss_sd)+1;


radar_heights_filename = 'NEXRAD_radar_heights_wrf_grid_trim.mat';
wrf_heights_filename = 'wrf_base_refl_heights.mat';

max_dist = 460000; %(m): radar radius

%% 0B. General Settings

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
end

addpath(path_to_extra_code);

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

% Subdomain settings
bound_lon_d1 = [-78,-78,-74.5,e_lim,e_lim,-75,-78];
bound_lat_d1 = [40,44,n_lim,n_lim,41.5,38,40];
bound_lon_d2 = [-77.5,-77.5,-74,e_lim,e_lim,-73.5,-77.5];
bound_lat_d2 = [42,44,n_lim,n_lim,42.5,42,42];


Re = 6.3781e6; % Radius pf Earth in m
g = 9.80665; % m/s^2 standard gravity at sea level
R = 287.0600676; % J/kg/K, specific gas constant

ecmwf_filename = 'ecmwf_mslp_d01_20200207.nc';
       
beam_width = 0.9; % Degrees, DIAMETER
beam_angle = 0.5;

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

datetime_format_file =  '%04.f%02.f%02.f%02.f00'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f00'; % year, month, day, hour [as numbers]

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

dimensions_wrf = size(lon_wrf);
dims = size(lon_wrf);
ylen_wrf = dimensions_wrf(1);
xlen_wrf = dimensions_wrf(2);
zlen_wrf = 50;
clearvars dimensions_wrf lon_wrf_trim lat_wrf_trim;



%% MAIN


% Set timestamps and experiment names
timestamp_title = sprintf(datetime_format_title,year,month,day,hour);
timestamp_file = sprintf(datetime_format_file,year,month,day,hour);
timestamp_file_background = sprintf(datetime_format_file,year,month,day,hour-1);
data_src = 'nex';


%% Plot vertical slices: x and y-aligned only

% Existing model height field name is "model_height" (above ground)
% Load NEXRAD beam heights and WRF Base Refl height map
load(sprintf('%s/%s',intermediate_path,radar_heights_filename),'radar_height_grid');
radar_height_grid_full = radar_height_grid;
load(sprintf('%s/%s',intermediate_path,'nexrad_stations.mat')); % wban, station_ids,_station_names,lon_deg,lat_deg,elevations,tower_heights
num_radars = length(station_ids);

min_heights = zeros(ylen_wrf,xlen_wrf);
min_pressures = zeros(ylen_wrf,xlen_wrf);

for y = 1:ylen_wrf
    for x = 1:xlen_wrf
        min_heights(y,x) = min(radar_height_grid(:,y,x));
    end
end


load(sprintf('%s/%s',intermediate_path,'mean_storage_AIR_202002071500.mat'),'data_p','model_height');
pressures = data_p;

for y = 1:ylen_wrf
    for x = 1:xlen_wrf
        for z = 1:zlen_wrf
            if(model_height(y,x,z) > min_heights(y,x))
                min_pressures(y,x) = linear_interpolate(model_height(y,x,max(1,z-1)),model_height(y,x,z),min_heights(y,x),pressures(y,x,(max(1,z-1))),pressures(y,x,z));
                break;
            end
        end
    end
end

load(sprintf('%s/%s',intermediate_path,'nexrad_stations.mat')); % wban, station_ids,_station_names,lon_deg,lat_deg,elevations,tower_heights

% Find the radars that are in the domain, plot only those labels
num_radars_all = size(station_ids,1);
num_radars = 0;
radar_checklist = zeros(num_radars_all,1);
for rad_idx = 1:num_radars_all
    if(lon_deg(rad_idx) > w_lim && lon_deg(rad_idx) < e_lim && lat_deg(rad_idx) > s_lim && lat_deg(rad_idx) < n_lim)
        num_radars = num_radars+1;
        radar_checklist(rad_idx) = 1;
    end
end

radars_used = find(radar_checklist == 1);

elevations_rounded = string(zeros(size(radars_used)));
for rad_idx = 1:num_radars
    rad_cur = radars_used(rad_idx);
    elevations_rounded(rad_idx) = sprintf("%3.1f",(elevations(rad_cur)+tower_heights(rad_cur)));
end


%% Plot vertical slices

% Lay out title and filename formats
title_format =  '[NEX|d0%d] Minimum Radar Return Altitude (m above ground level)'; % domain
title_format_mb = '[AIR|d0%d|%s] Minimum Radar Return Pressure (mb)|%s'; % domain, member #/mean, timestamp
output_format = '%s/nex_air_d0%d_%s_%s_%s.png'; % output path, domain, mem/mean, plot_type, timestamp

beam_color = "#000000";

plot_type = 'mha';
data_to_plot = min_heights;
cmap = 'jet';

member_string = "mean";

plot_filename = sprintf(output_format,output_path_large,domain,member_string,plot_type,timestamp_file);

if(isfile(plot_filename) && ~remake_plots) % If shouldn't override existing plots, don't
    %continue;
end

% Plot LARGE
f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure

if(show_progress && announce_plots)
    toc
    fprintf('Plotting %s...\n',plot_type); % debug statement
end

h = pcolor(lon_wrf,lat_wrf,data_to_plot); % Plot the data

set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
shading interp; % Smooth out plot from grid-boxes
c = colorbar('FontSize',axes_font_size); % Make colorbar
colormap(flipud(jet)); % Set colors

hold on;
% Plot state borders
borders('continental us','black','linewidth',1); 
hold off;

xlim([w_lim e_lim]);
ylim([s_lim n_lim]);

% Apply labels
title(sprintf(title_format,domain),'FontSize',title_font_size); 
xlabel('Longitude','FontSize',label_font_size);
ylabel('Latitude','FontSize',label_font_size);
set(gca,'Fontsize',axes_font_size);
%set (c,'YDir','reverse')

if(limit_raw_colorbar)
    caxis([clim_lower clim_upper]);
end

saveas(gcf,plot_filename); % Save as .png



%%

plot_type = 'mhp';
data_to_plot = min_pressures;
cmap = 'jet';

plot_filename = sprintf(output_format,output_path_large,domain,member_string,plot_type,timestamp_file);

if(isfile(plot_filename) && ~remake_plots) % If shouldn't override existing plots, don't
    %continue;
end

% Plot LARGE
f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure

if(show_progress && announce_plots)
    toc
    fprintf('Plotting %s...\n',plot_type); % debug statement
end

h = pcolor(lon_wrf,lat_wrf,data_to_plot); % Plot the data

set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
shading interp; % Smooth out plot from grid-boxes
c = colorbar('FontSize',axes_font_size); % Make colorbar
colormap(cmap); % Set colors

hold on;
% Plot state borders
borders('continental us','black','linewidth',1); 




xlim([w_lim e_lim]);
ylim([s_lim n_lim]);

% Apply labels
title(sprintf(title_format_mb,domain,member_string,timestamp_title),'FontSize',title_font_size); 
xlabel('Longitude','FontSize',label_font_size);
ylabel('Latitude','FontSize',label_font_size);
set(gca,'Fontsize',axes_font_size);

if(limit_raw_colorbar)
    caxis([clim_lower_mb clim_upper_mb]);
end

saveas(gcf,plot_filename); % Save as .png



% Add labels and save again

h = scatter(lon_deg,lat_deg,100,"white","filled","pentagram");
h.MarkerFaceColor = [1 1 1];
h.MarkerEdgeColor = [0 0 0];
labelpoints(lon_deg(radars_used),lat_deg(radars_used),elevations_rounded+" m",'position','N','FontSize',axes_font_size-5,'FontWeight','bold','Color','white');

hold off;
plot_type = 'mhp_label';
data_to_plot = min_pressures;
cmap = 'jet';

plot_filename = sprintf(output_format,output_path_large,domain,member_string,plot_type,timestamp_file);
saveas(gcf,plot_filename); % Save as .png

clearvars f h data_to_plot;
    
%%

% Get total runtime and print to stdout
runtime = toc;
hours = floor(runtime/3600);
mins = floor((runtime/60) - (hours*60));
secs = toc - (hours*3600) - (mins*60);
fprintf('Done. Total script runtime = %02.f:%02.f:%02.f\n',hours,mins,secs)

% END

