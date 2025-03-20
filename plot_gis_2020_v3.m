% ptp_analysis.m
% Purpose: Generate point-to-point (PTP) diagnostic comparison plots of outputs
% from IMPACTS WRF to EXRAD and NEXRAD observations
% 
% Author(s): Jon Seibert
% Last updated: 26 January 2023
% 
% Inputs: Various
% Outputs: *.png [ ]
%
% Process: (For each comparable timestep)
%   1) Read in obs, background, and analysis steps for a given datetime
%   2) Calculate basic ensemble mean, use mean for comparisons
%   3) Regrid WRF to each obs for direct comparison
%   4) Compare background-analysis, background-obs and analysis-obs and dif between
%   5) Calculate RMSE and Bias (mean of dif) for each
%   6) Plot all
%
% TODO:
%   - ADD BACKGROUND INPUT TO COMPARE IT TO ANALYSIS AND/OR OBS (previous
%   timestep) Format: wrfinput_d02_2022-01-29_13:00:00_001, in 1200 folder,
%   for 1300
%
% NOTES:
%  - TIME_COUNT is not robust when crossing month boundaries- replace with
%    manual time-tick count if using it to do so.
%  - NOT ROBUST TO USE IN THE SOUTHERN HEMISPHERE LATITUDES!

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
    input_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/';
    input_path_gis = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/gis_data/2020';
    intermediate_path = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/intermediate';
    output_path_base = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/';
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

%% 0. General Settings

% Date & time
% Range: 2020-02-07-1400 : 2020-02-07-1800
year = 2020;
month = 2;
day = 7;
%ex_hour_list = ["1412","1443","1501","1516"];
%gis_hour_list = ["1410","1445","1500","1515"];
%gis_hour_list = ["1200","1300","1400","1500","1600","1700","1800","1900","2000","2100","2200","2300"];
gis_hour_list = ["0000","0100","0200","0300","0400","0500","0600","0700","0800","0900","1000","1100","1200","1300","1400","1500","1600","1700","1800","1900","2000","2100","2200","2300"];
%gis_hour_list = ["0000","0100","0200","0300","0400"];

% Figure specs
fig_x = 100;
fig_y = 100;
fig_width = 925;
fig_height = 900;
%fig_width_small = 350;
%fig_height_small = 350;
fig_width_small = 525;
fig_height_small = 500;

limit_borders = true;

w_lim = -79;
e_lim = -69.75;
s_lim = 36;
n_lim = 46;

%w_lim = -115; % Degrees lat/lon
%e_lim = -68;
%s_lim = 24;
%n_lim = 48;

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

data_type = "refl"; % Variable being analyzed
units = "dBZ";

cmap = 'reflmap'; % Choice of colorbar colormap
clim_lower = -30; % Colorbar limits
clim_upper = 75;

dbz_nodata = -30;
dbz_nodata_gis = -66;

plot_boxes = false;

% Progress messages
show_plots = 0; % 0: Suppress plot popups [REQUIRES MATLAB RESTART TO UNDO]
show_progress = 1; % If true, prints out progress bar/debug messages
announce_plots = 0; % If true, prints out individual plot messages

%% 1B. Settings- GIS Data (NEXRAD, but in images)

%input_path_gis = strcat(input_path_gis,'/gis_data/2020');

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

% Lay out title and filename formats
datetime_format_file = '%04.f%02.f%02.f%02.f%02.f'; % year, month, day, hour [as numbers]
datetime_format_title = '%04.f-%02.f-%02.f-%02.f%02.f'; % year, month, day, hour [as numbers]

time_count = length(gis_hour_list);

load(sprintf('%s/%s',intermediate_path,'latlon_gis_trim.mat'));


%% For each timestamp being analyzed:
for time_idx = 1:time_count
    
    gis_hour = double(extractBetween(gis_hour_list(time_idx),1,2));
    gis_minute = double(extractBetween(gis_hour_list(time_idx),3,4));
   
    timestamp_title = sprintf(datetime_format_title,year,month,day,gis_hour,gis_minute);
    timestamp_file = sprintf(datetime_format_file,year,month,day,gis_hour,gis_minute);

    % Read in GIS data
    % Convert from greyscale to dBZ values: dBZ = grey*0.5 - 32.5; (Still a bit unclear on WHY this is the conversion.)
    data_gis = double(imread(sprintf(file_format_gis,input_path_gis,year,month,day,gis_hour,gis_minute)))*0.5 - 32.5; 
    %data_gis_trim = data_gis(s_idx:1:n_idx,w_idx:1:e_idx); % Trim down to WRF domain
    %dimensions_gis_trim = size(data_gis_trim);
    data_gis_trim = data_gis; % Trim down to WRF domain
    dimensions_gis_trim = size(data_gis_trim);
    lon_gis_trim = lon_gis;
    lat_gis_trim = lat_gis;
   

    data_to_plot = data_gis_trim;
    cmap = 'reflmap';

    % Plot LARGE
    f = figure('Position',[fig_x fig_y fig_width fig_height]); % Create initial blank figure
    h = pcolor(lon_gis_trim,lat_gis_trim,data_to_plot); % Plot the data
    set(h, 'EdgeColor', 'none'); % Remove weird grid-boxes from pcolor
    shading interp; % Smooth out plot from grid-boxes
    c = colorbar('FontSize',axes_font_size); % Make colorbar
    colormap(cmap); % Set colors
    caxis([clim_lower clim_upper]);

    % Plot state borders
    hold on;
    borders('continental us','black','linewidth',1); 
    
    % Focus on desired area, remove whitespace
    %if(limit_borders)
        xlim([w_lim e_lim]);
        ylim([s_lim n_lim]);
    %end

    
    hold off;
    % Apply labels
    title(sprintf('[NEXRAD] 0.5%s Base Reflectivity Mosaic (dBZ)|%s',char(176),timestamp_title),'FontSize',title_font_size);
    xlabel('Longitude','FontSize',label_font_size);
    ylabel('Latitude','FontSize',label_font_size);
    set(gca,'Fontsize',axes_font_size);
    saveas(h,sprintf('%s/nex/nex_trim_v2_%s.png',output_path_base,timestamp_file)); % Save as .png

    f.Position = [fig_x fig_y fig_width_small fig_height_small]; % Shrink figure
    saveas(h,sprintf('%s/nex/small/nex_trim_v2_%s_small.png',output_path_base,timestamp_file)); % Save as .png
    
    close('all');
    
end